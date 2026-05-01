library(tidyverse)
library(arrow)
library(data.table)
library(zoo)
library(fst)
library(survival)
library(gnm)
library(ggplot2)
library(broom)
library(sf)
library(purrr)
library(doParallel)
library(progress)


aggregate_data <- readRDS("/n/dominici_nsaph_l3/Lab/projects/garambyun_pm25-longterm-mortality-urbanrural/Data/aggregate_data.rds") %>%
  mutate(across(c(year, region, sex, race, dual, entry_age_break, followup_year), as.factor))
exposure <- readRDS("/n/dominici_nsaph_l3/Lab/projects/garambyun_pm25-longterm-mortality-urbanrural/Data/exposure.rds")

##########################################################
#################### Loop - components ###################
##########################################################

##### All area #####
exp_name <- c("pm25_ensemble","oc","ec","so4","no3","nh4",
              "br","ca","cu","fe","k","ni","pb","si","v","z")
              
aggregate_data.list <- split(aggregate_data, aggregate_data$zip)
num_uniq_zip <- length(aggregate_data.list)
registerDoParallel(cores = 8)

n_per_batch <- 50
n_batch <- 10


start_time <- Sys.time()

for (i in seq_along(exp_name)) {
  tryCatch({
    mexp <- exp_name[i]
    mcov <- paste0("PM25_", mexp)
    iqr <- quantile(get(mexp, exposure), 0.75, na.rm = TRUE) - quantile(get(mexp, exposure), 0.25, na.rm = TRUE)
    
    result_df_list <- list()
    
    
    for (model_type in 1:2) {
      
      # 1. Full model
      formula_base <- if (model_type == 1) paste0("dead ~ ", mexp) else paste0("dead ~ ", mexp, " + ", mcov)
      
      formula_full <- as.formula(paste0(formula_base,
                                        " + mean_bmi + smoke_rate + hispanic + pct_blk +",
                                        " medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ +",
                                        " summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +",
                                        " year + region + offset(log(time_count))"))
      
      full_model <- tryCatch({
        gnm(formula_full,
                        eliminate = (sex:race:dual:entry_age_break:followup_year),
                        data = aggregate_data,
                        family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (!inherits(full_model, "gnm")) next
      
      full_coef <- coef(full_model)[1]
      
      # 2. Bootstrap
      boot_coefs_all <- c()
      
      total_iter <- n_per_batch * n_batch
      pb <- txtProgressBar(min = 0, max = total_iter, style = 3)
      iter_count <- 0 
      
      for (j in 1:n_batch) {
        boot_coefs_part <- foreach(b = 1:n_per_batch, .combine = c, .packages = "gnm", .inorder = FALSE) %dopar% {
          global_b <- (j - 1) * n_per_batch + b
          set.seed(global_b)
          zip_sample <- sample(1:num_uniq_zip, floor(2 * sqrt(num_uniq_zip)), replace = TRUE)
          aggregate_data_boots <- do.call(rbind, aggregate_data.list[zip_sample])
          
          model <- tryCatch({
            gnm(formula_full,
                eliminate = (sex:race:dual:entry_age_break:followup_year),
                data = aggregate_data_boots,
                family = poisson(link = "log"))
          }, error = function(e) NA)
          
          if (inherits(model, "gnm")) coef(model)[1] else NA
        }
        boot_coefs_all <- c(boot_coefs_all, boot_coefs_part)
        
        iter_count <- iter_count + n_per_batch
        setTxtProgressBar(pb, iter_count)
      }
      
      close(pb)
      
      # 3. CI computation
      boot_coefs_all <- boot_coefs_all[!is.na(boot_coefs_all)]
      
      boot_se_adj <- sd(boot_coefs_all) * sqrt(2 * sqrt(num_uniq_zip)) / sqrt(num_uniq_zip)
      rr <- exp(full_coef * iqr)
      ci_low <- exp((full_coef - 1.96 * boot_se_adj) * iqr)
      ci_high <- exp((full_coef + 1.96 * boot_se_adj) * iqr)
      zval <- full_coef / boot_se_adj
      pval <- 2 * pnorm(-abs(zval))
      
      result_df_list[[length(result_df_list) + 1]] <- data.frame(
        model = model_type,
        urban_var = "all",
        urban_code = 0,
        exp = mexp,
        IQR = "all",
        IQR_val = iqr,
        coef = full_coef,
        se = boot_se_adj,
        RR = rr,
        CI_low = ci_low,
        CI_high = ci_high,
        pval = pval
      )
    }
    
    if (length(result_df_list) > 0) {
      result_df <- do.call(rbind, result_df_list)
      write.csv(result_df,
                paste0("result_boot_", mexp, ".csv"),
                row.names = FALSE)
      cat("Saved: result_boot_", mexp, ".csv\n")
    }
    
  }, error = function(e) { cat("ERROR:", conditionMessage(e), "\n") })
}

stopImplicitCluster()

end_time <- Sys.time()
print(end_time - start_time)


# Combine
file_list <- list.files(pattern = "^result_boot_.*\\.csv$")
result_all <- do.call(rbind, lapply(file_list, read.csv))



###### Urban vs rural #####
registerDoParallel(cores = 8)

urban_var <- c("uapop_2g", "ualand_2g", "poptot_2g", "popden_2g")
urban_code <- c(1, 2)
exp_name <- c("pm25_ensemble","oc","ec","so4","no3","nh4",
              "br","ca","cu","fe","k","ni","pb","si","v","z")

n_per_batch <- 50
n_batch <- 10


start_time <- Sys.time()

for (mvar in urban_var) {
  for (murban in urban_code) {
    tryCatch({
      result_df_list <- list()
      data_subset <- filter(aggregate_data, !!sym(mvar) == murban)
      exposure_subset <- filter(exposure, !!sym(mvar) == murban)
      
      # urban/rural stratification
      data_list_ur <- split(data_subset, data_subset$zip)
      num_uniq_zip <- length(data_list_ur)
      
      total_iter <- n_per_batch * n_batch
      cat("\nBootstrapping", mvar, "group", murban, "\n")
      pb <- txtProgressBar(min = 0, max = total_iter * length(exp_name), style = 3)
      iter_count <- 0
      
      for (mexp in exp_name) {
        mcov <- paste0("PM25_", mexp)
        
        iqr_all <- quantile(get(mexp, exposure), 0.75, na.rm=T) - quantile(get(mexp, exposure), 0.25, na.rm=T)
        iqr_ur <- quantile(get(mexp, exposure_subset), 0.75, na.rm = TRUE) - quantile(get(mexp, exposure_subset), 0.25, na.rm = TRUE)
        
        # Model
        formula_base <- if (mexp == "pm25_ensemble") {
          paste0("dead ~ ", mexp)
        } else {
          paste0("dead ~ ", mexp, " + ", mcov)
        }
        
        formula_full <- as.formula(paste0(formula_base,
                                          " + mean_bmi + smoke_rate + hispanic + pct_blk +",
                                          " medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ +",
                                          " summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +",
                                          " year + region + offset(log(time_count))"))
        
        full_model <- tryCatch({
          gnm(formula_full,
              eliminate = (sex:race:dual:entry_age_break:followup_year),
              data = data_subset,
              family = poisson(link = "log"))
        }, error = function(e) NA)
        
        if (!inherits(full_model, "gnm")) {
          cat("Full model failed:", mexp, mvar, murban, "\n")
          iter_count <- iter_count + total_iter
          setTxtProgressBar(pb, iter_count)
          next
        }
        
        full_coef <- coef(full_model)[1]
        
        boot_coefs_all <- c()
        for (j in 1:n_batch) {
          boot_coefs_part <- foreach(b = 1:n_per_batch, .combine = c,
                                     .packages = "gnm", .inorder = FALSE) %dopar% {
                                       global_b <- (j - 1) * n_per_batch + b
                                       set.seed(global_b)
                                       zip_sample <- sample(1:num_uniq_zip, floor(2 * sqrt(num_uniq_zip)), replace = TRUE)
                                       data_boots <- do.call(rbind, data_list_ur[zip_sample])
                                       
                                       model <- tryCatch({
                                         gnm(formula_full,
                                             eliminate = (sex:race:dual:entry_age_break:followup_year),
                                             data = data_boots,
                                             family = poisson(link = "log"))
                                       }, error = function(e) NA)
                                       
                                       if (inherits(model, "gnm")) coef(model)[1] else NA
                                     }
          
          boot_coefs_all <- c(boot_coefs_all, boot_coefs_part)
          iter_count <- iter_count + n_per_batch
          setTxtProgressBar(pb, iter_count)
        }
        
        boot_coefs_all <- boot_coefs_all[!is.na(boot_coefs_all)]
        
        boot_se_adj <- sd(boot_coefs_all) * sqrt(2 * sqrt(num_uniq_zip)) / sqrt(num_uniq_zip)
        zval <- full_coef / boot_se_adj
        pval <- 2 * pnorm(-abs(zval))
        
        for (iqr_type in c("all", "ur")) {
          iqr_val <- if (iqr_type == "all") iqr_all else iqr_ur
          
          rr <- exp(full_coef * iqr_val)
          ci_low <- exp((full_coef - 1.96 * boot_se_adj) * iqr_val)
          ci_high <- exp((full_coef + 1.96 * boot_se_adj) * iqr_val)
          
          
          result_df_list[[length(result_df_list) + 1]] <- data.frame(
            model = 2,
            urban_var = mvar,
            urban_code = murban,
            exp = mexp,
            IQR = iqr_type,
            IQR_val = iqr_val,
            coef = full_coef,
            se = boot_se_adj,
            RR = rr,
            CI_low = ci_low,
            CI_high = ci_high,
            pval = pval
          )
        }
      }
      
      close(pb)
      
      if (length(result_df_list) > 0) {
        result_df <- do.call(rbind, result_df_list)
        file_name <- paste0("result_boot_", mvar, "_", murban, ".csv")
        write.csv(result_df, file_name, row.names = FALSE)
        cat("Saved:", file_name, "\n")
      }
      
    }, error = function(e) {
      cat("ERROR:", mvar, murban, "-", conditionMessage(e), "\n")
    })
  }
}

stopImplicitCluster()
end_time <- Sys.time()
print(end_time - start_time)


file_list <- list.files(pattern = "^result_boot_.*_\\d\\.csv$")
result_urban <- do.call(rbind, lapply(file_list, read.csv))



##########################################################
############ Loop - individual characteristics ##########
##########################################################

##### All area #####
aggregate_data <- aggregate_data %>%
  mutate(age = ifelse(entry_age_break %in% c(1,2),1, 
                      ifelse(entry_age_break %in% c(3:8), 2, NA))) 

registerDoParallel(cores = 8)

exp_name <- c("pm25_ensemble")
pop_var <- c("age", "sex", "race", "dual")
strata <- c(0, 1, 2, 4, 5, 6, 9)

n_per_batch <- 50
n_batch <- 10
total_iter <- n_per_batch * n_batch

# IQR
iqr_values <- setNames(
  map(exp_name, ~ {
    if (.x %in% colnames(exposure)) {
      quantile(exposure[[.x]], 0.75, na.rm = TRUE) -
        quantile(exposure[[.x]], 0.25, na.rm = TRUE)
    } else {
      NA
    }
  }),
  exp_name
)

# bootstrap function
bootstrap_se <- function(formula_full, subset_data, n_per_batch, n_batch) {
  subset_data.list <- split(subset_data, subset_data$zip)
  num_zip_subset <- length(subset_data.list)
  
  boot_coefs_all <- c()
  for (j in 1:n_batch) {
    boot_part <- foreach(b = 1:n_per_batch, .combine = c, .packages = "gnm", .inorder = FALSE) %dopar% {
      set.seed((j - 1) * n_per_batch + b)
      zip_sample <- sample(1:num_zip_subset, floor(2 * sqrt(num_zip_subset)), replace = TRUE)
      data_boot <- do.call(rbind, subset_data.list[zip_sample])
      
      model <- tryCatch({
        gnm(formula_full,
            eliminate = (sex:race:dual:entry_age_break:followup_year),
            data = data_boot,
            family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (inherits(model, "gnm")) coef(model)[1] else NA
    }
    boot_coefs_all <- c(boot_coefs_all, boot_part)
  }
  
  boot_coefs_all <- boot_coefs_all[!is.na(boot_coefs_all)]
  sd(boot_coefs_all) * sqrt(2 * sqrt(num_zip_subset)) / sqrt(num_zip_subset)
}

start_time <- Sys.time()
# Main loop
for (current_pop_var in pop_var) {
  result_list <- list()
  pb <- txtProgressBar(min = 0, max = length(exp_name) * length(strata), style = 3)
  iter <- 0
  
  for (mexp in exp_name) {
    iqr_val <- iqr_values[[mexp]]
    
    for (current_pop_strata in strata) {
      subset_data <- filter(aggregate_data, !!sym(current_pop_var) == current_pop_strata)
      if (nrow(subset_data) < 100) {
        iter <- iter + 1; setTxtProgressBar(pb, iter); next
      }
      
      # Model
      formula_full <- as.formula(paste0(
        "dead ~ ", mexp,
        " + mean_bmi + smoke_rate + hispanic + pct_blk +",
        " medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ +",
        " summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +",
        " year + region + offset(log(time_count))"
      ))
      
      full_model <- tryCatch({
        gnm(formula_full,
            eliminate = (sex:race:dual:entry_age_break:followup_year),
            data = subset_data,
            family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (!inherits(full_model, "gnm")) {
        iter <- iter + 1; setTxtProgressBar(pb, iter); next
      }
      
      coef_est <- coef(full_model)[1]
      se_boot <- bootstrap_se(formula_full, subset_data, n_per_batch, n_batch)
      
      if (is.na(se_boot)) {
        iter <- iter + 1; setTxtProgressBar(pb, iter); next
      }
      
      rr <- exp(coef_est * iqr_val)
      ci_low <- exp((coef_est - 1.96 * se_boot) * iqr_val)
      ci_high <- exp((coef_est + 1.96 * se_boot) * iqr_val)
      pval <- 2 * pnorm(-abs(coef_est / se_boot))
      
      result_list[[length(result_list) + 1]] <- tibble(
        model = 1,
        urban_var = "all",
        urban_code = 0,
        exp = mexp,
        pop_var = current_pop_var,
        pop_strata = current_pop_strata,
        IQR = iqr_val,
        coef = coef_est,
        se = se_boot,
        RR = rr,
        CI_low = ci_low,
        CI_high = ci_high,
        pval = pval
      )
      
      iter <- iter + 1
      setTxtProgressBar(pb, iter)
    }
  }
  
  close(pb)
  
  if (length(result_list) > 0) {
    df_out <- bind_rows(result_list)
    write.csv(df_out, paste0("result_boot_by_", current_pop_var, ".csv"), row.names = FALSE)
    cat("Saved: result_boot_by_", current_pop_var, ".csv\n")
  }
}

stopImplicitCluster()
end_time <- Sys.time()
print(end_time - start_time)


file_list <- list.files(pattern = "^result_boot_by_.*\\.csv$")
result_all <- do.call(rbind, lapply(file_list, read.csv))






###### Urban vs rural #####
registerDoParallel(cores = 8)

exp_name <- c("pm25_ensemble")
pop_var <- c("age", "sex", "race", "dual")
strata <- c(0, 1, 2, 4, 5, 6, 9)

urban_var <- c("RUCA_2g", "uapop_2g", "ualand_2g", "poptot_2g", "popden_2g")
ur_codes <- c(1, 2)

n_per_batch <- 50
n_batch <- 10

# IQR
iqr_values <- setNames(
  map(exp_name, ~ {
    if (.x %in% colnames(exposure)) {
      quantile(exposure[[.x]], 0.75, na.rm = TRUE) -
        quantile(exposure[[.x]], 0.25, na.rm = TRUE)
    } else {
      NA
    }
  }),
  exp_name
)

# bootstrap function
bootstrap_se <- function(formula_full, subset_data, n_per_batch, n_batch) {
  subset_data.list <- split(subset_data, subset_data$zip)
  num_zip_subset <- length(subset_data.list)
  
  boot_coefs_all <- c()
  for (j in 1:n_batch) {
    boot_part <- foreach(b = 1:n_per_batch, .combine = c, .packages = "gnm", .inorder = FALSE) %dopar% {
      set.seed((j - 1) * n_per_batch + b)
      zip_sample <- sample(1:num_zip_subset, floor(2 * sqrt(num_zip_subset)), replace = TRUE)
      data_boot <- do.call(rbind, subset_data.list[zip_sample])
      
      model <- tryCatch({
        gnm(formula_full,
            eliminate = (sex:race:dual:entry_age_break:followup_year),
            data = data_boot,
            family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (inherits(model, "gnm")) coef(model)[1] else NA
    }
    boot_coefs_all <- c(boot_coefs_all, boot_part)
  }
  
  boot_coefs_all <- boot_coefs_all[!is.na(boot_coefs_all)]
  sd(boot_coefs_all) * sqrt(2 * sqrt(num_zip_subset)) / sqrt(num_zip_subset)
}

start_time <- Sys.time()

# Main loop (by urban_var)
for (uvar in urban_var) {
  result_list <- list()
  pb <- txtProgressBar(min = 0, max = length(exp_name) * length(ur_codes) * length(pop_var) * length(strata), style = 3)
  iter <- 0
  
  for (mexp in exp_name) {
    iqr_val <- iqr_values[[mexp]]
    
    for (ucode in ur_codes) {
      for (current_pop_var in pop_var) {
        for (current_pop_strata in strata) {
          
          subset_data <- filter(aggregate_data, !!sym(uvar) == ucode, !!sym(current_pop_var) == current_pop_strata)
          if (nrow(subset_data) < 100) {
            iter <- iter + 1; setTxtProgressBar(pb, iter); next
          }
          
          formula_full <- as.formula(paste0(
            "dead ~ ", mexp,
            " + mean_bmi + smoke_rate + hispanic + pct_blk +",
            " medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ +",
            " summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +",
            " year + region + offset(log(time_count))"
          ))
          
          full_model <- tryCatch({
            gnm(formula_full,
                eliminate = (sex:race:dual:entry_age_break:followup_year),
                data = subset_data,
                family = poisson(link = "log"))
          }, error = function(e) NA)
          
          if (!inherits(full_model, "gnm")) {
            iter <- iter + 1; setTxtProgressBar(pb, iter); next
          }
          
          coef_est <- coef(full_model)[1]
          se_boot <- bootstrap_se(formula_full, subset_data, n_per_batch, n_batch)
          
          if (is.na(se_boot)) {
            iter <- iter + 1; setTxtProgressBar(pb, iter); next
          }
          
          rr <- exp(coef_est * iqr_val)
          ci_low <- exp((coef_est - 1.96 * se_boot) * iqr_val)
          ci_high <- exp((coef_est + 1.96 * se_boot) * iqr_val)
          pval <- 2 * pnorm(-abs(coef_est / se_boot))
          
          result_list[[length(result_list) + 1]] <- tibble(
            model = 1,
            urban_var = uvar,
            urban_code = ucode,
            exp = mexp,
            pop_var = current_pop_var,
            pop_strata = current_pop_strata,
            IQR = iqr_val,
            coef = coef_est,
            se = se_boot,
            RR = rr,
            CI_low = ci_low,
            CI_high = ci_high,
            pval = pval
          )
          
          iter <- iter + 1
          setTxtProgressBar(pb, iter)
        }
      }
    }
  }
  
  close(pb)
  
  if (length(result_list) > 0) {
    df_out <- bind_rows(result_list)
    write.csv(df_out, paste0("result_boot_by_", uvar, ".csv"), row.names = FALSE)
    cat("Saved: result_boot_by_", uvar, ".csv\n")
  }
}

stopImplicitCluster()
end_time <- Sys.time()
print(end_time - start_time)


file_list <- list.files(pattern = "^result_boot_by_.*\\.csv$")
result_all <- do.call(rbind, lapply(file_list, read.csv))




############################################
##### Loop - population characteristics ####
############################################

aggregate_data <- aggregate_data %>%
  mutate(popdensity_2g=ntile(popdensity,2),
         hispanic_2g=ntile(hispanic,2),
         pct_blk_2g=ntile(pct_blk,2),
         medhouseholdincome_2g=ntile(medhouseholdincome,2),
         medianhousevalue_2g=ntile(medianhousevalue,2),
         poverty_2g=ntile(poverty,2),
         education_2g=ntile(education,2),
         pct_owner_occ_2g=ntile(pct_owner_occ,2),
         mean_bmi_2g=ntile(mean_bmi,2),
         smoke_rate_2g=ntile(smoke_rate,2))


##### All area #####

registerDoParallel(cores = 8)

exp_name <- c("pm25_ensemble")
pop_var <- c("hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue",
             "poverty", "education", "pct_owner_occ", "mean_bmi", "smoke_rate")
strata <- c(1, 2)

n_per_batch <- 50
n_batch <- 10

# IQR
iqr_values <- setNames(
  map(exp_name, ~ {
    if (.x %in% colnames(exposure)) {
      quantile(exposure[[.x]], 0.75, na.rm = TRUE) -
        quantile(exposure[[.x]], 0.25, na.rm = TRUE)
    } else {
      NA
    }
  }),
  exp_name
)

# bootstrap function
bootstrap_se <- function(formula_full, subset_data, n_per_batch, n_batch) {
  subset_data.list <- split(subset_data, subset_data$zip)
  num_zip_subset <- length(subset_data.list)
  
  boot_coefs_all <- c()
  for (j in 1:n_batch) {
    boot_part <- foreach(b = 1:n_per_batch, .combine = c, .packages = "gnm", .inorder = FALSE) %dopar% {
      set.seed((j - 1) * n_per_batch + b)
      zip_sample <- sample(1:num_zip_subset, floor(2 * sqrt(num_zip_subset)), replace = TRUE)
      data_boot <- do.call(rbind, subset_data.list[zip_sample])
      
      model <- tryCatch({
        gnm(formula_full,
            eliminate = (sex:race:dual:entry_age_break:followup_year),
            data = data_boot,
            family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (inherits(model, "gnm")) coef(model)[1] else NA
    }
    boot_coefs_all <- c(boot_coefs_all, boot_part)
  }
  
  boot_coefs_all <- boot_coefs_all[!is.na(boot_coefs_all)]
  sd(boot_coefs_all) * sqrt(2 * sqrt(num_zip_subset)) / sqrt(num_zip_subset)
}

start_time <- Sys.time()

# Main loop (by pop_var)
for (current_pop_var in pop_var) {
  result_list <- list()
  pb <- txtProgressBar(min = 0, max = length(exp_name) * length(strata), style = 3)
  iter <- 0
  
  for (mexp in exp_name) {
    iqr_val <- iqr_values[[mexp]]
    
    for (current_pop_strata in strata) {
      subset_data <- filter(aggregate_data, !!sym(paste0(current_pop_var, "_2g")) == current_pop_strata)
      if (nrow(subset_data) < 100) {
        iter <- iter + 1; setTxtProgressBar(pb, iter); next
      }
      
      formula_full <- as.formula(paste0(
        "dead ~ ", mexp,
        " + mean_bmi + smoke_rate + hispanic + pct_blk +",
        " medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ +",
        " summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +",
        " year + region + offset(log(time_count))"
      ))
      
      # base model to drop current covariate
      formula_out <- update(formula_full, paste0(". ~ . - ", current_pop_var))
      
      full_model <- tryCatch({
        gnm(formula_out,
            eliminate = (sex:race:dual:entry_age_break:followup_year),
            data = subset_data,
            family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (!inherits(full_model, "gnm")) {
        iter <- iter + 1; setTxtProgressBar(pb, iter); next
      }
      
      coef_est <- coef(full_model)[1]
      se_boot <- bootstrap_se(formula_out, subset_data, n_per_batch, n_batch)
      
      if (is.na(se_boot)) {
        iter <- iter + 1; setTxtProgressBar(pb, iter); next
      }
      
      rr <- exp(coef_est * iqr_val)
      ci_low <- exp((coef_est - 1.96 * se_boot) * iqr_val)
      ci_high <- exp((coef_est + 1.96 * se_boot) * iqr_val)
      pval <- 2 * pnorm(-abs(coef_est / se_boot))
      
      result_list[[length(result_list) + 1]] <- tibble(
        model = 1,
        urban_var = "all",
        urban_code = 0,
        exp = mexp,
        pop_var = current_pop_var,
        pop_strata = current_pop_strata,
        IQR = iqr_val,
        coef = coef_est,
        se = se_boot,
        RR = rr,
        CI_low = ci_low,
        CI_high = ci_high,
        pval = pval
      )
      
      iter <- iter + 1
      setTxtProgressBar(pb, iter)
    }
  }
  
  close(pb)
  
  if (length(result_list) > 0) {
    df_out <- bind_rows(result_list)
    write.csv(df_out, paste0("result_boot_by_", current_pop_var, ".csv"), row.names = FALSE)
    cat("Saved: result_boot_by_", current_pop_var, ".csv\n")
  }
}

stopImplicitCluster()
end_time <- Sys.time()
print(end_time - start_time)


file_list <- list.files(pattern = "^result_boot_by_.*\\.csv$")
result_all <- do.call(rbind, lapply(file_list, read.csv))



##### Urban vs Rural #####

registerDoParallel(cores = 8)

exp_name <- c("pm25_ensemble")
pop_var <- c("hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue",
             "poverty", "education", "pct_owner_occ", "mean_bmi", "smoke_rate")
strata <- c(1, 2)

urban_var <- c("RUCA_2g","uapop_2g", "ualand_2g","poptot_2g","popden_2g")
urban_code <- c(1, 2)


n_per_batch <- 50
n_batch <- 10

# IQR
iqr_values <- setNames(
  map(exp_name, ~ {
    if (.x %in% colnames(exposure)) {
      quantile(exposure[[.x]], 0.75, na.rm = TRUE) -
        quantile(exposure[[.x]], 0.25, na.rm = TRUE)
    } else {
      NA
    }
  }),
  exp_name
)

# bootstrap function
bootstrap_se <- function(formula_full, subset_data, n_per_batch, n_batch) {
  subset_data.list <- split(subset_data, subset_data$zip)
  num_zip_subset <- length(subset_data.list)
  
  boot_coefs_all <- c()
  for (j in 1:n_batch) {
    boot_part <- foreach(b = 1:n_per_batch, .combine = c, .packages = "gnm", .inorder = FALSE) %dopar% {
      set.seed((j - 1) * n_per_batch + b)
      zip_sample <- sample(1:num_zip_subset, floor(2 * sqrt(num_zip_subset)), replace = TRUE)
      data_boot <- do.call(rbind, subset_data.list[zip_sample])
      
      model <- tryCatch({
        gnm(formula_full,
            eliminate = (sex:race:dual:entry_age_break:followup_year),
            data = data_boot,
            family = poisson(link = "log"))
      }, error = function(e) NA)
      
      if (inherits(model, "gnm")) coef(model)[1] else NA
    }
    boot_coefs_all <- c(boot_coefs_all, boot_part)
  }
  
  boot_coefs_all <- boot_coefs_all[!is.na(boot_coefs_all)]
  sd(boot_coefs_all) * sqrt(2 * sqrt(num_zip_subset)) / sqrt(num_zip_subset)
}

start_time <- Sys.time()

# Main loop (by urban_var)
for (uvar in urban_var) {
  result_list <- list()
  pb <- txtProgressBar(min = 0, max = length(exp_name) * length(urban_code) * length(pop_var) * length(strata), style = 3)
  iter <- 0
  
  for (mexp in exp_name) {
    iqr_val <- iqr_values[[mexp]]
    
    for (ucode in urban_code) {
      for (current_pop_var in pop_var) {
        for (current_pop_strata in strata) {
          
          subset_data <- filter(aggregate_data, !!sym(uvar) == ucode, !!sym(paste0(current_pop_var, "_2g")) == current_pop_strata)
          if (nrow(subset_data) < 100) {
            iter <- iter + 1; setTxtProgressBar(pb, iter); next
          }
          
          formula_full <- as.formula(paste0(
            "dead ~ ", mexp,
            " + mean_bmi + smoke_rate + hispanic + pct_blk +",
            " medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ +",
            " summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +",
            " year + region + offset(log(time_count))"
          ))
          
          formula_out <- update(formula_full, paste0(". ~ . - ", current_pop_var))
          
          full_model <- tryCatch({
            gnm(formula_out,
                eliminate = (sex:race:dual:entry_age_break:followup_year),
                data = subset_data,
                family = poisson(link = "log"))
          }, error = function(e) NA)
          
          if (!inherits(full_model, "gnm")) {
            iter <- iter + 1; setTxtProgressBar(pb, iter); next
          }
          
          coef_est <- coef(full_model)[1]
          se_boot <- bootstrap_se(formula_out, subset_data, n_per_batch, n_batch)
          
          if (is.na(se_boot)) {
            iter <- iter + 1; setTxtProgressBar(pb, iter); next
          }
          
          rr <- exp(coef_est * iqr_val)
          ci_low <- exp((coef_est - 1.96 * se_boot) * iqr_val)
          ci_high <- exp((coef_est + 1.96 * se_boot) * iqr_val)
          pval <- 2 * pnorm(-abs(coef_est / se_boot))
          
          result_list[[length(result_list) + 1]] <- tibble(
            model = 1,
            urban_var = uvar,
            urban_code = ucode,
            exp = mexp,
            pop_var = current_pop_var,
            pop_strata = current_pop_strata,
            IQR = iqr_val,
            coef = coef_est,
            se = se_boot,
            RR = rr,
            CI_low = ci_low,
            CI_high = ci_high,
            pval = pval
          )
          
          iter <- iter + 1
          setTxtProgressBar(pb, iter)
        }
      }
    }
  }
  
  close(pb)
  
  if (length(result_list) > 0) {
    df_out <- bind_rows(result_list)
    write.csv(df_out, paste0("result_boot_by_", uvar, ".csv"), row.names = FALSE)
    cat("Saved: result_boot_by_", uvar, ".csv\n")
  }
}

stopImplicitCluster()
end_time <- Sys.time()
print(end_time - start_time)


file_list <- list.files(pattern = "^result_boot_by_.*\\.csv$")
result_all <- do.call(rbind, lapply(file_list, read.csv))
