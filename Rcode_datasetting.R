library(tidyverse)
library(arrow)
library(data.table)
library(zoo)
library(fst)
library(survival)
library(gnm)
library(ggplot2)
library(sf)
library(purrr)
library(data.table)


### DATA: PM25 ###
years <- 2000:2016
pm25 <- map_df(years, ~ {
  file_path <- file.path("/n/dominici_nsaph_l3/Lab/exposure/pm25/PM25_v2/annual/", paste0(.x, ".rds"))
  readRDS(file_path) %>% select(-STATE) %>%
    mutate(year = .x) # Add the year column
})


### DATA: denominator ###
file_list <- list.files(path="/n/dominici_nsaph_l3/Lab/data/data_warehouse/dw_legacy_medicare_00_16",
                        pattern="bene_.*\\.parquet", full.names=TRUE)
data_list <- lapply(file_list, read_parquet)
denom <- do.call(rbind, data_list) %>% arrange(bene_id, year)


denom <- as.data.table(denom)
pm25 <- as.data.table(pm25)

# Perform operations
denom_clean <- denom[
  !is.na(bene_id) &                                 # Filter: `bene_id` is not NA
    !state %in% c("AK", "HI", "PR") &               # Exclude non-contiguous states
    sprintf("%05d", zip) %in% pm25$ZIP &            # Keep ZIPs present in `pm25` (zero-padded)
    age_dob >= 65 & age_dob <= 150                  # Filter age between 65 and 150
][
  , .SD[!duplicated(year)], by = bene_id            # Remove duplicates by `year` within each `bene_id`
][
  , `:=`(
    race = if (n_distinct(race) > 1) 0 else first(race), # Assign `race = 0` if inconsistent within the group
    race = if (any(race %in% c(0, 3))) 9 else race      # Assign `race = 9` if any `race` is 0 or 3
  ), by = bene_id
][
  , if (
    n_distinct(sex) == 1 &&                         # `sex` is consistent
    n_distinct(dod[!is.na(dod)]) <= 1               # `dod` is unique or all are NA
  ) .SD, by = bene_id
][
  , `:=`(
    entry_year = min(year),                         # Earliest year for each `bene_id`
    followup_year = year - min(year),               # Follow-up year calculation
    followup_year_plus_one = year - min(year) + 1,  # Follow-up year + 1
    entry_age = age_dob[which.min(year)]            # Entry age at earliest year
  ), by = bene_id
][
  , entry_age_break := fifelse(entry_age >= 65 & entry_age <= 69, 1,
                               fifelse(entry_age >= 70 & entry_age <= 74, 2,
                                       fifelse(entry_age >= 75 & entry_age <= 79, 3,
                                               fifelse(entry_age >= 80 & entry_age <= 84, 4,
                                                       fifelse(entry_age >= 85 & entry_age <= 89, 5,
                                                               fifelse(entry_age >= 90 & entry_age <= 94, 6,
                                                                       fifelse(entry_age >= 95 & entry_age <= 99, 7, 8)))))))
][
  , dead := fifelse(is.na(dod), 0, 1)               # Assign `dead` flag (1 = dead, 0 = alive)
]


saveRDS(denom_clean, "/n/dominici_nsaph_l3/Lab/projects/garambyun_pm25-longterm-mortality-urbanrural/Data/denom_clean.rds")




################# Exposure ###################
# PM25
years <- 2000:2016
pm25 <- map_df(years, ~ {
  file_path <- file.path("/n/dominici_nsaph_l3/Lab/exposure/pm25/PM25_v2/annual/", paste0(.x, ".rds"))
  readRDS(file_path) %>% select(-STATE) %>%
    mutate(year = .x) %>% # Add the year column
    rename(pm25_ensemble = pm25)
})

# Ozone
ozone <- map_df(years, ~ {
  file_path <- file.path("/n/dominici_nsaph_l3/Lab/exposure/ozone/O3_v2/annual/", paste0(.x, ".rds"))
  readRDS(file_path) %>% select(-STATE) %>%
    mutate(year = .x) # Add the year column
})

# Temperature
temp <- read.csv("/n/dominici_nsaph_l3/Lab/projects/analytic/temperature_seasonal_zipcode/temperature_seasonal_zipcode_combined.csv") %>%
  mutate(ZIP=sprintf("%05d", ZIP))


# PM25 component
data_list <- list()

for (year in 2000:2016) {
  file_path <- paste0("/n/dominici_nsaph_l3/Lab/exposure/pm25_components/pm25_components_v2/PM25_components/", year, ".rds")
  data_list[[as.character(year)]] <- readRDS(file_path) %>%
    mutate(year = year)
}

pmcom <- bind_rows(data_list)


# Area variables
covar <- read.csv("/n/dominici_nsaph_l3/Lab/projects/analytic/merged_covariates_pm_census_temp/merged_covariates.csv") %>%
  select("zip", "year","mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue",
         "poverty", "education", "popdensity", "pct_owner_occ") %>%
  mutate(ZIP=sprintf("%05d", zip))


# Merge
exposure <- pm25 %>%
  left_join(ozone, by = c("ZIP","year")) %>%
  left_join(temp, by = c("ZIP","year")) %>%
  left_join(covar, by = c("ZIP","year")) %>%
  left_join(pmcom, by = c("ZIP","year"))
  
exposure <- exposure[complete.cases(exposure), ]


################# Urban #################
RUCA <- readRDS("/n/dominici_nsaph_l3/Lab/projects/garambyun_pm25-longterm-mortality-urbanrural/Data/RUCA_ZIP.rds")
Census <- readRDS("/n/dominici_nsaph_l3/Lab/projects/garambyun_pm25-longterm-mortality-urbanrural/Data/Census_zcta.rds")

ztca <- read.csv("/n/dominici_nsaph_l3/Lab/exposure/zip2zcta_master_xwalk/zip2zcta_master_xwalk.csv") %>%
  mutate(ZIP=sprintf("%05d", zip)) %>%
  select(ZIP,zip,zcta)

Urban_indicator <- merge(Census, ztca, by="zcta", all=F) %>%
  left_join(RUCA, by="ZIP") %>% select(ZIP,uapop_2g,uapop_5g,ualand_2g,ualand_5g,
                                       poptot_2g,poptot_5g,popden_2g,popden_5g,RUCA_2g,RUCA_5g) %>%
  mutate(
    RUCA_2g = ifelse(is.na(RUCA_2g), 2, RUCA_2g),
    RUCA_5g = ifelse(is.na(RUCA_5g), 5, RUCA_5g)
  )




################# Aggregate ###################
denom_clean <- readRDS("/n/dominici_nsaph_l3/Lab/projects/garambyun_pm25-longterm-mortality-urbanrural/Data/denom_clean.rds")

denom_clean$time_count <- denom_clean$followup_year_plus_one - denom_clean$followup_year
dead_personyear <- aggregate(cbind(denom_clean$dead,
                                   denom_clean$time_count), 
                             by=list(denom_clean$zip,
                                     denom_clean$year,
                                     denom_clean$sex,
                                     denom_clean$race,
                                     denom_clean$dual,
                                     denom_clean$entry_age_break,
                                     denom_clean$followup_year),
                             FUN = sum)

colnames(dead_personyear)[8:10] <- c("dead", "time_count")
colnames(dead_personyear)[1:7] <- c("zip", "year", "sex", "race", "dual", "entry_age_break", "followup_year")

states_zips <- unique(denom_clean[, c("state", "zip")])

aggregate_data <- merge(dead_personyear, exposure, by = c("zip","year"), all=FALSE) %>%
  left_join(states_zips, by = c("zip")) %>%
  left_join(Urban_indicator, by = c("ZIP")) 



NORTHEAST <- c("NY", "MA", "PA", "RI", "NH", "ME", "VT", "CT", "NJ")  
SOUTH <- c("DC", "VA", "NC", "WV", "KY", "SC", "GA", "FL", "AL", "TN", "MS", "AR", "MD", "DE", "OK", "TX", "LA")
MIDWEST <- c("OH", "IN", "MI", "IA", "MO", "WI", "MN", "SD", "ND", "IL", "KS", "NE")
WEST <- c("MT", "CO", "WY", "ID", "UT", "NV", "CA", "OR", "WA", "AZ", "NM")

aggregate_data$region <- ifelse(aggregate_data$state %in% NORTHEAST, "NORTHEAST",
                                ifelse(aggregate_data$state %in% SOUTH, "SOUTH",
                                       ifelse(aggregate_data$state %in% MIDWEST, "MIDWEST",
                                              ifelse(aggregate_data$state %in% WEST, "WEST",
                                                     NA))))