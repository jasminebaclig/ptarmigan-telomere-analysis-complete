#Combines RTL data with ptarmigan data from Molly and collaborators
#Author: Jasmine Baclig

library(tidyverse)



#Adds sex of birds to relative telomere length data (Took sex data from Molly's scripts and CSVs)
lm_data <- left_join(calculated_data %>% select(id, rtl), all_samples %>% select(value, sex, age), by = join_by(id == value))



#Imports data provided by Kristinn
kpm_data_1 <- read_csv("./lm_info_files/RP_data.csv")
kpm_data_2 <- read_csv("./lm_info_files/RP_data2.csv")

#Reformats bird ID on each table
kpm_data_1 <- mutate(kpm_data_1, id = substr(id, 4, 9))
kpm_data_2 <- mutate(kpm_data_2, BirdID = substr(BirdID, 4, 9))

#Selects desired columns from data
kpm_data_1 <- select(kpm_data_1, -sex, -age, -year, -densj)
kpm_data_2 <- select(kpm_data_2, -Latitude, -Longitude, -CollYear, -CollMonth, -CollDay, -DatePrec, -CollHour, -CollMinute, -Sex, -Age)

#Combines Kristinn's data for selected birds
combined_kpm <- full_join(kpm_data_1, kpm_data_2, by = join_by(id == BirdID))



#Imports data provided by Amanda
amanda_data_1 <- read_csv("./lm_info_files/amanda_data_1.csv")
amanda_data_2 <- read_csv("./lm_info_files/Diet_diversity.csv")
amanda_data_3 <- read_csv("./lm_info_files/bci.csv")
amanda_data_4 <- read_csv("./lm_info_files/Ptarmigan_micro_morpho_diversity_AMS.csv") %>% rename(id = "Bird ID")

#Reformats bird ID on each table
amanda_data_1 <- mutate(amanda_data_1, BirdID = substr(BirdID, 4, 9))
amanda_data_2 <- mutate(amanda_data_2, bird_ID = substr(bird_ID, 4, 9))
amanda_data_3 <- mutate(amanda_data_3, id = substr(id, 4, 9))
amanda_data_4 <- mutate(amanda_data_4, id = substr(id, 4, 9))

#Selects desired columns from data
amanda_data_1 <- select(amanda_data_1, -"Subsample Barcode", -"Notes...2", -Project, -Fate, -Phase, -CollYear, -CollMonth, -CollDay, -Sex, -Age, -"Notes...18")
amanda_data_2 <- select(amanda_data_2, -Year)
amanda_data_3 <- select(amanda_data_3, id, Wing, Head, TarsMT, PC1, "BC-index", "Standardized BC-index")
amanda_data_4 <- select(amanda_data_4, -"...1", -Weight, -"Ceca:Weight")

#Combines Amanda's data for selected birds
combined_amanda <- full_join(amanda_data_1, amanda_data_2, by = join_by(BirdID == bird_ID)) %>%
                   full_join(amanda_data_3, by = join_by(BirdID == id)) %>%
                   full_join(amanda_data_4, by = join_by(BirdID == id)) %>%
                   filter(!is.na(BirdID))



#Adds year to data
get_year <- function(prefix) {
  year <- switch(prefix, "06" = "2006", "07" = "2007", "08" = "2008", "09" = "2009", "10" = "2010", "11" = "2011",
                 "12" = "2012", "13" = "2013", "14" = "2014", "15" = "2015", "16" = "2016", "17" = "2017", "18" = "2018")
  return(year)
}

lm_data <- transform(lm_data, year = sapply(substr(id, 1, 2), get_year))
combined_kpm <- transform(combined_kpm, year = sapply(substr(id, 1, 2), get_year))
combined_amanda <- transform(combined_amanda, year = sapply(substr(BirdID, 1, 2), get_year))

#Reformats year column
lm_data$year = as.factor(lm_data$year)
combined_kpm$year = as.factor(combined_kpm$year)
combined_amanda$year = as.factor(combined_amanda$year)



#Imports population density data
pop_density <- read_csv("./lm_info_files/population_over_time.csv")



#Cleans up unwanted data in memory
rm(get_year)
rm(kpm_data_1, kpm_data_2)
rm(amanda_data_1, amanda_data_2, amanda_data_3, amanda_data_4)
rm(all_samples, calculated_data)
