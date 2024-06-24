# Takes exported data from LightCycler and calculates relative telomere length for each bird
# Author: Jasmine Baclig

library(tidyverse)

# Stores path for required files
txt_path <- "./txt_files/"
effic_path <- "./indiv_effic_files/"

# Gets file names for all required files
txt_files <- list.files(path = txt_path) %>% sort()
effic_files <- list.files(path = effic_path) %>% sort()

all_data <- data.frame()

# Combines all data into one data frame
for(i in 1:length(txt_files)) {
  df_txt <- read_tsv(paste(txt_path, txt_files[[i]], sep = "")) %>% select("Position", "Sample Name", "Gene Name", "Cq Mean", "Cq Error")
  df_effic <- read_tsv(paste(effic_path, effic_files[[i]], sep = "")) %>% select("well", "indiv PCR eff")
  all_data <- rbind(all_data, left_join(df_txt, df_effic, by = join_by(x$"Position" == y$"well"))
                              %>% select(-"Position")
                              %>% mutate(run = str_sub(txt_files[[i]], end = -5)))
}

# Cleans up data frame
all_data <- rename(all_data, id = "Sample Name", gene = "Gene Name", sample_cq = "Cq Mean", cq_error = "Cq Error", sample_effic = "indiv PCR eff")
reduced_data <- filter(all_data, !is.na(as.numeric(cq_error)), cq_error < 0.51, !grepl("Negative", id)) %>% select(-cq_error)
avg_data <- group_by(reduced_data, run, id, gene, sample_cq) %>% summarise_at(vars(sample_effic), list(sample_effic = mean))
non_single_data <- count(reduced_data, run, id, gene) %>% filter(n != 1) %>%
                   select(-n) %>%
                   left_join(avg_data, by = join_by(run, id, gene))


## FOR DEBUGGING PURPOSES ##
temp <- filter(non_single_data, is.nan(sample_effic))
good_effic_data <- filter(non_single_data, !is.nan(sample_effic))


# Adds appropriate GB data to sample rows
gb_data <- filter(good_effic_data, id == "GB") %>% select(run, gene, sample_cq, sample_effic) %>% rename(gb_cq = sample_cq, gb_effic = sample_effic)
sample_data <- left_join(good_effic_data %>% filter(id != "GB"), gb_data, by = join_by(run, gene))


## FOR DEBUGGING PURPOSES ##
temp_2 <- count(sample_data, id, gene) %>% filter(n != 1)


# Combines TOX and TELO data for each sample
sample_id <- sample_data["id"] %>% distinct()
tox_telo_data <- data.frame()
temp_3 <- data.frame() ## FOR DEBUGGING PURPOSES ##

for(i in 1:length(sample_id[[1]])) {
  tox_data <- filter(sample_data, id == sample_id[[1]][i], gene == "TOX")
  telo_data <- filter(sample_data, id == sample_id[[1]][i], gene == "TELO")
  
  if(dim(tox_data)[[1]] != 0 && dim(telo_data)[[1]] != 0) {
    new_row <- data.frame(id = sample_id[[1]][i], tox_sample_cq = tox_data["sample_cq"][[1]], tox_sample_effic = tox_data["sample_effic"][[1]],
                                               tox_gb_cq = tox_data["gb_cq"][[1]], tox_gb_effic = tox_data["gb_effic"][[1]],
                                               telo_sample_cq = telo_data["sample_cq"][[1]], telo_sample_effic = telo_data["sample_effic"][[1]],
                                               telo_gb_cq = telo_data["gb_cq"][[1]], telo_gb_effic = telo_data["gb_effic"][[1]])
    tox_telo_data <- rbind(tox_telo_data, new_row)
  } else { ## FOR DEBUGGING PURPOSES ##
    temp_3 <- rbind(temp_3, sample_id[[1]][i])
  }
}


# #Calculates relative telomere length
# calculate_rtl <- function(efficiency_telo, gb_telo, cq_telo, efficiency_tox, gb_tox, cq_tox) {
#   telo_data <- efficiency_telo ^ (gb_telo - cq_telo)
#   tox_data <- efficiency_tox ^ (gb_tox - cq_tox)
#   return(telo_data / tox_data)
# }
# 
# combined_data <- transform(combined_data, rtl = calculate_rtl(efficiency_telo, gb_telo, cq_telo, efficiency_tox, gb_tox, cq_tox))
# 
# 
# 
# #Cleans up unwanted data in memory
# rm(female_adult, female_juvenile, male_adult, male_juvenile)
# rm(age_sex, alphanumeric, df, lm2006, lm2007, lm2008, lm2009, lm2010, lm2011, lm2012, lm2013, lm2014, lm2015, lm2016, lm2017, lm2018, lmxtra, raw_plate_data)
# rm(adult, df_name, female, juvenile, letters_a_h, male, numbers, plates_names)
# 
# rm(plate_names, efficiencies, i)
# rm(tox_1, tox_2, tox_3, telo_1, telo_2, telo_3, redo_1_tox, redo_1_telo, redo_2_tox, redo_2_telo)
# rm(plate_1, plate_2, plate_3, redo_1, redo_2)
# rm(combine_telo_tox, calculate_rtl)
