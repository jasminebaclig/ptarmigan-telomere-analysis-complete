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
all_data <- rename(all_data, id = "Sample Name", gene = "Gene Name", sample_cq = "Cq Mean", cq_error = "Cq Error", effic = "indiv PCR eff")
reduced_data <- filter(all_data, !is.na(as.numeric(cq_error)), cq_error <= 0.5, !grepl("Negative", id))
avg_data <- group_by(reduced_data, run, id, gene, sample_cq, cq_error) %>% summarise_at(vars(effic), list(effic = mean))
non_single_data <- count(reduced_data, run, id, gene) %>% filter(n != 1) %>%
                   select(-n) %>%
                   left_join(avg_data, by = join_by(run, id, gene))


temp <- filter(reduced_data, is.nan(effic))


# #Goes through every plate
# for(i in 1:length(plate_names)) {
#   df <- get(plate_names[[i]])
#   
#   #Selects desired rows and columns
#   df <- filter(df, !is.na(as.numeric(df$"Cq Error")), #Selects non-excluded samples, excluding...
#                df$"Cq Error" <= 0.2, #...samples with high Cq error,...
#                !grepl("S", df$"Sample Name"), #...standard, and...
#                !grepl("N", df$"Sample Name")) #...negatives marked as positive
#   df <- select(df, "Sample Name", "Cq Mean") %>%
#         rename(id = "Sample Name", sample_cq = "Cq Mean") %>%
#         distinct() #Removes two of triplicate data
#   df$sample_cq = as.numeric(df$sample_cq)
#   
#   #Adds plate names as a column that depends on position in list
#   df$plate_name = c(plate_names[i])
#   
#   #Adds plate efficiency as a column that depends on position in list
#   df$plate_efficiency = as.numeric(c(efficiencies[i]))
#   
#   #Adds GB Cq as a column
#   df$gb_cq = c(df$sample_cq[which(df$id == "GD")]) #Finds sample_cq of GB and puts as value in new column
#   df <- filter(df, id != "GD") #Removes GD as entry in table
#   
#   #Saves modified data frame to original name
#   assign(plate_names[[i]], df)
# }
# 
# 
# 
# #Left joins corresponding TELO and TOX data
# combine_telo_tox <- function(telo_plate, tox_plate) {
#   df <- left_join(telo_plate, tox_plate, by = join_by(id)) %>%
#         select(id, cq_telo = sample_cq.x, cq_tox = sample_cq.y,
#                gb_telo = gb_cq.x, gb_tox = gb_cq.y,
#                efficiency_telo = plate_efficiency.x, efficiency_tox = plate_efficiency.y)
#   return(df)
# }
# 
# plate_1 <- combine_telo_tox(telo_1, tox_1)
# plate_2 <- combine_telo_tox(telo_2, tox_2)
# plate_3 <- combine_telo_tox(telo_3, tox_3)
# redo_1 <- combine_telo_tox(redo_1_telo, redo_1_tox)
# redo_2 <- combine_telo_tox(redo_2_telo, redo_2_tox)
# 
# #Combines all plates
# combined_data <- bind_rows(plate_1, plate_2, plate_3, redo_1, redo_2) %>%
#                  filter(!grepl("12-109", id)) %>% #Removes all 12-109 samples (see Rotating Book #1 p. 143)
#                  arrange(id)
# 
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
