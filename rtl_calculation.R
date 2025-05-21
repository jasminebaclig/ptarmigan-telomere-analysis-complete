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

## Cleans up data frame
all_data <- rename(all_data, id = "Sample Name", gene = "Gene Name", sample_cq = "Cq Mean", cq_error = "Cq Error", sample_effic = "indiv PCR eff")
# Removes empty wells, samples with bad Cq error, and negative control
reduced_data <- filter(all_data, !is.na(as.numeric(cq_error)), cq_error < 0.51, !grepl("Negative", id)) %>% select(-cq_error)
# Creates one entry for each unique sample with the data from the replicate runs averaged
avg_data <- group_by(reduced_data, run, id, gene, sample_cq) %>% summarise_at(vars(sample_effic), list(sample_effic = mean))
# Removes samples with only one data point (i.e., the other replicate reactions failed)
non_single_data <- count(reduced_data, run, id, gene) %>% filter(n != 1) %>%
                   select(-n) %>%
                   left_join(avg_data, by = join_by(run, id, gene)) %>%
                   mutate(sample_effic = round(sample_effic, 1))
# Selects samples with a good efficiency
good_effic_data <- filter(non_single_data, !is.nan(sample_effic), sample_effic <= 2.2, (gene == "TOX" & sample_effic >= 1.7) | (gene == "TELO" & sample_effic >= 1.4))
# Selects samples with a bad efficiency to track what samples to redo (Note: if a GB sample appears here, the whole plate should be repeated)
bad_effic_data <- filter(non_single_data, is.nan(sample_effic) | sample_effic > 2.2 | (gene == "TOX" & sample_effic < 1.7) | (gene == "TELO" & sample_effic < 1.4))


# Adds appropriate GB data to sample rows
gb_data <- filter(good_effic_data, id == "GB") %>% select(run, gene, sample_cq, sample_effic) %>% rename(gb_cq = sample_cq, gb_effic = sample_effic)
sample_data <- left_join(good_effic_data %>% filter(id != "GB"), gb_data, by = join_by(run, gene))


# Selects samples that have more than one instance where the TOX and/or TELO assay passed.
# One of these instances should be excluded in the LightCycler file. Use sample_data to know which run the duplicate samples are in.
duplicates <- count(sample_data, id, gene) %>% filter(n != 1)


# Combines TOX and TELO data for each sample
sample_id <- sample_data["id"] %>% distinct()
tox_telo_data <- data.frame()

for(i in 1:length(sample_id[[1]])) {
  tox_data <- filter(sample_data, id == sample_id[[1]][i], gene == "TOX")
  telo_data <- filter(sample_data, id == sample_id[[1]][i], gene == "TELO")
  
  if(dim(tox_data)[[1]] != 0 && dim(telo_data)[[1]] != 0) {
    new_row <- data.frame(id = sample_id[[1]][i], tox_sample_cq = as.numeric(tox_data["sample_cq"][[1]]), tox_sample_effic = as.numeric(tox_data["sample_effic"][[1]]),
                                                  tox_gb_cq = as.numeric(tox_data["gb_cq"][[1]]), tox_gb_effic = as.numeric(tox_data["gb_effic"][[1]]),
                                                  telo_sample_cq = as.numeric(telo_data["sample_cq"][[1]]), telo_sample_effic = as.numeric(telo_data["sample_effic"][[1]]),
                                                  telo_gb_cq = as.numeric(telo_data["gb_cq"][[1]]), telo_gb_effic = as.numeric(telo_data["gb_effic"][[1]]))
    tox_telo_data <- rbind(tox_telo_data, new_row)
  }
}


# Indicates which samples are not yet done given a list of all samples available
missing_data <- data.frame(setdiff(pull(age_sex, value), pull(tox_telo_data, id)))


# Calculates relative telomere length
calculate_rtl <- function(tox_sample_cq, tox_sample_effic, tox_gb_cq, tox_gb_effic, telo_sample_cq, telo_sample_effic, telo_gb_cq, telo_gb_effic) {
  telo_gb_calc <- telo_gb_effic ^ telo_gb_cq
  telo_sample_calc <- telo_sample_effic ^ telo_sample_cq
  telo_calc <- telo_gb_calc / telo_sample_calc
  
  tox_gb_calc <- tox_gb_effic ^ tox_gb_cq
  tox_sample_calc <- tox_sample_effic ^ tox_sample_cq
  tox_calc <- tox_gb_calc / tox_sample_calc
  
  return(telo_calc / tox_calc)
}

calculated_data <- transform(tox_telo_data, rtl = calculate_rtl(tox_sample_cq, tox_sample_effic, tox_gb_cq, tox_gb_effic,
                                                                telo_sample_cq, telo_sample_effic, telo_gb_cq, telo_gb_effic))

#Cleans up unwanted data in memory
rm(age_sex, all_data, avg_data, df_effic, df_txt, gb_data, good_effic_data, new_row, non_single_data, reduced_data, sample_data, sample_id, telo_data, tox_data, tox_telo_data)
rm(effic_files, effic_path, i, txt_files, txt_path, calculate_rtl)
