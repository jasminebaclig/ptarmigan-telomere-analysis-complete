library(tidyverse)

raw_plate_data <- read.csv(file= './lm_info_files/rock_ptarmigan_DNA_stock_plates.xlsx.csv', header = FALSE)

lm2006 <- slice(raw_plate_data, 3:10)
lm2008 <- slice(raw_plate_data, 14:21)
lm2009 <- slice(raw_plate_data, 26:33)
lm2010 <- slice(raw_plate_data, 38:45)
lm2011 <- slice(raw_plate_data, 50:57)
lm2012 <- slice(raw_plate_data, 62:69)
lm2013 <- slice(raw_plate_data, 74:81)
lm2014 <- slice(raw_plate_data, 86:93)
lm2015 <- slice(raw_plate_data, 98:105)
lm2016 <- slice(raw_plate_data, 110:117)
lm2017 <- slice(raw_plate_data, 122:129)
lm2018 <- slice(raw_plate_data, 134:141)
lmxtra <- slice(raw_plate_data, 149:156)
lm2007 <- slice(raw_plate_data, 161:168)

plates_names <- c('lm2006', 'lm2007', 'lm2008', 'lm2009', 'lm2010', 'lm2011', 'lm2012', 'lm2013', 'lm2014', 'lm2015',
                  'lm2016', 'lm2017', 'lm2018', 'lmxtra') 

for (df_name in plates_names) {
  df <- get (df_name) %>% select(2: 13)
  assign(df_name, df)
} #takes all plate names in vector and removes column 1, which didnt have anything useful in it anyway

rm(df)

letters_a_h <- letters[1:8]

numbers <- rep(1:12, each = 8)

alphanumeric <- tibble(letters = rep(letters_a_h, length.out = length(numbers)), numbers = numbers)

#loop to pivot each 8X19 into a 96X1 to match alphanumeric
for (df_name in plates_names) {
  df <- get(df_name) %>% 
  pivot_longer(cols = everything(), cols_vary = "slowest", names_to = NULL)
  assign(df_name, df)}

for (df_name in plates_names) {
  df <- get(df_name) %>% 
    bind_cols(alphanumeric)
  assign(df_name, df)
}

all_samples <- bind_rows(lm2006, lm2007, lm2008, lm2009, lm2010, lm2011, lm2012, lm2013, lm2014, lm2015,
                         lm2016, lm2017, lm2018, lmxtra)

all_samples <- as_tibble(all_samples) %>% drop_na() 
  
all_samples <- left_join(all_samples, age_sex, by = "value")

rm(age_sex, alphanumeric, df, lm2006, lm2007, lm2008, lm2009, lm2010, lm2011, lm2012, lm2013, lm2014, lm2015, lm2016, lm2017, lm2018, lmxtra, raw_plate_data)
rm(df_name, letters_a_h, numbers, plates_names)
