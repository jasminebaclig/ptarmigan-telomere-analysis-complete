library(tidyverse)

male_juvenile <- read.csv(file= './lm_info_files/male_juvenile.csv', header = FALSE) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), cols_vary = "fastest", names_to = NULL) %>% 
  slice(1:403)

female_juvenile <- read.csv(file= './lm_info_files/female_juvenile.csv', header = FALSE) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), cols_vary = "fastest", names_to = NULL) %>% 
  slice(1:(385-6))

male_adult <- read.csv(file= './lm_info_files/male_adult.csv', header = FALSE) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), cols_vary = "fastest", names_to = NULL) %>% 
  slice(1:313)

female_adult <- read.csv(file= './lm_info_files/female_adult.csv', header = FALSE) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), cols_vary = "fastest", names_to = NULL) %>% 
  slice(1:199)


female <- c('female')
male <- c('male')
juvenile <- c('juvenile')
adult <- c('adult')

female_adult <- female_adult %>% mutate(sex = (rep(female, length.out = length(value))))
female_adult <- female_adult %>% mutate(age = (rep(adult, length.out = length(value))))

female_juvenile <- female_juvenile %>% mutate(sex = (rep(female, length.out = length(value))))
female_juvenile <- female_juvenile %>% mutate(age = (rep(juvenile, length.out = length(value))))

male_adult <- male_adult %>% mutate(sex = (rep(male, length.out = length(value))))
male_adult <- male_adult %>% mutate(age = (rep(adult, length.out = length(value))))

male_juvenile <- male_juvenile %>% mutate(sex = (rep(male, length.out = length(value))))
male_juvenile <- male_juvenile %>% mutate(age = (rep(juvenile, length.out = length(value))))

age_sex <- bind_rows(female_adult, male_adult, female_juvenile, male_juvenile, .id = NULL)

rm(female_adult, female_juvenile, male_adult, male_juvenile)
rm(adult, female, juvenile, male)
