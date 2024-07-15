#Performs statistical tests and creates graphs for ptarmigan RTL data
#Author: Jasmine Baclig

library(tidyverse)
library(ggbeeswarm)
library(coin)
library(car)
library(outliers)



hist(lm_data$rtl)

no_outlier_data <- filter(lm_data, id != "12-072", id != "12-036")



male_data <- filter(lm_data, sex == "male")
female_data <- filter(lm_data, sex == "female")

hist(filter(male_data, rtl < 5)$rtl)
hist(filter(female_data, rtl < 10)$rtl)

#Compares variances of RTL between sexes using F-test
var.test(rtl ~ sex, data = no_outlier_data) #F = 8.9967e-07, p < 2.2e-16 (Variances are not equal)

#Performs Welch two-sample t-test
t.test(rtl ~ sex, data = no_outlier_data, var.equal = FALSE) #t = -1.1101, p = 0.2674 (Means are not different???)

#Creates box plot of RTL values by sex each year
ggplot(no_outlier_data, aes(x = sex, y = rtl)) + geom_boxplot()

var.test(rtl ~ age, data = no_outlier_data)
t.test(rtl ~ age, data = no_outlier_data, var.equal = FALSE) #Different



y2006_data <- filter(lm_data, year == "2006")
y2007_data <- filter(lm_data, year == "2007")
y2008_data <- filter(lm_data, year == "2008")
y2009_data <- filter(lm_data, year == "2009")
y2010_data <- filter(lm_data, year == "2010")
y2011_data <- filter(lm_data, year == "2011")
y2012_data <- filter(lm_data, year == "2012")
y2013_data <- filter(lm_data, year == "2013")
y2014_data <- filter(lm_data, year == "2014")
y2015_data <- filter(lm_data, year == "2015")
y2016_data <- filter(lm_data, year == "2016")
y2017_data <- filter(lm_data, year == "2017")
y2018_data <- filter(lm_data, year == "2018")

hist(y2006_data$rtl)
hist(y2007_data$rtl)
hist(y2008_data$rtl)
hist(y2009_data$rtl)
hist(y2010_data$rtl)
hist(y2011_data$rtl)
hist(y2012_data$rtl)
hist(y2013_data$rtl)
hist(y2014_data$rtl)
hist(y2015_data$rtl)
hist(y2016_data$rtl)
hist(y2017_data$rtl)
hist(y2018_data$rtl)

#Compares variances of RTL between years using Bartlett's test
bartlett.test(rtl ~ year, no_outlier_data) #K-squared = 15617, p < 2.2e-16 (Variances are not equal)

#Performs Mood's median test
median_test(rtl ~ year, data = no_outlier_data) #chi-squared = 167.44, p < 2.2e-16 (Medians are different)

#Creates violin plot of RTL values by year
ggplot(no_outlier_data, aes(x = year, y = rtl)) + geom_boxplot()



combined_data <- left_join(lm_data, select(combined_kpm, -year), by = join_by(id == id)) %>%
                 left_join(select(combined_amanda, -year, -WeightNE, -Wing, -Head), by = join_by(id == BirdID))
write.csv(combined_data, "lm_data.csv")

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = BC.index, y = rtl)) + geom_point()
summary(lm(rtl ~ BC.index, data = combined_data)) #R-squared = 0.001656, p = 0.2125

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = TotalFat, y = rtl)) + geom_point()
summary(lm(rtl ~ TotalFat, data = combined_data)) #R-squared = 0.0005338, p = 0.4415

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = WeightNE, y = rtl)) + geom_point()
summary(lm(rtl ~ WeightNE, data = combined_data)) #R-squared = 0.003785, p = 0.03345

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = Clostridia.UCG.014.Relative.Abundance, y = rtl)) + geom_point()
summary(lm(rtl ~ Clostridia.UCG.014.Relative.Abundance, data = combined_data)) #R-squared = 0.002426, p = 0.2005

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = Shuttleworthia.Relative.Abundance, y = rtl)) + geom_point()
summary(lm(rtl ~ Shuttleworthia.Relative.Abundance, data = combined_data)) #R-squared = 0.0001117, p = 0.7837

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = Actinomyces.Relative.Abundance, y = rtl)) + geom_point()
summary(lm(rtl ~ Actinomyces.Relative.Abundance, data = combined_data)) #R-squared = 0.0001794, p = 0.7279

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = Colidextribacter.Relative.Abundance, y = rtl)) + geom_point()
summary(lm(rtl ~ Colidextribacter.Relative.Abundance, data = combined_data)) #R-squared = 8.337e-06, p = 0.9402

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = ceca.weight, y = rtl)) + geom_point()
summary(lm(rtl ~ ceca.weight, data = combined_data)) #R-squared = 0.002473, p = 0.1362

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = ceca.Gut.length, y = rtl)) + geom_point()
summary(lm(rtl ~ ceca.Gut.length, data = combined_data)) #R-squared = 0.00191, p = 0.1917

ggplot(combined_data %>% filter(id != "12-072", id != "12-036"), aes(x = Shannon.Diversity.Index, y = rtl)) + geom_point()
summary(lm(rtl ~ Shannon.Diversity.Index, data = combined_data)) #R-squared = 1.895e-05, p = 0.91



#Calculates mean values by year
create_mean_table <- function(df, col) {
  mean_table <- data.frame(year = unique(df$year)) %>%
                transform(mean = sapply(year, calculate_mean, df = df, col = col))
}

calculate_mean <- function(curr_year, df, col) {
  return(mean(as.numeric(unlist((filter(df, df$year == curr_year) %>% select(col))))))
}

mean_data <- create_mean_table(no_outlier_data, "rtl") %>% rename(rtl = mean)
mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(Standardized.BC.index)), "Standardized.BC.index") %>%
                                  rename(std_bci = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_kpm, !is.na(TotalFat)), "TotalFat") %>%
                                  rename(total_fat = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_kpm, !is.na(WeightNE)), "WeightNE") %>%
                                  rename(weight = mean), by = join_by(year))

mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(ceca.weight)), "ceca.weight") %>%
                         rename(ceca_weight = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(ceca.Gut.length)), "ceca.Gut.length") %>%
                         rename(ceca_gut_length = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(Shannon.Diversity.Index)), "Shannon.Diversity.Index") %>%
                         rename(shannon = mean), by = join_by(year))

mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(Clostridia.UCG.014.Relative.Abundance)), "Clostridia.UCG.014.Relative.Abundance") %>%
                                  rename(clostridia = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(Shuttleworthia.Relative.Abundance)), "Shuttleworthia.Relative.Abundance") %>%
                                  rename(shuttleworthia = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(Actinomyces.Relative.Abundance)), "Actinomyces.Relative.Abundance") %>%
                                  rename(actinomyces = mean), by = join_by(year))
mean_data <- left_join(mean_data, create_mean_table(filter(combined_amanda, !is.na(Colidextribacter.Relative.Abundance)), "Colidextribacter.Relative.Abundance") %>%
                                  rename(colidexteribacter = mean), by = join_by(year))




#Combines mean values of physiological data with relative telomere length data for graphing
rtl_bci <- merge(rtl_mean, bci_mean, by = "year", all = TRUE)
rtl_weight <- merge(rtl_mean, weight_mean, by = "year", all = TRUE)
rtl_fat <- merge(rtl_mean, fat_mean, by = "year", all = TRUE)
rtl_gut <- merge(rtl_mean, clostridia_mean, by = "year", all = TRUE) %>%
           merge(shuttleworthia_mean, by = "year", all = TRUE) %>%
           merge(actinomyces_mean, by = "year", all = TRUE) %>%
           merge(colidextribacter_mean, by = "year", all = TRUE)
rtl_ceca <- merge(rtl_mean, ceca_weight_mean, by = "year", all = TRUE) %>%
            merge(ceca_gut_mean, by = "year", all = TRUE)
rtl_shannon <- merge(rtl_mean, shannon_mean, by = "year", all = TRUE)
rtl_pop <- merge(rtl_mean, pop_density, by = "year", all = TRUE) %>% filter(!is.na(mean), year != 2008)



#Combines RTL and physiological data for each bird for statistical analysis
rtl_corr <- left_join(select(lm_data, id, rtl), select(combined_amanda, BirdID, Standardized.BC.index, Clostridia.UCG.014.Relative.Abundance, Shuttleworthia.Relative.Abundance, Actinomyces.Relative.Abundance, Colidextribacter.Relative.Abundance, ceca.weight, ceca.Gut.length, Shannon.Diversity.Index), by = join_by(id == BirdID)) %>%
            left_join(select(combined_kpm, id, WeightNE, TotalFat), by = join_by(id == id))

#Calculates correlation values
summary(lm(rtl ~ Standardized.BC.index, data = rtl_corr)) #R-squared = -0.006178, p = 0.4132
summary(lm(rtl ~ WeightNE, data = rtl_corr)) #R-squared = 0.01096, p = 0.202
summary(lm(rtl ~ TotalFat, data = rtl_corr)) #R-squared = 0.004121, p = 0.2684
summary(lm(rtl ~ Clostridia.UCG.014.Relative.Abundance, data = rtl_corr)) #R-squared = 0.02013, p = 0.2015
summary(lm(rtl ~ Shuttleworthia.Relative.Abundance, data = rtl_corr)) #R-squared = 0.09589, p = 0.0393
summary(lm(rtl ~ Actinomyces.Relative.Abundance, data = rtl_corr)) #R-squared = 0.03583, p = 0.142
summary(lm(rtl ~ Colidextribacter.Relative.Abundance, data = rtl_corr)) #R-squared = 0.08064, p = 0.0543
summary(lm(rtl ~ ceca.weight, data = rtl_corr)) #R-squared = 0.1202, p = 0.007861
summary(lm(rtl ~ ceca.Gut.length, data = rtl_corr)) #R-squared = 0.05568, p = 0.05438
summary(lm(rtl ~ Shannon.Diversity.Index, data = rtl_corr)) #R-squared = -0.01098, p = 0.4328

summary(lm(mean ~ rtl_pop$"density/km2", data = rtl_pop)) #R-squared = 0.09964, p = 0.2535



#Cleans up unwanted data in memory
rm(pop_density)
rm(create_mean_table, calculate_mean)
rm(rtl_mean, bci_mean, fat_mean, weight_mean, clostridia_mean, shuttleworthia_mean, actinomyces_mean, colidextribacter_mean, ceca_weight_mean, ceca_gut_mean, shannon_mean)
