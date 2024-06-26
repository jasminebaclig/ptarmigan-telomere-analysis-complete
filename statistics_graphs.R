#Performs statistical tests and creates graphs for ptarmigan RTL data
#Author: Jasmine Baclig

library(tidyverse)
library(ggbeeswarm)
library(coin)
library(car)
library(outliers)



#Performs Grubbs' Test on potential outliers and removes them from data frame
grubbs.test(lm_data$rtl) #G = 34.328657, p < 2.2e-16 (12-072 is an outlier)
no_outlier_data <- filter(lm_data, id != "12-072")
grubbs.test(no_outlier_data$rtl) #G = 3.4509e+01, p < 2.2e-16 (12-036 is an outlier)
no_outlier_data <- filter(no_outlier_data, id != "12-036")
grubbs.test(no_outlier_data$rtl) #G = 28.04866, p < 2.2e-16 (09-002 is an outlier)
no_outlier_data <- filter(no_outlier_data, id != "09-002")
grubbs.test(no_outlier_data$rtl) #G = 22.52881, p < 2.2e-16 (09-010 is an outlier)
no_outlier_data <- filter(no_outlier_data, id != "09-010")
grubbs.test(no_outlier_data$rtl) #G = 13.49175, p < 2.2e-16 (08-059 is an outlier)
no_outlier_data <- filter(no_outlier_data, id != "08-059")
grubbs.test(no_outlier_data$rtl) #G = 14.13032, p < 2.2e-16 (09-072 is an outlier)
no_outlier_data <- filter(no_outlier_data, id != "09-072")
grubbs.test(no_outlier_data$rtl) #G = 10.38871, p < 2.2e-16 (09-016 is an outlier)
no_outlier_data <- filter(no_outlier_data, id != "09-016")
grubbs.test(no_outlier_data$rtl) #G = 10.38871, p < 2.2e-16 (09-016 is an outlier)

male_data <- filter(lm_data, sex == "male")
female_data <- filter(lm_data, sex == "female")

grubbs.test(male_data$rtl)



#Compares variances of RTL between sexes using F-test
var.test(rtl ~ sex, data = lm_data) #F = 8.9967e-07, p < 2.2e-16 (Variances are not equal)

#Performs Welch two-sample t-test
t.test(rtl ~ sex, data = lm_data, var.equal = FALSE) #t = -1.1101, p = 0.2674 (Means are not different???)

#Creates violin plot of RTL values by sex each year
rtl_by_sex <- ggplot(lm_data, aes(x = sex, y = rtl)) + 
              geom_violin(aes(fill = sex)) +
              geom_beeswarm(size = 2, cex = 2, alpha = 0.7) +
              guides(fill = FALSE) +
              scale_fill_manual(values = c("#796C9D", "#F5AD29")) +
              labs(title = "Relative Telomere Lengths Are Higher Overall\nin Female Ptarmigans", x = NULL, y = "Relative Telomere Length") +
              theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 25), axis.title = element_text(size = 15), axis.text = element_text(size = 12)) +
              scale_x_discrete(labels = c("Female", "Male")) +
              scale_y_continuous(breaks = seq(0, 4, 1), minor_breaks = seq(0, 4, 0.5))
rtl_by_sex



#Compares variances of RTL between years using Bartlett's test
bartlett.test(rtl ~ year, lm_data) #K-squared = 20.208, p = 0.005137 (Variances are not equal)

#Performs Mood's median test
median_test(rtl ~ year, data = lm_data) #chi-squared = 26.548, p = 0.0004018 (Medians are different)

#Creates violin plot of RTL values by year
rtl_by_year <- ggplot(lm_data, aes(x = year, y = rtl)) +
               geom_violin(width = 1.25, aes(fill = year)) +
               geom_beeswarm(size = 1.5, cex = 1.5, alpha = 0.7) +
               guides(fill = "none") +
               scale_fill_manual(values = c("#796C9D", "#8B758C", "#9C7F7C", "#AE886B", "#C0915B", "#D29A4A", "#E3A43A", "#F5AD29")) +
               labs(title = "Relative Telomere Lengths Are Different\nThroughout 2006 to 2016", x = NULL, y = "Relative Telomere Length") +
               theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 25), axis.title = element_text(size = 15), axis.text = element_text(size = 12)) +
               scale_y_continuous(breaks = seq(0, 4, 1), minor_breaks = seq(0, 4, 0.5))
rtl_by_year



#Calculates mean values by year
create_mean_table <- function(df, col) {
  mean_table <- data.frame(year = unique(df$year)) %>%
                transform(mean = sapply(year, calculate_mean, df = df, col = col))
}

calculate_mean <- function(curr_year, df, col) {
  return(mean(as.numeric(unlist((filter(df, df$year == curr_year) %>% select(col))))))
}

rtl_mean <- create_mean_table(lm_data, "rtl")

bci_mean <- create_mean_table(filter(combined_amanda, !is.na(Standardized.BC.index)), "Standardized.BC.index")
fat_mean <- create_mean_table(filter(combined_kpm, !is.na(TotalFat)), "TotalFat")
weight_mean <- create_mean_table(filter(combined_kpm, !is.na(WeightNE)), "WeightNE")

clostridia_mean <- create_mean_table(filter(combined_amanda, !is.na(Clostridia.UCG.014.Relative.Abundance)), "Clostridia.UCG.014.Relative.Abundance")
shuttleworthia_mean <- create_mean_table(filter(combined_amanda, !is.na(Shuttleworthia.Relative.Abundance)), "Shuttleworthia.Relative.Abundance")
actinomyces_mean <- create_mean_table(filter(combined_amanda, !is.na(Actinomyces.Relative.Abundance)), "Actinomyces.Relative.Abundance")
colidextribacter_mean <- create_mean_table(filter(combined_amanda, !is.na(Colidextribacter.Relative.Abundance)), "Colidextribacter.Relative.Abundance")

ceca_weight_mean <- create_mean_table(filter(combined_amanda, !is.na(ceca.weight)), "ceca.weight")
ceca_gut_mean <- create_mean_table(filter(combined_amanda, !is.na(ceca.Gut.length)), "ceca.Gut.length")

shannon_mean <- create_mean_table(filter(combined_amanda, !is.na(Shannon.Diversity.Index)), "Shannon.Diversity.Index")

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
