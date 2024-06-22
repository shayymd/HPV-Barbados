
library(magrittr)
library(tidyverse)
library(readxl)
library(binom)
library(epiR)

#' Read HPV dataset 
data = readxl::read_xlsx("raw_hpv.xlsx",
                         col_names = T, na = c("", "9999")) #Mark 9999 as a missing value
readxl::read
#' Create tibble with desired age groups, 18-25, 25-35, 35+
age_groups = tibble(
  age_low = c(15, 25, 35),
  age_high = c(24, 34, 75)) %>%
  mutate(name = paste0(age_low, "-", age_high)) %>%
  rows_update(tibble(age_low=35, name="35+"))

#' Edit dataset to automatically allocate age groups (in new variable age_group)
data = data %>%
  mutate(age_group = NA_character_)
for(i in 1:nrow(age_groups)){
  data = data %>% rows_update(
    data %>%
      subset(age >= age_groups[[i, "age_low"]] & age <= age_groups[[i, "age_high"]]) %>%
      select(id) %>%
      mutate(age_group = age_groups[[i, "name"]]))
}

#' Summarise the data, only where hpv_16 or hpv_18 is 1
summary_data = data %>%
  select(id, age_group, hpv_16, hpv_18) %>%
  rowwise(id) %>%
  mutate(hpv = (hpv_16+hpv_18) > 0) %>%
  group_by(age_group) %>%
  summarise(total=sum(hpv), n=n())

#' Plot data by age, also calculate binomial confidence intervals around values
summary_data %>%
  cbind(binom::binom.confint(x = summary_data$total, summary_data$n, method="exact") %>%
          as_tibble() %>%
          select(mean, lower, upper)) %>%
  subset(!is.na(age_group)) %>%
  ggplot(aes(x=age_group, y=mean, ymin=lower, ymax=upper))+
  geom_errorbar()+
  geom_point()+
  scale_y_continuous(labels=scales::percent, limits = c(0, NA))+
  theme_minimal()+
  labs(x="Age group", y="HPV16/18 carrier", title="HPV 16/18 carriers by age")

# Drop missing values and rename groups to create graph for sexual partners by age group (all 5 sexual partner groups)
partner_list = data %>%
  select(partner_year, age_group) %>% 
  drop_na() %>% 
  mutate(partner_year = recode(partner_year, 
                               '0' = "0",
                               '1' = "1",
                               '2' = "2-3",
                               '3' = "4-5",
                               '4' = "6-10" ))

partner_list2 = data %>%
  select(partner_year, age_group) %>% 
  drop_na() %>% 
  mutate(partner_year_cat = recode(partner_year, 
                               '0' = "Low",
                               '1' = "Low",
                               '2' = "High",
                               '3' = "High",
                               '4' = "High")) %>%
  mutate(partner_year_med = recode(partner_year, 
                               '0' = "0",
                               '1' = "1",
                               '2' = "2.5",
                               '3' = "4.5",
                               '4' = "8")) %>%
  mutate(partner_year_med = as.numeric(partner_year_med)) %>%
    group_by(partner_year_cat) %>%
    summarise(partner_year_mean = mean(partner_year_med)) %>% saveRDS("./partnerships_year.RDS")

#' graph for sexual partners by age group for all 5 sexual partner groups
ggplot(data = partner_list) +
  geom_bar(mapping = aes(x = partner_year, y = stat(prop), group = 1)) +
  facet_wrap(~ age_group, nrow = 1,) +
  scale_y_continuous(labels=scales::percent) +
  theme_minimal() +
  labs(x="Number of sexual partners", y="Proportion", title="Sexual partners by age group")

#' Regroup sexual partners into Low and High with graph based on  0-1=low, 2-10=high
#' Consider (if time left) exploring a situation where high = 4+ sexual partners 
partner_groups = data %>%
  select(partner_year, age_group) %>% 
  drop_na() %>% 
  mutate(partner_year = recode(partner_year, 
                               '0' = "Low",
                               '1' = "Low",
                               '2' = "High",
                               '3' = "High",
                               '4' = "High" ) %>% fct_infreq())
  
#' graph for sexual partners by age group and sexual groups 
ggplot(data = partner_groups) +
  geom_bar(mapping = aes(x = partner_year, y = stat(prop), group = 1)) +
  facet_wrap(~ age_group, nrow = 1,) +
  scale_y_continuous(labels=scales::percent) +
  theme_minimal() +
  labs(x="Risk group", y="Proportion", title="Sexual partner risk group by age group")

#' crosstabulation of partner_groups and age_group
with(partner_groups, table(age_group, partner_year))

#' percentages of partner_group by age_group
#' Limitation: we only have data available for girls. We will assume that boys have similar rates
#' Make this assumption explicit in methods section
#' when 14yo girls age into the 15-24 yo age-group, a proportion 0.339 will be in the high-risk group,
#'   and (1-0.339) in the low-risk group. When 15-24yo age into the 25-34 yo group, a proportion (0.252/0.339) remains
#'   in the high-risk group, while (1 - 0.252/0.339) will move to the low-risk group
with(partner_groups, prop.table(table(age_group, partner_year), margin = 1))

# prevalence of HPV by age group, risk group
prev_list = data %>% 
  select(id, age_group, hpv_16, hpv_18, partner_year) %>%
  rowwise(id) %>%
  mutate(hpv = (hpv_16+hpv_18) > 0) %>%
  mutate(hpv = 1*hpv) %>% 
  mutate(hpv = recode(hpv, '0' = "No HPV", '1' = "HPV")) %>% 
  group_by(age_group) %>%
  drop_na() %>% 
  mutate(partner_year = recode(partner_year, 
                               '0' = "Low",
                               '1' = "Low",
                               '2' = "High",
                               '3' = "High",
                               '4' = "High" ) %>% fct_infreq()) 

#' Prevalence by age and risk-group
prev_list %>%
  group_by(age_group, partner_year) %>%
  summarise(n=n(), total=sum(hpv == "HPV")) %>%
  (function(x){
    x %>%
      cbind(x %>%
              ungroup %>%
              select(total, n) %>%
              as.matrix %>%
              epi.conf(ctype="prevalence", method="exact") %>%
              as_tibble() %>%
              select(est, lower, upper)) %>%
      return
  })() %>% saveRDS("./prevalence.RDS")

# age group prevalence
with(prev_list, prop.table(table(age_group, hpv), margin = 1))
with(prev_list, table(age_group, hpv))
6+47 # n=53
prev_age_1 = 6; agepop = 53 
  age_p_1 = as.matrix(cbind(prev_age_1, agepop))
(epi.conf(matrix(c(6, 53), ncol=2), ctype = "prevalence", method = "exact", N = Inf, design = 1, 
           conf.level = 0.95) * 100)
6+113 # n=119
prev_age_2 = 6; agepop = 119 
age_p_2 = as.matrix(cbind(prev_age_2, agepop))
epi.conf(age_p_2, ctype = "prevalence", method = "exact", N = 426, design = 1, 
         conf.level = 0.95) * 100
7+195 # n=202
prev_age_3 = 7; agepop = 202 
age_p_3 = as.matrix(cbind(prev_age_3, agepop))
epi.conf(age_p_3, ctype = "prevalence", method = "exact", N = 426, design = 1, 
         conf.level = 0.95) * 100

age_groups_prevalence = tibble(
  prevalence = c(11.32075, 5.042017, 2.970297), 
  lower = c(4.269639, 1.872551, 1.09769),
    upper = c(23.02899, 10.6515, 6.352666))

# risk group prevalence
with(prev_list, prop.table(table(partner_year, hpv), margin = 1))  
with(prev_list, table(partner_year, hpv))
12+288 # n=300
prev_risk_l = 12; agepop = 300 
risk_p_l = as.matrix(cbind(prev_age, agepop))
epi.conf(risk_p_l, ctype = "prevalence", method = "exact", N = 426, design = 1, 
         conf.level = 0.95) * 100
7+67 # n=74
prev_risk_h = 7; agepop = 74 
risk_p_h = as.matrix(cbind(prev_age, agepop))
epi.conf(risk_p_h, ctype = "prevalence", method = "exact", N = 426, design = 1, 
         conf.level = 0.95) * 100

sex_risk_prevalence = tibble(
  prevalence = c(2, 8.108108), 
  lower = c(0.7374162, 3.033757),
  upper = c(4.302108, 16.81712))


