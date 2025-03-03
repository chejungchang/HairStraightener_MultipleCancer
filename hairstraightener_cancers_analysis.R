# Hair Straightener and Cancers (analysis) -- 

library(tidyverse)
library(survival)
library(broom)

# Load data ----
dat <- read_csv("clean_data/clean_data.csv")

# Check the cancer numbers ----
order <- c("event_HdNk2", "event_CrCNS2", "event_Mel2", "event_Thy2", "event_Lung2", 
           "event_LymHem2", "event_NHLym2", "event_HLym2", "event_MMyel2", "event_Leuk2", 
           "event_Dig2", "event_Pan2","event_CR2", 
           "event_Uri2", "event_Ki2", "event_Bldr2")

dat %>% mutate(FT = eof - age) %>%
  summarise(medianFT = median(FT, na.rm = TRUE)) # medican follow up year 13.1 

# Case numbers 
case <- dat %>% 
  select(ex_12mo_ever, starts_with("event_") & ends_with("2")) %>%
  group_by(ex_12mo_ever) %>%
  summarise(
    across(
      starts_with("event_") & ends_with("2"),
      ~ sum(.x, na.rm = TRUE)  # or sum(.x) if no NA issues
    )
  ) %>% 
  select(all_of(order)) %>%
  t()


case2 <- dat %>% 
  select(ex_12mo_freq, starts_with("event_") & ends_with("2")) %>%
  group_by(ex_12mo_freq) %>%
  summarise(
    across(
      starts_with("event_") & ends_with("2"),
      ~ sum(.x, na.rm = TRUE)  # or sum(.x) if no NA issues
    )
  ) %>% 
  select(all_of(order)) %>%
  t()


# Table 2 ever never ---- 

mycox <- function(data, Event_var) {

  # number of cases and non cases
  cancer_name <- data %>% select(Event_var) %>% colnames() 
  all_case <- data %>% filter(.[[Event_var]]==1) %>% nrow()
  all_noncase <- data %>% filter(.[[Event_var]]==0) %>% nrow()
  ever_noncase <- data %>% filter(ex_12mo_ever==1, .[[Event_var]]==0) %>% nrow()  
  ever_noncase_per <- ever_noncase/all_noncase * 100
  
  ever_case <- data %>% filter(ex_12mo_ever==1, .[[Event_var]]==1) %>% nrow()  
  ever_case_per <- ever_case/all_case * 100
  
  # cox model 1 - age 
  mod1 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_ever) + factor(race1), ties = 'breslow', data = data)
  mod1.res <- mod1 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(estimate, conf.low, conf.high)
  
  # cox model 2 - age race edu smoking
  mod2 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_ever) + factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = data)
  mod2.res <- mod2 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(estimate, conf.low, conf.high)
  
  # cox model 3 - age race edu smoking BMI
  mod3 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_ever) + factor(race1) + factor(educ) + factor(smoking1) + BMI_ex, ties = 'breslow', data = data)
  mod3.res <- mod3 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(estimate, conf.low, conf.high)
  
  # Column names for consistency
  col_names <- c("cancer_type", 
                 "all_noncase", "ever_noncase_per",
                 "all_case", "ever_case_per", 
                 "mod1_est", "mod1_low", "mod1_high", 
                 "mod2_est", "mod2_low", "mod2_high", 
                 "mod3_est", "mod3_low", "mod3_high")
  
  # For ever group (combine model results)
  res <- cbind(cancer_name, all_noncase, ever_noncase_per, all_case, ever_case_per, mod1.res, mod2.res, mod3.res)
  colnames(res) <- col_names
  
  return(res)
}



results_list <- list(
  Mel <- mycox(dat, "event_Mel2"),
  Thy <- mycox(dat, "event_Thy2"),
  Lung <- mycox(dat, "event_Lung2"),
  LymHem <- mycox(dat, "event_LymHem2"),
  NHLym <- mycox(dat, "event_NHLym2"),
  Leuk <- mycox(dat, "event_Leuk2"),
  Dig <- mycox(dat, "event_Dig2"),
  Pan <- mycox(dat, "event_Pan2"),
  CR <- mycox(dat, "event_CR2"),
  Uri <- mycox(dat, "event_Uri2"),
  Ki <- mycox(dat, "event_Ki2")
)


combined_results <- do.call(rbind, results_list)

plot1 <- combined_results 

combined_results <- combined_results %>% mutate(
                            noncase = format(round(all_noncase, digits = 0), nsmall = 0), 
                            ever_noncase_per = format(round(ever_noncase_per, digits = 1), nsmall = 1), 
                            number1 = str_c(noncase, " (", ever_noncase_per, "%)"),
  
                            case = format(round(all_case, digits = 0), nsmall = 0), 
                            ever_case_per = format(round(ever_case_per, digits = 1), nsmall = 1), 
                            number2 = str_c(case, " (", ever_case_per, "%)"),
                            
                            estimate1 = format(round(mod1_est, digits = 2), nsmall = 1), 
                            lower1 = format(round(mod1_low, digits = 2), nsmall = 1), 
                            upper1 = format(round(mod1_high, digits = 2), nsmall = 1),
                            est1 = str_c(estimate1, " (", lower1, ", ", upper1, ")"),
                            
                            estimate2 = format(round(mod2_est, digits = 2), nsmall = 1), 
                            lower2 = format(round(mod2_low, digits = 2), nsmall = 1), 
                            upper2 = format(round(mod2_high, digits = 2), nsmall = 1),
                            est2 = str_c(estimate2, " (", lower2, ", ", upper2, ")"),
                            
                            estimate3 = format(round(mod3_est, digits = 2), nsmall = 1), 
                            lower3 = format(round(mod3_low, digits = 2), nsmall = 1), 
                            upper3 = format(round(mod3_high, digits = 2), nsmall = 1),
                            est3 = str_c(estimate3, " (", lower3, ", ", upper3, ")")) |>
  select(cancer_type, number1, number2, est1, est2, est3) 



# Table 3 frequency ---- 

mycoxF <- function(data, Event_var) {
  
  # number of cases and non cases
  cancer_name <- data %>% select(Event_var) %>% colnames() 
  
  case <- data %>% filter(.[[Event_var]]==1) %>% nrow()
  f1_case <- data %>% filter(ex_12mo_freq==1, .[[Event_var]]==1) %>% nrow()  
  f1_case_per <- f1_case/all_case * 100
  f2_case <- data %>% filter(ex_12mo_freq==2, .[[Event_var]]==1) %>% nrow()  
  f2_case_per <- f2_case/all_case * 100
  
  # cox model 1 - age
  mod1 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_freq) , ties = 'breslow', data = data)
  mod1.res <- mod1 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1:2) %>% select(estimate, conf.low, conf.high)
  
  # cox model 2 - age race edu smoking
  mod2 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_freq) + factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = data)
  mod2.res <- mod2 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1:2) %>% select(estimate, conf.low, conf.high)
  
  # cox model 3 - age race edu smoking BMI
  mod3 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_freq) + factor(race1) + factor(educ) + factor(smoking1) + BMI, ties = 'breslow', data = data)
  mod3.res <- mod3 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1:2) %>% select(estimate, conf.low, conf.high)
  
  # Column names for consistency
  col_names <- c("cancer_type", "all_case", "ever_case_per", 
                 "mod1_est", "mod1_low", "mod1_high", 
                 "mod2_est", "mod2_low", "mod2_high", 
                 "mod3_est", "mod3_low", "mod3_high")
  
  # For never group (use 1 for all model values as a placeholder)
  case <- cbind(rbind(cancer_name, cancer_name), rbind(case, case), rbind(f1_case_per, f2_case_per))
  modall.res <- cbind(mod1.res, mod2.res, mod3.res)
  
  # Combine the two groups
  res <- cbind(case, modall.res)
  colnames(res) <- col_names 
  
  return(res)
}


results_list2 <- list(
  Mel <- mycoxF(dat, "event_Mel2"),
  Thy <- mycoxF(dat, "event_Thy2"),
  Lung <- mycoxF(dat, "event_Lung2"),
  LymHem <- mycoxF(dat, "event_LymHem2"),
  NHLym <- mycoxF(dat, "event_NHLym2"),
  Leuk <- mycoxF(dat, "event_Leuk2"),
  Dig <- mycoxF(dat, "event_Dig2"),
  Pan <- mycoxF(dat, "event_Pan2"),
  CR <- mycoxF(dat, "event_CR2"),
  Uri <- mycoxF(dat, "event_Uri2"),
  Ki <- mycoxF(dat, "event_Ki2")
)

combined_results2 <- do.call(rbind, results_list2) %>% 
  mutate(all_case = as.numeric(all_case),
         ever_case_per = as.numeric(ever_case_per))
plot2 <- combined_results2 

combined_results2 <- combined_results2 %>% mutate(
  case = format(round(all_case, digits = 0), nsmall = 0), 
  ever_case_per = format(round(ever_case_per, digits = 1), nsmall = 1), 
  number1 = str_c(case),
  number2 = str_c(ever_case_per, "%"),
  
  estimate1 = format(round(mod1_est, digits = 2), nsmall = 1), 
  lower1 = format(round(mod1_low, digits = 2), nsmall = 1), 
  upper1 = format(round(mod1_high, digits = 2), nsmall = 1),
  est1 = str_c(estimate1, " (", lower1, ", ", upper1, ")"),
  
  estimate2 = format(round(mod2_est, digits = 2), nsmall = 1), 
  lower2 = format(round(mod2_low, digits = 2), nsmall = 1), 
  upper2 = format(round(mod2_high, digits = 2), nsmall = 1),
  est2 = str_c(estimate2, " (", lower2, ", ", upper2, ")"),
  
  estimate3 = format(round(mod3_est, digits = 2), nsmall = 1), 
  lower3 = format(round(mod3_low, digits = 2), nsmall = 1), 
  upper3 = format(round(mod3_high, digits = 2), nsmall = 1),
  est3 = str_c(estimate3, " (", lower3, ", ", upper3, ")")) |>
  select(cancer_type, number1, number2, est1, est2, est3) 

# P for trend -----
 
mycox_ptrend<- function(data, Event_var) {
  
  cancer_name <- data %>% select(Event_var) %>% colnames() 
  
  # cox model 2 - age race edu smoking
  mod2 <- coxph(Surv(age, eof, data[[Event_var]]) ~ as.numeric(ex_12mo_freq) + factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = data)
  mod2.res <- mod2 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(p.value)
  
  # Column names for consistency
  col_names <- c("cancer_type", "ptrend")
  
  # Combine the two groups
  res <- cbind(cancer_name, mod2.res)
  colnames(res) <- col_names 
  
  return(res)
}

coxph(Surv(age, eof, event_Thy2) ~ as.numeric(ex_12mo_freq) + factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = dat)


results_list3 <- list(
  Mel <- mycox_ptrend(dat, "event_Mel2"),
  Thy <- mycox_ptrend(dat, "event_Thy2"),
  Lung <- mycox_ptrend(dat, "event_Lung2"),
  LymHem <- mycox_ptrend(dat, "event_LymHem2"),
  NHLym <- mycox_ptrend(dat, "event_NHLym2"),
  Leuk <- mycox_ptrend(dat, "event_Leuk2"),
  Dig <- mycox_ptrend(dat, "event_Dig2"),
  Pan <- mycox_ptrend(dat, "event_Pan2"),
  CR <- mycox_ptrend(dat, "event_CR2"),
  Uri <- mycox_ptrend(dat, "event_Uri2"),
  Ki <- mycox_ptrend(dat, "event_Ki2")
)

combined_results3 <- do.call(rbind, results_list3)


# Table 4 Black and non-Hispanic White ---- 

case3_1 <- dat %>% 
  filter(race1 == 1) %>%
  select(ex_12mo_ever, starts_with("event_") & ends_with("2")) %>%
  group_by(ex_12mo_ever) %>%
  summarise(
    across(
      starts_with("event_") & ends_with("2"),
      ~ sum(.x, na.rm = TRUE)  # or sum(.x) if no NA issues
    )
  ) %>% 
  select(all_of(order)) %>%
  t() %>% as.data.frame() %>%
  mutate(sum = V1 + V2, 
         per = V2/sum*100, 
         case = format(round(sum, digits = 0), nsmall = 0), 
         ever_case_per = format(round(per, digits = 1), nsmall = 1), 
         number = str_c(case, " (", ever_case_per, "%)")
  ) %>% select(V1,V2, sum, ever_case_per, number)

case3_2 <- dat %>% 
  filter(race1 == 2) %>%
  select(ex_12mo_ever, starts_with("event_") & ends_with("2")) %>%
  group_by(ex_12mo_ever) %>%
  summarise(
    across(
      starts_with("event_") & ends_with("2"),
      ~ sum(.x, na.rm = TRUE)  # or sum(.x) if no NA issues
    )
  ) %>% 
  select(all_of(order)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sum = V1 + V2, 
         per = V2/sum*100, 
         case = format(round(sum, digits = 0), nsmall = 0), 
         ever_case_per = format(round(per, digits = 1), nsmall = 1), 
         number = str_c(case, " (", ever_case_per, "%)")
  ) %>% select(V1, V2, sum, ever_case_per, number)

library(emmeans)

mycoxBW <- function(data, Event_var) {
  
  datain <- data[data$race1 %in% c(1,2),]
  # number of cases and non cases
  cancer_name <- data %>% select(Event_var) %>% colnames() 
  
  mod <- coxph(Surv(age, eof, datain[[Event_var]]) ~ factor(ex_12mo_ever) * factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = datain)
  
  em_cat <- emmeans(mod, ~ factor(ex_12mo_ever) | factor(race1), type = "response") %>% 
    pairs(reverse = TRUE, infer = TRUE) %>% as.data.frame() %>% 
    select(ratio, asymp.LCL, asymp.UCL)
  
  # Column names for consistency
  col_names <- c("cancer_type", "est", "low", "high")
  
  # For ever group (combine model results)
  res <- cbind(rbind(cancer_name, cancer_name), em_cat)
  colnames(res) <- col_names
  
  return(res)
}


results_list4 <- list(
  Mel <- mycoxBW(dat, "event_Mel2"),
  Thy <- mycoxBW(dat, "event_Thy2"),
  Lung <- mycoxBW(dat, "event_Lung2"),
  LymHem <- mycoxBW(dat, "event_LymHem2"),
  NHLym <- mycoxBW(dat, "event_NHLym2"),
  Leuk <- mycoxBW(dat, "event_Leuk2"),
  Dig <- mycoxBW(dat, "event_Dig2"),
  Pan <- mycoxBW(dat, "event_Pan2"),
  CR <- mycoxBW(dat, "event_CR2"),
  Uri <- mycoxBW(dat, "event_Uri2"),
  Ki <- mycoxBW(dat, "event_Ki2")
)

combined_results4 <- do.call(rbind, results_list4)
combined_results4 <- combined_results4 %>% 
  mutate(
  estimate = format(round(est, digits = 2), nsmall = 1), 
  lower = format(round(low, digits = 2), nsmall = 1), 
  upper = format(round(high, digits = 2), nsmall = 1),
  est = str_c(estimate, " (", lower, ", ", upper, ")")) |>
  select(cancer_type, est) 

# p interaction ----

coxph(Surv(age, eof, event_Lung2) ~ factor(ex_12mo_ever) * factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = dat[dat$race1 %in% c(1,2),]) %>% summary()
coxph(Surv(age, eof, event_LymHem2) ~ factor(ex_12mo_ever) * factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = dat[dat$race1 %in% c(1,2),]) %>% summary()
coxph(Surv(age, eof, event_Dig2) ~ factor(ex_12mo_ever) * factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = dat[dat$race1 %in% c(1,2),]) %>% summary()
coxph(Surv(age, eof, event_Uri2) ~ factor(ex_12mo_ever) * factor(race1) + factor(educ) + factor(smoking1), ties = 'breslow', data = dat[dat$race1 %in% c(1,2),]) %>% summary()


# BMI missing numbers -----
case3 <- dat %>% 
  filter(!is.na(BMI)) %>%
  select(ex_12mo_ever, starts_with("event_") & ends_with("2")) %>%
  group_by(ex_12mo_ever) %>%
  summarise(
    across(
      starts_with("event_") & ends_with("2"),
      ~ sum(.x, na.rm = TRUE)  # or sum(.x) if no NA issues
    )
  ) %>% 
  select(all_of(order)) %>%
  t() %>% as.data.frame() %>%
  mutate(sum = V1 + V2, 
         per = V2/sum*100,
         case = format(round(sum, digits = 0), nsmall = 0), 
         ever_case_per = format(round(per, digits = 1), nsmall = 1), 
         number = str_c(case, " (", ever_case_per, "%)"))

dat %>% filter(!is.na(BMI)) %>% nrow() #46157


# Which variable changes the estimate the most ----
library(huxtable)

mod1 <- coxph(Surv(age, eof, event_NHLym2) ~ factor(ex_12mo_ever) , ties = 'breslow', data = dat)
mod2 <- coxph(Surv(age, eof, event_NHLym2) ~ factor(ex_12mo_ever) + factor(race1) , ties = 'breslow', data = dat)
mod3 <- coxph(Surv(age, eof, event_NHLym2) ~ factor(ex_12mo_ever) + factor(educ) , ties = 'breslow', data = dat)
mod4 <- coxph(Surv(age, eof, event_NHLym2) ~ factor(ex_12mo_ever) + factor(smoking1) , ties = 'breslow', data = dat)

huxreg("Age" = mod1, "+ Race" = mod2, "+ Education" = mod3, "+ Smoking" = mod4,  tidy_args = list(exponentiate = TRUE))



mod1 <- coxph(Surv(age, eof, event_Thy2) ~ factor(ex_12mo_ever) , ties = 'breslow', data = dat)
mod2 <- coxph(Surv(age, eof, event_Thy2) ~ factor(ex_12mo_ever) + factor(race1) , ties = 'breslow', data = dat)
mod3 <- coxph(Surv(age, eof, event_Thy2) ~ factor(ex_12mo_ever) + factor(educ) , ties = 'breslow', data = dat)
mod4 <- coxph(Surv(age, eof, event_Thy2) ~ factor(ex_12mo_ever) + factor(smoking1) , ties = 'breslow', data = dat)

huxreg("Age" = mod1, "+ Race" = mod2, "+ Education" = mod3, "+ Smoking" = mod4,  tidy_args = list(exponentiate = TRUE))


mod1 <- coxph(Surv(age, eof, event_Pan2) ~ factor(ex_12mo_ever) , ties = 'breslow', data = dat)
mod2 <- coxph(Surv(age, eof, event_Pan2) ~ factor(ex_12mo_ever) + factor(race1) , ties = 'breslow', data = dat)
mod3 <- coxph(Surv(age, eof, event_Pan2) ~ factor(ex_12mo_ever) + factor(educ) , ties = 'breslow', data = dat)
mod4 <- coxph(Surv(age, eof, event_Pan2) ~ factor(ex_12mo_ever) + factor(smoking1) , ties = 'breslow', data = dat)

huxreg("Age" = mod1, "+ Race" = mod2, "+ Education" = mod3, "+ Smoking" = mod4,  tidy_args = list(exponentiate = TRUE))


mod1 <- coxph(Surv(age, eof, event_CR2) ~ factor(ex_12mo_ever) , ties = 'breslow', data = dat)
mod2 <- coxph(Surv(age, eof, event_CR2) ~ factor(ex_12mo_ever) + factor(race1) , ties = 'breslow', data = dat)
mod3 <- coxph(Surv(age, eof, event_CR2) ~ factor(ex_12mo_ever) + factor(educ) , ties = 'breslow', data = dat)
mod4 <- coxph(Surv(age, eof, event_CR2) ~ factor(ex_12mo_ever) + factor(smoking1) , ties = 'breslow', data = dat)

huxreg("Age" = mod1, "+ Race" = mod2, "+ Education" = mod3, "+ Smoking" = mod4,  tidy_args = list(exponentiate = TRUE))


# Subset to white women ----
# Table 2 ever never ---- 

mycox_w <- function(data, Event_var) {
  
  # number of cases and non cases
  cancer_name <- data %>% select(Event_var) %>% colnames() 
  all_case <- data %>% filter(.[[Event_var]]==1) %>% nrow()
  ever_case <- data %>% filter(ex_12mo_ever==1, .[[Event_var]]==1) %>% nrow()  
  ever_case_per <- ever_case/all_case * 100
  
  # cox model 1 - age
  mod1 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_ever) , ties = 'breslow', data = data)
  mod1.res <- mod1 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(estimate, conf.low, conf.high)
  
  # cox model 2 - age edu smoking
  mod2 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_ever) +  factor(educ) + factor(smoking1), ties = 'breslow', data = data)
  mod2.res <- mod2 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(estimate, conf.low, conf.high)
  
  # cox model 3 - age edu smoking BMI
  mod3 <- coxph(Surv(age, eof, data[[Event_var]]) ~ 
                  factor(ex_12mo_ever) + factor(educ) + factor(smoking1) + BMI, ties = 'breslow', data = data)
  mod3.res <- mod3 %>% tidy(., exponentiate = TRUE, conf.int = TRUE) %>% slice(1) %>% select(estimate, conf.low, conf.high)
  
  # Column names for consistency
  col_names <- c("cancer_type", "all_case", "ever_case_per", "mod1_est", "mod1_low", "mod1_high", 
                 "mod2_est", "mod2_low", "mod2_high", 
                 "mod3_est", "mod3_low", "mod3_high")
  
  # For ever group (combine model results)
  res <- cbind(cancer_name, all_case, ever_case_per, mod1.res, mod2.res, mod3.res)
  colnames(res) <- col_names
  
  return(res)
}


dat_w <- dat %>% filter(race1 == 1)

results_list <- list(
  Mel <- mycox_w(dat_w, "event_Mel2"),
  Thy <- mycox_w(dat_w, "event_Thy2"),
  Lung <- mycox_w(dat_w, "event_Lung2"),
  LymHem <- mycox_w(dat_w, "event_LymHem2"),
  NHLym <- mycox_w(dat_w, "event_NHLym2"),
  Leuk <- mycox_w(dat_w, "event_Leuk2"),
  Dig <- mycox_w(dat_w, "event_Dig2"),
  Pan <- mycox_w(dat_w, "event_Pan2"),
  CR <- mycox_w(dat_w, "event_CR2"),
  Uri <- mycox_w(dat_w, "event_Uri2"),
  Ki <- mycox_w(dat_w, "event_Ki2")
)

combined_results <- do.call(rbind, results_list)
combined_results <- combined_results %>% mutate(
  case = format(round(all_case, digits = 0), nsmall = 0), 
  ever_case_per = format(round(ever_case_per, digits = 1), nsmall = 1), 
  number = str_c(case, " (", ever_case_per, "%)"),
  
  estimate1 = format(round(mod1_est, digits = 2), nsmall = 1), 
  lower1 = format(round(mod1_low, digits = 2), nsmall = 1), 
  upper1 = format(round(mod1_high, digits = 2), nsmall = 1),
  est1 = str_c(estimate1, " (", lower1, ", ", upper1, ")"),
  
  estimate2 = format(round(mod2_est, digits = 2), nsmall = 1), 
  lower2 = format(round(mod2_low, digits = 2), nsmall = 1), 
  upper2 = format(round(mod2_high, digits = 2), nsmall = 1),
  est2 = str_c(estimate2, " (", lower2, ", ", upper2, ")"),
  
  estimate3 = format(round(mod3_est, digits = 2), nsmall = 1), 
  lower3 = format(round(mod3_low, digits = 2), nsmall = 1), 
  upper3 = format(round(mod3_high, digits = 2), nsmall = 1),
  est3 = str_c(estimate3, " (", lower3, ", ", upper3, ")")) |>
  select(cancer_type, number, est1, est2, est3) 



# Making plots for Lexie  ---- 
plot1 <- plot1 |> 
  mutate(cancer = case_when(
    cancer_type == "event_Mel2" ~ "Melanoma",
    cancer_type == "event_Thy2" ~ "Thyroid Cancer",
    cancer_type == "event_Lung2" ~ "Lung Cancer",
    cancer_type == "event_LymHem2" ~ "Lymphoma/Hematologic Cancer",
    cancer_type == "event_NHLym2" ~ "Non-Hodgkin Lymphoma",
    cancer_type == "event_Leuk2" ~ "Leukemia",
    cancer_type == "event_Dig2" ~ "Digestive Cancer",
    cancer_type == "event_Pan2" ~ "Pancreatic Cancer",
    cancer_type == "event_CR2" ~ "Colorectal Cancer",
    cancer_type == "event_Uri2" ~ "Urinary Cancer",
    cancer_type == "event_Ki2" ~ "Kidney Cancer"
  )) 

plot2 <- plot2 |> 
  mutate(cancer = case_when(
    cancer_type == "event_Mel2" ~ "Melanoma",
    cancer_type == "event_Thy2" ~ "Thyroid Cancer",
    cancer_type == "event_Lung2" ~ "Lung Cancer",
    cancer_type == "event_LymHem2" ~ "Lymphoma/Hematologic Cancer",
    cancer_type == "event_NHLym2" ~ "Non-Hodgkin Lymphoma",
    cancer_type == "event_Leuk2" ~ "Leukemia",
    cancer_type == "event_Dig2" ~ "Digestive Cancer",
    cancer_type == "event_Pan2" ~ "Pancreatic Cancer",
    cancer_type == "event_CR2" ~ "Colorectal Cancer",
    cancer_type == "event_Uri2" ~ "Urinary Cancer",
    cancer_type == "event_Ki2" ~ "Kidney Cancer"
  )) |> 
  mutate(Frequency = rep(c("Less frequent", "More frequent"), each = 1, times = 11))

library(ggplot2)
library(forcats)

plot1 |> 
  mutate(cancer = factor(cancer, 
                         levels = c("Melanoma", "Thyroid Cancer", "Lung Cancer", 
                                    "Lymphoma/Hematologic Cancer", "Non-Hodgkin Lymphoma", "Leukemia", 
                                    "Digestive Cancer", "Pancreatic Cancer", "Colorectal Cancer", 
                                    "Urinary Cancer", "Kidney Cancer")), 
         cancer = fct_rev(cancer)) |>
  ggplot(aes(x=mod2_est, y = cancer)) + 
  geom_point(aes(x=mod2_est, y = cancer), size = 2) + 
  geom_errorbar(aes(xmin = mod2_low, xmax = mod2_high), width = 0.2)  + 
  geom_vline(xintercept = 1,  linetype = "dashed") + 
  scale_x_log10() + 
  xlab("Adjusted HR (95% CI)") + 
  ylab("") + 
  theme_bw()


plot2 |> 
  mutate(cancer = factor(cancer, 
                         levels = c("Melanoma", "Thyroid Cancer", "Lung Cancer", 
                                    "Lymphoma/Hematologic Cancer", "Non-Hodgkin Lymphoma", "Leukemia", 
                                    "Digestive Cancer", "Pancreatic Cancer", "Colorectal Cancer", 
                                    "Urinary Cancer", "Kidney Cancer")), 
         cancer = fct_rev(cancer)) |>
  mutate(Frequency = factor(Frequency, levels = c( "More frequent", "Less frequent"))) |> 
  ggplot(aes(x=mod2_est, y = cancer, color=Frequency, shape = Frequency)) + 
  geom_point(aes(x=mod2_est, y = cancer), size = 2, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(xmin = mod2_low, xmax = mod2_high), width = 0.2, position = position_dodge(width = 0.5))  + 
  geom_vline(xintercept = 1,  linetype = "dashed") + 
  scale_color_manual(values = c("Less frequent" = "#303bd1", "More frequent" = "#D13030")) +  # Assign colors manually
  scale_shape_manual(values = c("Less frequent" = 17, "More frequent" = 16)) +  # Assign shape manually
  scale_x_log10() + 
  xlab("Adjusted HR (95% CI)") + 
  ylab("") + 
  theme_bw()
