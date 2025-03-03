# Hair Straightener and Cancers (clean) -- 

library(tidyverse)
library(haven)

# Load data ----
data1 <- read_sas('raw_data/dr00340_01_03.sas7bdat')
formats <- read_sas("raw_data/sisformats_data.sas7bdat")

set <- data1 %>% mutate(age_diff = abs(HR_MenopauseAgeExact_T0-AgeExact_T0)) %>%
   filter(age_diff <= 5) %>% 
   mutate(income = case_when(SE19Impute %in% c(1,2) ~ 0,  # < $50,000
                             SE19Impute %in% c(3) ~ 1, 
                             SE19Impute %in% c(4,5) ~ 2),
          educ = case_when(SE18 %in% c(1,2,3,4,5) ~ 0,  # HS or less 
                           SE18 %in% c(6,7) ~ 1,    # some college
                           SE18 %in% c(8,9,10) ~ 2))  # college or higher 
  
  
set %>% summarise(age_median = median(AgeExact_T0, na.rm = TRUE), 
            age_min = min(AgeExact_T0, na.rm = TRUE),
            age_max = max(AgeExact_T0, na.rm = TRUE)) 

round(prop.table(table(set$SE_RACE4))*100, 2)
round(prop.table(table(set$income))*100, 2)
round(prop.table(table(set$educ))*100, 2)

# Remove pre-baseline cancer cases ----
# 1. Based on the FU_Ca_Event (use Ca instead of Ca_Inv because we want to remove all the cancer outcomes)
data2 <- data1 %>% 
  filter(!is.na(FU_Ca_Event)) # 50884-47770 = 3114
# b    u    w 
# 2847  258    9 

# 2. Based on all the individual cancer event variables
data3 <- data2 %>%
  # based on the cancer variables because we want to exclude pre-baseline in situ as well 
  filter(if_all(c(FU_BCInvD_Event, 
                  FU_ReproCa_Event, FU_OvCa_Event, FU_UECa_Event, #Cervical, vulvar, vaginal
                  FU_ThyCa_Event, 
                  FU_DigesCa_Event, FU_CRCa_Event, FU_PanCa_Event, #Small intestinal, anal, biliary tract, liver, stomach and esophageal
                  FU_LungCa_Event, #Bronchus
                  FU_HdNkCa_Event, 
                  DR340_UrinCa_Event, DR340_BldrCa_Event, DR340_KiCa_Event, #Renal pelvis, ureter, urethra
                  DR340_MelCa_Event,
                  DR340_LymHemCa_Event, DR340_NHLym_Event, DR340_HLym_Event, DR340_MMyel_Event, DR340_Leuk_Event, #Myelodysplastic syndrome, Polycythemia vera
                  FU_CrCNSTum_Event #and other cancers
  ), ~ !is.na(.))) #47770-47660 = 110

data4 <- data3 %>%
  # based on the invasive cancer variables again because there might be inconsistency 
  filter(if_all(c(FU_BCInvNoD_Event, 
                  FU_ReproCaInv_Event, 
                  
                  FU_DigesCaInv_Event, 
                  
                  FU_HdNkCaInv_Event, 
                  DR340_UrinCaInv_Event, DR340_BldrCaInv_Event, DR340_KiCaInv_Event, 
                  DR340_MelCaInv_Event,
                  DR340_LymHemCaInv_Event, 
                  FU_CrCNSTumInv_Event #and other cancers
  ), ~ !is.na(.))) #47660-47657=3; 3 additional people being excluded using the na information in the invasive event variables

setdiff(data3$PSID, data4$PSID)
na_tag(data3[data3$PSID == '00340_212143',]$DR340_MelCaInv_Event) # one additional uncertain 
na_tag(data3[data3$PSID == '00340_220021',]$DR340_MelCaInv_Event) # one additional uncertain 
na_tag(data3[data3$PSID == '00340_229915',]$DR340_BldrCaInv_Event) # one additional uncertain 


# Change the in situ cancer to event = 0 ----
data5 <- data4 %>% 
  mutate(event_BC = FU_BCInvD_Event, # Include the breast cancer and DCIS cases 
         
         event_Repro = ifelse(FU_ReproCaInv_Event == 1 & FU_ReproCaInv_DxStage == 0 , 0, FU_ReproCaInv_Event), #0
         event_Ov = ifelse(FU_OvCa_Event == 1 & FU_OvCa_DxStage == 0, 0, FU_OvCa_Event), #0
         event_UE = ifelse(FU_UECa_Event == 1 & FU_UECa_DxStage == 0, 0, FU_UECa_Event), #0
         
         event_Thy = ifelse(FU_ThyCa_Event == 1 & FU_ThyCa_DxStage == 0, 0, FU_ThyCa_Event), #0
         
         event_Dig = ifelse(FU_DigesCaInv_Event == 1 & FU_DigesCaInv_DxType == 0 , 0, FU_DigesCaInv_Event), #0
         event_CR = ifelse(FU_CRCa_Event == 1 & FU_CRCa_DxStage == 0, 0, FU_CRCa_Event), #10
         event_Pan = ifelse(FU_PanCa_Event == 1 & FU_PanCa_DxStage == 0, 0, FU_PanCa_Event), #2
         
         event_Lung = ifelse(FU_LungCa_Event == 1 & FU_LungCa_DxStage == 0, 0, FU_LungCa_Event), #4
         event_HdNk = ifelse(FU_HdNkCaInv_Event == 1 & FU_HdNkCaInv_DxStage == 0, 0, FU_HdNkCaInv_Event), #0
         
         event_Uri = DR340_UrinCaInv_Event, 
         event_Bldr = DR340_BldrCaInv_Event, 
         event_Ki = DR340_KiCaInv_Event, 
         
         event_Mel = DR340_MelCaInv_Event, 
         
         event_LymHem = DR340_LymHemCaInv_Event, 
         event_NHLym = DR340_NHLym_Event,
         event_HLym = DR340_HLym_Event,
         event_MMyel = DR340_MMyel_Event,
         event_Leuk = DR340_Leuk_Event,
         
         event_CrCNS = FU_CrCNSTumInv_Event, 
         
         event_InvCa = FU_Ca_Inv_Event) 

# Update the EOF age and and based on EOF update event variables -----
data6 <- data5 %>% 
  mutate(eof = pmin(FU_BCInvD_EOFAgeExact,
                    FU_ReproCa_EOFAgeExact, FU_OvCa_EOFAgeExact, FU_UECa_EOFAgeExact, 
                    FU_ThyCa_EOFAgeExact, 
                    FU_DigesCa_EOFAgeExact, FU_CRCa_EOFAgeExact, FU_PanCa_EOFAgeExact, 
                    FU_LungCa_EOFAgeExact, 
                    FU_HdNkCa_EOFAgeExact, 
                    DR340_UrinCa_EOFAgeExact, DR340_BldrCa_EOFAgeExact, DR340_KiCa_EOFAgeExact, 
                    DR340_MelCa_EOFAgeExact,
                    DR340_LymHemCa_EOFAgeExact, DR340_NHLym_EOFAgeExact, DR340_HLym_EOFAgeExact, DR340_MMyel_EOFAgeExact, DR340_Leuk_EOFAgeExact, 
                    FU_CrCNSTum_EOFAgeExact,
                    FU_Ca_EOFAgeExact),
         FUtime = eof - AgeExact_T0)  

data7 <- data6 %>% filter(FUtime>0)  #4 7657-47388 = 269

data8 <- data7 %>%
  mutate(event_BC2 =  ifelse(FU_BCInvD_EOFAgeExact > eof & event_BC == 1, 0, event_BC), 
         
         event_Repro2 = ifelse(FU_ReproCa_EOFAgeExact > eof & event_Repro == 1, 0, event_Repro), 
         event_Ov2 = ifelse(FU_OvCa_EOFAgeExact > eof & event_Ov == 1, 0, event_Ov), 
         event_UE2 = ifelse(FU_UECa_EOFAgeExact > eof & event_UE == 1, 0, event_UE), 
         
         event_Thy2 = ifelse(FU_ThyCa_EOFAgeExact > eof & event_Thy == 1, 0, event_Thy), 
         
         event_Dig2 = ifelse(FU_DigesCa_EOFAgeExact > eof & event_Dig == 1, 0, event_Dig), 
         event_CR2 = ifelse(FU_CRCa_EOFAgeExact > eof & event_CR == 1, 0, event_CR), 
         event_Pan2 = ifelse(FU_PanCa_EOFAgeExact > eof & event_Pan == 1, 0, event_Pan), 
         
         event_Lung2 = ifelse(FU_LungCa_EOFAgeExact > eof & event_Lung == 1, 0, event_Lung),  
         event_HdNk2 = ifelse(FU_HdNkCa_EOFAgeExact > eof & event_HdNk == 1, 0, event_HdNk), 
         
         event_Uri2 = ifelse(DR340_UrinCa_EOFAgeExact > eof & event_Uri == 1, 0, event_Uri),   
         event_Bldr2 = ifelse(DR340_BldrCa_EOFAgeExact > eof & event_Bldr == 1, 0, event_Bldr), 
         event_Ki2 = ifelse(DR340_KiCa_EOFAgeExact > eof & event_Ki == 1, 0, event_Ki), 
         
         event_Mel2 = ifelse(DR340_MelCa_EOFAgeExact > eof & event_Mel == 1, 0, event_Mel),  
         
         event_LymHem2 = ifelse(DR340_LymHemCa_EOFAgeExact > eof & event_LymHem == 1, 0, event_LymHem), 
         event_NHLym2 = ifelse(DR340_NHLym_EOFAgeExact > eof & event_NHLym == 1, 0, event_NHLym), 
         event_HLym2 = ifelse(DR340_HLym_EOFAgeExact > eof & event_HLym == 1, 0, event_HLym), 
         event_MMyel2 = ifelse(DR340_MMyel_EOFAgeExact > eof & event_MMyel == 1, 0, event_MMyel), 
         event_Leuk2 = ifelse(DR340_Leuk_EOFAgeExact > eof & event_Leuk == 1, 0, event_Leuk), 
         
         event_CrCNS2 = ifelse(FU_CrCNSTum_EOFAgeExact > eof & event_CrCNS == 1, 0, event_CrCNS))


# Align combined cancer event with individual cancer events ----

data9 <- data8 %>% 
  mutate(event_Repro2 = ifelse(event_Ov2 == 1 | event_UE2 == 1, 1, event_Repro2), 
         event_Dig2 = ifelse(event_CR2 == 1 | event_Pan2 == 1, 1, event_Dig2),
         event_LymHem2 = ifelse(event_NHLym2 == 1 | event_HLym2 == 1 | event_MMyel2 == 1 | event_Leuk2 == 1, 1, event_LymHem2))

# Remove participants who did not answer any hair straightener use questions ----

data10 <- data9 %>% 
  # clean the vanguard data 
  mutate(ex_vanguard_12mo = case_when(P142 == 0 ~ 1, # never use
                                      P142 == 1 & is.na(P143) ~ 2, # 1-2 times; those who say ever use in P142 but missing P143, assuming they did not use the product very often
                                      P142 == 1 & !is.na(P143) ~ P143, # those who say ever use in P142 but not missing P143, use information in P143 
                                      is.na(P142) ~ P143, # those who missing P142 use information in P143
                                      TRUE ~ P143)) %>% # other conditions, use P143
  mutate(ex_main_12mo = PC103) %>%
  # combine the vanguard and main study data
  mutate(ex_12mo = case_when(is.na(ex_main_12mo) ~ ex_vanguard_12mo, 
                             TRUE ~ ex_main_12mo)) %>% # if missing main study data, use vanguard for the 12 monthes variables (not teen and 10-13)
  # create new frequency variable
  mutate(ex_12mo_freq = case_when(ex_12mo == 1 ~ 0, # never use
                                  ex_12mo %in% c(2,3) ~ 1, #<=4 times
                                  ex_12mo %in% c(4,5,6) ~ 2)) %>% #>4 times
  # create ever and never variable
  mutate(ex_12mo_ever = case_when(ex_12mo_freq == 0 ~ 0, # never use
                                  ex_12mo_freq %in% c(1,2) ~ 1)) # ever use

data11 <- data10 %>% 
  filter(!is.na(ex_12mo_ever)) #47388-46311 = 1077

# Remove participants who have any missing covariates ----
data12 <- data11 %>% 
  mutate(age = AgeExact_T0,
         race = SE_RACE4, 
         race1 = case_when(race == 1 ~ 1, 
                           race == 2 ~ 2, 
                           race %in% c(3,4) ~ 3),
         educ = case_when(SE18 %in% c(1,2,3,4,5) ~ 0,  # HS or less 
                          SE18 %in% c(6,7) ~ 1,    # some college
                          SE18 %in% c(8,9,10) ~ 2),  # college or higher 
         smoking = SM_Smoke_T0, 
         smoking1 = case_when(smoking == 0 ~ 0, 
                              smoking %in% c(1,2) ~ 1),
         BMI = BMI_T0,
         BMI_ex = EX_BMI_final, 
         BMI_cat = BMI_CDC_T0) %>%
  filter(if_all(c(age, race, educ, smoking), ~!is.na(.))) # 46311-46287 = 24

# Table 1 -----
library(table1)
table1(~ age + factor(race) + factor (educ) + factor(smoking1) | factor(ex_12mo_ever), data = data12)
data12 %>% group_by(ex_12mo_ever) %>% summarise(age= IQR(age))
data12 %>% summarise(age= IQR(age))

# Save the data ----
write_csv(data12, "clean_data/clean_data.csv")



