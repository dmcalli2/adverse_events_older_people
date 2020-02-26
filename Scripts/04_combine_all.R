# 04_combine_all.R

library(tidyverse)
trials_over <- read_csv("Supporting/Overview_trial_extaction.csv")

## Trials sources ----
## Jo
jo <- read_tsv("Data/jo_trials.tsv")
names(jo) <- names(jo) %>% str_to_lower()
jo <- jo %>% 
  rename(nct_id = id)
## Neave and Guy ----
ng <- read_csv("Data/final_nv_guy.csv")


## Combine Jo and ng
jo <- jo %>% 
  mutate(arm_name = "Total",
         female = numberpart - male,
         total_derived = "jo",
         fu_days = followupweeks*7) %>% 
  select(nct_id, 
         arm_name,
         ae = anyae,
         sae = sae,
         total = numberpart,
         female,
         male,
         age_m = meanage,
         age_s = sdmeanage,
         bmi_m = meanbmi,
         bmi_s = sdmeanbmi,
         total_derived,
         fu_days)

ngj <- bind_rows(ng = ng, jo = jo, .id = "source") %>% 
  rename(subjects_at_risk = total)
rm(jo, ng)

## ctg ----
# add in ctg results to 97 trials where had manual extraction of results
ngj %>% distinct(nct_id) %>% filter(!is.na(nct_id)) %>% arrange(nct_id)
## 58trials in this category
anti_join(trials_over, ngj %>% select(-source))
## Time frame. Note where time frame for SAE ws different number of weeks took this rather than the primary outcome one.
## However, did not do so wherethere was a complex instruction "eg 14 days after last dose, or 30 days after last dose for SAE"
ctg_aes <- read_csv("Data_extraction_david/06_ae_sae_results.csv")
ctg_aes_tf <- read_csv("Data_extraction_david/07_ae_sae_results_time_frame.csv") %>% 
  select(nct_id, weeks)
ctg_aes <- ctg_aes %>% 
  inner_join(ctg_aes_tf)
rm(ctg_aes_tf)

ctg_aes_wide <- ctg_aes %>% 
  select(nct_id, ctgov_group_code, arm_name_ae, subjects_at_risk, adverse_event_term, subjects_affected, weeks) %>%
  group_by(nct_id, ctgov_group_code, arm_name_ae) %>% 
  mutate(subjects_at_risk = max(subjects_at_risk, na.rm = TRUE)) %>% 
  ungroup() %>% 
  spread(adverse_event_term, subjects_affected) %>% 
  mutate(sae = `Total, serious adverse events`,
         ae = sae + `Total, other adverse events`) %>% 
  select(-`Total, serious adverse events`, -`Total, other adverse events` ) %>% 
  rename(arm_name = arm_name_ae) %>% 
  mutate(fu_days = weeks*7) %>% 
  select(-weeks)
  
ctg_aes_wide_total <- ctg_aes_wide %>% 
  group_by(nct_id) %>% 
  summarise_at(vars(subjects_at_risk, sae, ae), sum) %>% 
  ungroup() %>% 
  inner_join(ctg_aes_wide %>% 
               select(nct_id, fu_days) %>% distinct()) %>% 
  mutate(arm_name = "Total",
         ctgov_group_code = "Total")

ctg_aes_wide <- bind_rows(ctg_aes_wide, ctg_aes_wide_total)
rm(ctg_aes_wide_total)

## Need same baseline data from ctg ----
ctg_base <- read_csv("Data_extraction_david/05_baseline_results.csv")
## 55 with mean age, includes 55 with Total age
ctg_age <- ctg_base %>% 
  filter(title == "Age", param_type == "Mean", dispersion_type == "Standard Deviation", units %in% c("years", "Years")) %>% 
  select(nct_id, ctgov_group_code, arm_name_bas, age_m = param_value, age_s = dispersion_value)
  
## 56 with gender in male female counts, includes 56 with total each sex
ctg_sex <- ctg_base %>% 
  filter(title == "Gender", classification %in% c("Male", "Female"), param_type == "Number") %>% 
  select(nct_id, ctgov_group_code, arm_name_bas, classification, param_value) %>% 
  mutate(classification = str_to_lower(classification)) %>% 
  spread(classification, param_value)
ctg_sex %>% filter(arm_name_bas == "Total") %>% distinct(nct_id)

## bmi, only 8 include BMI data, class already included in one of other two bmi terms
ctg_bmi <- ctg_base %>% 
  filter(title %in% c("Body Mass Index", "Body Mass Index (BMI)", 
                      "Body Mass Index Class"))
## Only 7 with mean and sd, all of which have a total
ctg_bmi <- ctg_bmi %>% 
  filter(title %in% c("Body Mass Index", "Body Mass Index (BMI)"), 
         param_type == "Mean", dispersion_type == "Standard Deviation",
         units %in% c("kg/m^2", "kilograms/square meter", "kilogram / square meter", 
                      "kilogram/meter^2", "kg/m2")) %>% 
  select(nct_id, ctgov_group_code, arm_name_bas, bmi_m = param_value, bmi_s = dispersion_value) 
ctg_bmi %>% 
  filter(arm_name_bas == "Total")

## Join all together, no increase in numbers compared for original data, ie no duplicates
ctg_bas2 <- bind_rows(ctg_base,
                      ctg_age,
                      ctg_sex,
                      ctg_bmi) %>% 
  distinct(nct_id, ctgov_group_code, arm_name_bas)
ctg_bas2 <- ctg_bas2 %>% 
  left_join(ctg_age) %>% 
  left_join(ctg_sex) %>% 
  left_join(ctg_bmi) %>% 
  arrange(nct_id, ctgov_group_code)
ctg_base %>% 
  distinct(nct_id, ctgov_group_code)

### get ones wiht age, sex and BMI in non-standard formats
no_age <- c("NCT00439738", "NCT00591578", "NCT00696241", "NCT00696436", 
            "NCT00706134")
no_sex <- c("NCT00439738", "NCT00546754", "NCT00591578", "NCT01785472")
no_bmi <- c("NCT00149227", "NCT00219141", "NCT00281580", "NCT00368277", 
            "NCT00413049", "NCT00413413", "NCT00425373", "NCT00439738", "NCT00529451", 
            "NCT00546754", "NCT00550953", "NCT00553267", "NCT00558064", "NCT00558428", 
            "NCT00591578", "NCT00614380", "NCT00624052", "NCT00631917", "NCT00649389", 
            "NCT00687973", "NCT00696241", "NCT00696436", "NCT00698646", "NCT00699192", 
            "NCT00705575", "NCT00706134", "NCT00739973", "NCT00760266", "NCT00765674", 
            "NCT00777946", "NCT00778921", "NCT00787605", "NCT00797316", "NCT00809926", 
            "NCT00841672", "NCT00853957", "NCT00865020", "NCT00902304", "NCT00923091", 
            "NCT00926289", "NCT00927394", "NCT00931710", "NCT00942994", "NCT01001572", 
            "NCT01042392", "NCT01167153", "NCT01237223", "NCT01368536", "NCT01599104", 
            "NCT01615198", "NCT01785472", "NCT01876368", "NCT01975246")

## These are all categorical so got from published papers identified from clinicaltrials.gov except
# NCT00706134 which got from systematic review https://www.nature.com/articles/hr2012185/tables/1
# For that one only presented age for selcted 300mg and placebo, so assume is same for all treatment arms
age_bmi_ctg_manual <- read_csv(
"nct_id,ctgov_group_code,arm_name_bas,age_m,age_s,bmi_m,bmi_s
NCT00439738,B1,Valsartan/HCTZ (Hydrochlorothiazide),55.4,8.5,35.2,7.3
NCT00439738,B2,HCTZ +Amlodipine,56.5,8.6,34.8,6.9
NCT00439738,B3,Total,55.95,8.55,35,7.1
NCT00591578,B1,Azilsartan Medoxomil 40 mg QD,57.8,12.1,30.8,5.7
NCT00591578,B2,Azilsartan Medoxomil 80 mg QD,56.8,10.7,30.7,5.3
NCT00591578,B3,Valsartan 320 mg QD,58.1,10.9,31.2,5.8
NCT00591578,B4,Total,57.56666667,11.23333333,30.9,5.6
NCT00696241,B1,Azilsartan Medoxomil 20 mg QD,57.1,11.02,30.4,5.67
NCT00696241,B2,Azilsartan Medoxomil 40 mg QD,57.4,9.62,30.6,5.94
NCT00696241,B3,Azilsartan Medoxomil 80 mg QD,58.1,11.56,30,5.48
NCT00696241,B4,Olmesartan 40 mg QD,58.9,11.57,29.8,5.25
NCT00696241,B5,Placebo QD,59.4,10.53,30,4.93
NCT00696241,B6,Total,58.04439216,10.89703529,30.17772549,5.51214902
NCT00696436,B1,Azilsartan Medoxomil 40 mg QD,57,12,31.7,6
NCT00696436,B2,Azilsartan Medoxomil 80 mg QD,56,11,30.7,5.9
NCT00696436,B3,Valsartan 320 mg QD,55,11,31.1,5.5
NCT00696436,B4,Olmesartan 40 mg QD,56,11,31.1,5.5
NCT00696436,B5,Placebo QD,56,11,30.5,5.4
NCT00696436,B6,Total,55.99845081,11.21688613,31.07025562,5.684817971
NCT00706134,B1,Placebo,72.3,5.25,,
NCT00706134,B2,Aliskiren 75 mg,72.1,5.47,,
NCT00706134,B3,Aliskiren 150 mg,72.1,5.47,,
NCT00706134,B4,Aliskiren 300 mg,72.1,5.47,,
NCT00706134,B5,Total,72.1,5.47,,"
)

ctg_age <- bind_rows(ctg_age,
                      age_bmi_ctg_manual %>% select(-bmi_m, -bmi_s))

## For sex, one trial is very unusally coded take this one separately, add in adverse event data at this point
NCT00591578 <- ctg_base %>% 
  filter(nct_id == "NCT00591578", title == "Gender, Customized", classification %in% c("Male (Double Blind Phase)",
                                                                                       "Female (Double Blind Phase)")) %>% 
  inner_join(ctg_aes_wide %>% select(nct_id, arm_name_bas = arm_name, subjects_at_risk, sae, ae, fu_days)) %>% 
  mutate(classification = if_else(classification == "Male (Double Blind Phase)", "male", "female")) %>% 
  select(nct_id, arm_name = arm_name_bas, subjects_at_risk, sae, ae, fu_days, classification, param_value) %>% 
  spread(classification, param_value, fill = 0L)
NCT00591578 <- NCT00591578 %>% 
  inner_join(ctg_age %>% filter(nct_id == "NCT00591578") %>% rename(arm_name = arm_name_bas))

## Remainder are very straightforward
no_sex2 <- ctg_base %>% 
  filter(nct_id %in% no_sex,
         !nct_id == "NCT00591578",
         title %in% c("Sex: Female, Male", 
                       "Gender")) %>% 
  select(nct_id, ctgov_group_code, arm_name_bas, category, param_value) %>% 
  spread(category, param_value)
names(no_sex2) <- str_to_lower(names(no_sex2))
ctg_sex <- bind_rows(ctg_sex,
                     no_sex2)
rm(no_sex2)

## check bmi, just categories, not feasible to find all
ctg_bmi <- bind_rows(ctg_bmi,
                     age_bmi_ctg_manual %>% select(-age_m, -age_s) %>% filter(!is.na(bmi_m)))

no_bmi2 <- ctg_base %>% 
  filter(nct_id %in% no_bmi, title == "Body Mass Index") 

## Combine data without missing rows
ctg_bas_final <- bind_rows(
                      ctg_age,
                      ctg_sex,
                      ctg_bmi) %>% 
  distinct(nct_id, ctgov_group_code, arm_name_bas)
ctg_bas_final <- ctg_bas_final %>% 
  left_join(ctg_age) %>% 
  left_join(ctg_sex) %>% 
  left_join(ctg_bmi) %>% 
  arrange(nct_id, ctgov_group_code) %>% 
  rename(arm_name = arm_name_bas)

## 60 to being with, NCT00591578 missing sex is in separate file NCT00591578
ctg_base %>% 
  distinct(nct_id)
# 60 at end, none missing for age, 
ctg_bas_final %>% 
  filter(arm_name == "Total") %>% 
  mutate_at(vars(age_m:bmi_s), is.na) %>% 
  summarise_at(vars(age_m:bmi_s), mean)
# age_m age_s female   male bmi_m bmi_s
# <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl>
#   1     0     0 0.0167 0.0167 0.817 0.817

## Join base and adverse events, one in bas that is not in aes, none in aes that is not in bas
setdiff(ctg_bas_final$nct_id, ctg_aes_wide$nct_id)
setdiff(ctg_aes_wide$nct_id, ctg_bas_final$nct_id)

ctg_aes_bas_mtch <- ctg_bas_final %>% 
  filter(!nct_id == "NCT00591578") %>% 
  rename(ctgov_group_code_bas = ctgov_group_code) %>% 
  inner_join(ctg_aes_wide %>% rename(ctgov_group_code_aes = ctgov_group_code))

ctg_aes_bas_no_mtch <- ctg_bas_final %>% 
  filter(!nct_id == "NCT00591578") %>% 
  rename(ctgov_group_code_bas = ctgov_group_code) %>% 
  anti_join(ctg_aes_bas_mtch) 

ctg_ae_rematch <- read_tsv("Created_metadata/arm_match_ctg.txt") %>% 
  filter(!is.na(match)) %>% 
  select(nct_id, arm_name, match)

ctg_aes_bas_no_mtch2 <- ctg_aes_bas_no_mtch %>% 
  inner_join(ctg_ae_rematch) %>% 
  select(-arm_name) %>% 
  rename(arm_name = match) %>% 
  inner_join(ctg_aes_wide %>% rename(ctgov_group_code_aes = ctgov_group_code))

ctg_aes_bas <- bind_rows(ctg_aes_bas_mtch,
                          ctg_aes_bas_no_mtch2)
rm(age_bmi_ctg_manual, ctg_ae_rematch, ctg_aes, ctg_aes_bas_mtch, ctg_aes_bas_no_mtch,
   ctg_aes_bas_no_mtch2, ctg_aes_wide, ctg_age, ctg_bas_final,
   ctg_bas2, ctg_base, ctg_bmi, ctg_sex, no_age, no_bmi, no_bmi2, no_sex, trials_over)

bas_aes <- bind_rows(ngj, 
                     ctg_aes_bas %>% mutate(source = "ctg"), 
                     NCT00591578 %>% rename(ctgov_group_code_bas = ctgov_group_code) %>% mutate(source = "ctg"))
rm(ngj, ctg_aes_bas, NCT00591578)
bas_aes %>% 
  # filter(arm_name == "Total") %>% 
  distinct(nct_id)

bas_aes %>%
  group_by(nct_id) %>% 
  mutate(total = any(arm_name == "Total")) %>% 
  ungroup() %>% 
  filter(!total)

## One trial without totals NCT00538486, calcualte this, same number in each arm
NCT00538486 <- bas_aes %>% 
  filter(nct_id == "NCT00538486") %>% 
  mutate_at(vars(subjects_at_risk, female, male, ae, sae), sum, na.rm = TRUE) %>% 
  mutate_at(vars(age_m, age_s, bmi_m, bmi_s, fu_days), mean, na.rm = TRUE) %>% 
  slice(1) %>% 
  mutate(arm_name = "Total")

bas_aes <- bind_rows(bas_aes, NCT00538486) %>% 
  arrange(nct_id, arm_name)
rm(NCT00538486)

## Final data with age, sex, bmi and adverse events for all trials
saveRDS(bas_aes, "Processed_data/age_sex_bmi_ae_sae.Rds")
write_tsv(bas_aes, "Processed_data/age_sex_bmi_ae_sae.csv")
