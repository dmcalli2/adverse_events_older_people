# 03b_gn_clean
source("Scripts/03a_gn_extract.R")
rm(list = ls())

library(tidyverse)
mydf_lst <- readRDS("Scratch_data/cleaned_guy_neave.Rds")
list2env(mydf_lst, envir = .GlobalEnv)
rm(mydf_lst)

## 
base <- base  %>% 
  mutate(measure_type = if_else(measure_type == "s", "sd", measure_type)) %>% 
  rename(arm_name_base = arm_name_definitive)

ae_sae <- ae_sae %>% 
  rename(arm_name_sae = arm_name)

## Clean baseline data further
prtcpts <- base %>% 
  filter(measure == "participants")

sex_n <- base %>% 
  filter(measure %in% c("male", "female"), measure_type == "nmbr") %>% 
  select(nct_id, arm_name_base, measure, value) %>% 
  spread(measure, value) %>% 
  left_join(prtcpts %>% select(nct_id, arm_name_base, total = value))

sex_n <- sex_n %>% 
  mutate(female = if_else(is.na(female), total - male, female),
         male = if_else(is.na(male), total - female, male),
         total = if_else(is.na(total), male+female, total))
## only two trials without complete
msng <- sex_n %>% 
  group_by(nct_id) %>% 
  summarise_at(vars(total, male, female), ~ any(is.na(.x))) %>% 
  ungroup() %>% 
  filter(total | male | female) %>% 
  distinct(nct_id)


## Only examine percent where dont have n
## Convert percentages to proportions
sex_prcnt <-  base %>% 
  filter(measure %in% c("male", "female"), measure_type == "prcnt") %>% 
  anti_join(sex_n) %>% 
  mutate(value = if_else(value <1, value, value /100))

sex_prcnt <-  sex_prcnt %>% 
  select(nct_id, arm_name_base, measure, value) %>% 
  spread(measure, value) %>% 
  mutate(total = 1,
         female = if_else(is.na(female), total - male, female),
         male = if_else(is.na(male), total - female, male)) %>% 
  select(-total)

sex_prcnt_n <- sex_prcnt %>% 
  inner_join(prtcpts %>% select(nct_id, arm_name_base, total = value)) %>% 
  mutate(male = male * total,
         female = female * total) %>% 
  select(nct_id, arm_name_base, male, female, total) 

## NCT00168857 gives overall percentage male and female, but gives numbers for arms, assume same proportion mena dn
## omwne in each arm
sex_prcnt_p <- sex_prcnt %>% 
  anti_join(prtcpts %>% select(nct_id, arm_name_base)) %>% 
  select(-arm_name_base) %>% 
  left_join(prtcpts %>% select(nct_id, arm_name_base, total = value)) %>% 
  mutate(male = male *total,
         female = female * total)

sex <- bind_rows(sex_n,
                         sex_prcnt_n,
                         sex_prcnt_p)
rm(sex_n, sex_prcnt, sex_prcnt_n, sex_prcnt_p, msng)

## No sex data for 14 trials
no_sex <- ae_sae %>% 
  anti_join(sex) %>% 
  distinct(nct_id)


## Partcipant numbers not in sex table
prtcpts <- prtcpts %>% 
  anti_join(sex) %>% 
  select(nct_id, arm_name_base, total = value)

prtcpts <- bind_rows(prtcpts,
                     sex) %>% 
  arrange(nct_id, arm_name_base)
rm(sex)

## age data
age <- base %>% 
  filter(measure == "age")


age_m_s <- age %>% 
  filter(measure_type %in% c("mean", "sd")) %>% 
  select(nct_id, arm_name_base, measure_type, value) %>% 
  spread(measure_type, value) %>% 
  arrange(nct_id, arm_name_base) %>% 
  select(nct_id, arm_name_base, mean, sd, everything())

## treat other age data as missing for now, only 3 trials
age_other <- age %>% 
  anti_join(age_m_s)

## bmi
bmi <- base %>% 
  filter(measure == "bmi")
bmi %>% distinct(measure_type)

bmi_m_s <- bmi %>% 
  filter(measure_type %in% c("mean", "sd")) %>% 
  select(nct_id, arm_name_base, measure_type, value) %>% 
  spread(measure_type, value) %>% 
  arrange(nct_id, arm_name_base) %>% 
  select(nct_id, arm_name_base, mean, sd, everything())

## treat other bmi data as missing for now, only one trial
bmi_other <- bmi %>% 
  anti_join(bmi_m_s)

## save participants (by sex), age and bmi data
all_base <- list(age = age_m_s, prtcpts = prtcpts, bmi = bmi_m_s)
rm(list = setdiff(ls(), c("all_base", "ae_sae")))

## Spread base to wide
all_base <- bind_rows(all_base) %>% 
  distinct(nct_id, arm_name_base) %>% 
  left_join(all_base$prtcpts) %>% 
  left_join(all_base$age %>% rename(age_m = mean, age_s = sd)) %>% 
  left_join(all_base$bmi %>% rename(bmi_m = mean, bmi_s = sd)) %>% 
  arrange(nct_id, arm_name_base)

## Add in adverse events data
ae_sae_mtch <- ae_sae %>% 
  left_join(all_base %>% mutate(arm_name_sae = arm_name_base))
ae_no_mtch <- ae_sae_mtch %>% 
  filter(is.na(arm_name_base)) %>% 
  distinct(nct_id) 
ae_no_mtch <- ae_sae %>% 
  semi_join(all_base %>% select(nct_id)) %>% 
  semi_join(ae_no_mtch) %>% 
  left_join(all_base %>% mutate(arm_name_sae = arm_name_base)) %>% 
  arrange(nct_id, arm_name_sae, arm_name_base) %>% 
  select(nct_id, arm_name_sae, arm_name_base)

write_csv(ae_no_mtch, "clipboard")
## reviewed do for each trial
mtch_rvd <- read_tsv("Created_metadata/arm_match.txt")
mtch_rvd %>% distinct(nct_id)
mtch_rvd %>% filter(nct_id == "NCT02269176")
## 
NCT00170989 <- ae_sae %>% 
  filter(nct_id == "NCT00170989") %>% 
  inner_join(mtch_rvd %>% filter(!comment == "exclude") %>% select(nct_id, arm_name_sae, arm_name_base_new) ) %>% 
  select(-arm_name_sae) %>% 
  rename(arm_name_sae = arm_name_base_new) %>% 
  inner_join(all_base %>% mutate(arm_name_sae = arm_name_base)) %>% 
  select(-arm_name_base)

NCT00219037 <- ae_sae %>% 
  filter(nct_id == "NCT00219037") %>% 
  mutate(arm_name_sae = if_else(arm_name_sae %in% c("Aliskiren/Hydrochlorothiazide 300/12.5mg",
                                                    "Aliskiren/Hydrochlorothiazide 300/25mg"),
                                "Combo therapy", arm_name_sae)) %>% 
  group_by(nct_id, arm_name_sae) %>% 
  summarise_all(sum) %>% 
  ungroup() %>% 
  inner_join(all_base %>% mutate(arm_name_sae = arm_name_base)) %>% 
  select(-arm_name_base)

## ask Neave to review NCT00220233
NCT00220233 <- ae_sae %>% 
  filter(nct_id == "NCT00220233") 

NCT00240448 <- ae_sae %>% 
  filter(nct_id == "NCT00240448") 
NCT00240448b <- all_base %>% 
  filter(nct_id == "NCT00240448") %>% 
  mutate(total = sum(total), arm_name_sae = "Total") %>% 
  select(-arm_name_base) %>% 
  slice(1)
NCT00240448 <- NCT00240448 %>% 
  inner_join(NCT00240448b)

NCT00400777 <- ae_sae %>% 
  filter(nct_id == "NCT00400777") %>% 
  mutate(arm_name_sae = "Total") %>% 
  group_by(nct_id, arm_name_sae) %>% 
  summarise_all(sum) %>% 
  ungroup()
NCT00400777b <- all_base %>% 
  filter(nct_id == "NCT00400777") %>% 
  rename(arm_name_sae = arm_name_base)
NCT00400777 <- NCT00400777 %>% 
  inner_join(NCT00400777b)

NCT00409760 <- ae_sae %>% 
  filter(nct_id == "NCT00409760") %>% 
  inner_join(all_base %>% rename(arm_name_sae = arm_name_base))

NCT00409760_tot <- NCT00409760 %>% 
  filter(arm_name_sae %in% c("Placebo",  "Valsartan/Amlodipine 160/2.5mg", 
                        "Valsartan/Amlodipine 160/5mg", "Valsartan/Amlodipine 320/2.5mg", 
                        "Valsartan/Amlodipine 320/5mg", "Valsartan/Amlodipine 40/2.5mg", 
                        "Valsartan/Amlodipine 40/5mg", "Valsartan/Amlodipine 80/2.5mg", 
                        "Valsartan/Amlodipine 80/5mg"))  %>% 
  mutate(arm_name_sae = "Total") %>% 
  group_by(nct_id, arm_name_sae) %>% 
  mutate(age_m = weighted.mean(age_m, total)) %>% 
  mutate_at(vars(ae, sae, total, female, male), sum) %>% 
  slice(1)

NCT00409760 <- bind_rows(NCT00409760,
                         NCT00409760_tot)

NCT01928628 <- ae_sae %>% 
  filter(nct_id == "NCT01928628") %>% 
  inner_join(all_base %>% 
               filter(nct_id == "NCT01928628") %>% 
               rename(arm_name_sae = arm_name_base))

NCT02177409 <- ae_sae %>% 
  filter(nct_id == "NCT02177409")
NCT02177409b <- all_base %>% 
  filter(nct_id == "NCT02177409") %>% 
  mutate(arm_name_sae = "Total") %>% 
  select(-arm_name_base) %>% 
  group_by(nct_id, arm_name_sae) %>% 
  summarise_all(sum)
NCT02177409 <- NCT02177409 %>% 
  inner_join(NCT02177409b)

## Note have numerator but not denominator for 100, so drop
NCT02200653 <- ae_sae %>% 
  filter(nct_id == "NCT02200653")
NCT02200653b <- all_base %>% 
  filter(nct_id == "NCT02200653")
NCT02200653 <- NCT02200653 %>% 
  inner_join(NCT02200653b %>% 
               rename(arm_name_sae = arm_name_base))

NCT02242318 <- ae_sae %>% 
  filter(nct_id == "NCT02242318") %>% 
  mutate(arm_name_sae = if_else(arm_name_sae == "Valsartan 160mg", "Valsartan", arm_name_sae)) %>% 
  inner_join(all_base %>% rename(arm_name_sae = arm_name_base)) 

NCT02269176 <- ae_sae %>% 
  filter(nct_id == "NCT02269176") %>% 
  mutate(arm_name_sae = "Total") %>% 
  group_by(nct_id, arm_name_sae) %>% 
  summarise_all(sum) %>% 
  ungroup() %>% 
  inner_join(all_base %>% rename(arm_name_sae = arm_name_base)) 

no_match_resolve <- bind_rows(NCT00170989, NCT00219037, NCT00220233, NCT00240448, 
                              NCT00400777, NCT00409760, NCT01928628, NCT02177409, 
                              NCT02200653, NCT02242318, NCT02269176)
setdiff(ae_no_mtch$nct_id, no_match_resolve$nct_id)
setdiff(no_match_resolve$nct_id, ae_no_mtch$nct_id)

ae_sae_all <- bind_rows(ae_sae_mtch %>% filter(!is.na(arm_name_base)) %>% rename(arm_name = arm_name_sae) %>% select(-arm_name_base),
                        no_match_resolve  %>% rename(arm_name = arm_name_sae))

## Next step need to calcualte total for events, (ns and %s)
saveRDS(ae_sae_all, "Scratch_data/cleaned_guy_neave_combined.Rds")

rm(list = ls())


