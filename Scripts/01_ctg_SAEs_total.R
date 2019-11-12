#01_ctg_SAEs_total
library(tidyverse)

aact <- readRDS("../Trial_identify/clinical_trials_august_2017/ctg/Output_data/all_aact_tables.Rds")
load("../Trial_identify/clinical_trials_august_2017/ctg/Output_data/trials_conditions_lookup_atc_code_lookup.Rdata")

## Make nrowreturn function
nrowR <- function(x) {
  print(nrow(x))
  x
}

#  Read in all trials with drugs that are in older people trials
drugsearch <- c("telmisartan", "valsartan", "olmesartan", "aliskiren", "irbesartan")

trials_htn <- intvn_lbl %>% 
  filter(unq_name %in% drugsearch) %>% 
  distinct(nct_id)

trials <- trials %>% 
  semi_join(trials_htn %>%  select(nct_id)) %>% 
  filter(str_detect(conditions, "Hypertension"))
rm(trials_htn)
arms <- arms %>% 
  semi_join(trials)

conditions_lkp <- conditions_lkp %>% 
  semi_join(trials_htn)

intvn_lbl <- intvn_lbl %>% 
  semi_join(trials)

aact <- map(aact, function(x) {
  if("nct_id" %in% names(x)) x %>% filter(nct_id %in% trials$nct_id)
})

## Trial treatment comparisons ----
trial_comparison <- intvn_lbl %>% 
  select( nct_id, ctg_str_original, ctg_n, ctg_str, type_label, unq_name) %>% 
  distinct()

## Note complete join of intervention ID to ctg_strings
trial_comparison <- trial_comparison %>% 
  inner_join(aact$interventions %>% 
             rename(ctg_str_original = name) %>% 
              mutate(ctg_str_original = str_to_lower(ctg_str_original)) %>% 
               select(nct_id, intervention_id = id, ctg_str_original) %>% 
               distinct())

trial_comparison <- trial_comparison %>% 
  left_join(aact$design_group_interventions %>% select(-id) %>% distinct()) %>% 
  left_join(aact$design_groups %>% select(arm_type = group_type, arm_name = title, design_group_id = id) %>% distinct()) %>% 
  select(nct_id, arm_id = design_group_id, arm_type, arm_name, intervention_id, intvn_text_original = ctg_str_original,  everything()) %>% 
  rename(intvn_txt_n = ctg_n,  intvn_txt_cmpnt = ctg_str, intvn_txt_lbl = type_label, who_atc_name = unq_name) %>% 
  arrange(nct_id, arm_id, intervention_id)

## Primary outcome ----
primary <- aact$design_outcomes %>% 
  filter(outcome_type == "primary") %>% 
  select(nct_id, measure, description) %>% 
  distinct()


## eligibilities ----
elig <- aact$eligibilities %>% 
  select(nct_id, minimum_age, maximum_age, gender) %>% 
  distinct()

## participants for all ----
participants <- trials %>% 
  select(nct_id, enrollment)

## Baseline results ----
baseline <- aact$baseline_measurements %>% 
  filter(title %in% c("Age", "Age, Customized", "Body Mass Index", "Body Mass Index (BMI)",  "Body Mass Index Class",
                      "Gender", "Gender, Customized", "Sex: Female, Male")) %>% 
  left_join(aact$result_groups %>% select(nct_id, ctgov_group_code, result_type, group_name = title)) %>% 
  select(nct_id, ctgov_group_code, group_name, everything()) %>% 
  arrange(nct_id, title, ctgov_group_code)


baseline_xmn <- baseline %>% 
  select(nct_id, group_name) %>% 
  distinct()

baseline <- baseline %>% 
  left_join(trial_comparison %>% select(nct_id, group_name = arm_name) %>% mutate(same_name_as_drug_compare = TRUE) %>% distinct()) %>% 
  mutate(same_name_as_drug_compare = if_else(is.na(same_name_as_drug_compare), FALSE, same_name_as_drug_compare),
         name_drug_compare = if_else(same_name_as_drug_compare, group_name, "" )) %>% 
  select(nct_id, ctgov_group_code, arm_name_bas = group_name, same_name_as_drug_compare, name_drug_compare,
         classification, category, title, description, units, 
         param_type, param_value, 
         dispersion_type, dispersion_value)

## Examine AEs/SAEs
sae <- aact$reported_events %>%
  filter(adverse_event_term %in% c("Total, serious adverse events", "Total, other adverse events"))
sae <- sae %>%
  select(id, nct_id, adverse_event_term, result_group_id, ctgov_group_code, time_frame,frequency_threshold,
         default_vocab, subjects_affected, subjects_at_risk) %>%
  as_tibble()
sae <- sae %>%
  inner_join(aact$result_groups %>% select(-id) %>% rename(group_name = title)) %>%
  select(nct_id, ctgov_group_code, arm_name_ae = group_name, everything()) %>% 
  arrange(nct_id, ctgov_group_code)  %>% 
  select(-id, -result_group_id, -result_type)

sae <- sae %>% 
  left_join(baseline %>% select(nct_id, arm_name_ae = arm_name_bas) %>% mutate(same_name_as_baseline = TRUE) %>% distinct()) %>% 
  mutate(same_name_as_baseline = if_else(is.na(same_name_as_baseline), FALSE, same_name_as_baseline),
         name_baseline = if_else(same_name_as_baseline, arm_name_ae, "" )) %>% 
  select(nct_id, ctgov_group_code, arm_name_ae, same_name_as_baseline, name_baseline, everything())

## Where time frame is missing pull this from outcomes data
outcomes <- aact$outcomes %>% 
  filter(outcome_type == "Primary") %>% 
  distinct(nct_id, time_frame) %>% 
  as_tibble()
## two trials both are from baseline to 8 weeks
outcomes_diff <- outcomes %>% 
  group_by(nct_id) %>% 
  mutate(n = length(time_frame)) %>% 
  ungroup() %>% 
  filter(n >=2)

outcomes <- outcomes %>% 
  filter(!nct_id %in% c("NCT00558428", "NCT00281580"))
outcomes2 <- tibble(nct_id = c("NCT00558428", "NCT00281580"), time_frame = "Baseline and Week 8")

outcomes <- bind_rows(outcomes, outcomes2)


sae_time_frame <- sae %>% 
  inner_join(outcomes %>% rename(time_frame_primary = time_frame)) %>% 
  distinct(nct_id, time_frame, time_frame_primary)

sae$time_frame <- NULL

###
csr_chk <- read_csv("Data/csr_status_current.csv") 


## Identify trial overview
trials_over <- trials %>% 
  select(nct_id, official_title, first_received_date) %>% 
  mutate(baseline_res = nct_id %in% baseline$nct_id,
         sae_res = nct_id %in% sae$nct_id) %>%
  arrange(first_received_date, nct_id) %>% 
  left_join(csr_chk %>% select(nct_id, source, csr_status)) %>% 
  mutate(csr_status = if_else(baseline_res & sae_res, "Not reqd, in clinicaltrials.gov", csr_status))

 
## write files
write_csv(trials_over, "Data_extraction/01_trials_overview.csv")
write_csv(primary, "Data_extraction/02_primary_outcomes.csv")
write_csv(elig, "Data_extraction/03_age_sex_criteria.csv")
write_csv(trial_comparison, "Data_extraction/04_drug_comparisons.csv")
write_csv(baseline, "Data_extraction/05_baseline_results.csv")
write_csv(sae, "Data_extraction/06_ae_sae_results.csv")
write_csv(sae_time_frame, "Data_extraction/07_ae_sae_results_time_frame.csv")



## Connect to CTG
## Checked on 29th October, no new AE reports
# library(RPostgreSQL)
# drv <- dbDriver('PostgreSQL')
# username <- read_lines("myuser.ignore")
# mypassword <- read_lines("mypass.ignore")
# con <- dbConnect(drv, dbname="aact",host="aact-db.ctti-clinicaltrials.org", 
#                  port=5432, user=username, password=mypassword)
# 
# dbListTables(con) %>%  sort()
# sae_new <- dbGetQuery(con, paste0("SELECT * FROM reported_events
#                          WHERE nct_id IN('",
#                                   paste(no_sae$nct_id, collapse = "', '"),
#                                      "')"))
# dbDisconnect(con)

