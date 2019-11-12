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

trial_comparison2 <- trial_comparison %>% 
  left_join(aact$design_group_interventions %>% select(-id) %>% distinct()) %>% 
  left_join(aact$design_groups %>% select(arm_type = group_type, arm_title = title, design_group_id = id) %>% distinct()) %>% 
  select(nct_id, arm_id = design_group_id, arm_type, arm_title, intervention_id, intvn_text_original = ctg_str_original,  everything()) %>% 
  rename(intvn_txt_n = ctg_n,  intvn_txt_cmpnt = ctg_str, intvn_txt_lbl = type_label, who_atc_name = unq_name) %>% 
  arrange(nct_id, arm_id, intervention_id)


intvn_arms <- aact$interventions %>% 
  filter(name %in% trial_comparison$ctg_str_original)


arms_known <- arms %>% 
  mutate(arms_compare = str_sub(arms_compare, 4)) %>% 
  separate(arms_compare, into = c("arm_a", "arm_b"), sep = "v")
arms_a <- arms_known %>% 
  select(nct_id, arm = arm_a, drug = a) 
arms_b <- arms_known %>% 
  select(nct_id, arm = arm_b, drug = b) 
arms_new <- bind_rows(a = arms_a, b = arms_b) %>% 
  arrange(nct_id, arm, drug) %>% 
  distinct() %>% 
  group_by(nct_id, drug) %>% 
  summarise(probable_arm = paste(arm, collapse = ", ")) %>% 
  ungroup()

trial_comparison <- trial_comparison %>% 
  left_join(arms_new %>% rename(unq_name = drug)) %>% 
  left_join(atc_lbl %>% rename(unq_name = title)) %>% 
  mutate(code_shrt = str_sub(code, 1, 5)) %>% 
  left_join(atc_lbl %>% rename(code_shrt = code)) %>% 
  select(-code_shrt)
write_csv(trial_comparison, "Data_extraction/01_trial_drug_comparisons.csv")


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
                      "Gender", "Gender, Customized")) %>% 
  left_join(aact$result_groups %>% select(nct_id, ctgov_group_code, result_type, group_name = title)) %>% 
  select(nct_id, ctgov_group_code, group_name, everything()) %>% 
  arrange(nct_id, title, ctgov_group_code)


## Examine AEs/SAEs
sae <- aact$reported_events %>%
  filter(adverse_event_term %in% c("Total, serious adverse events", "Total, other adverse events"))
sae <- sae %>%
  select(id, nct_id, adverse_event_term, result_group_id, ctgov_group_code, time_frame,frequency_threshold,
         default_vocab, subjects_affected, subjects_at_risk) %>%
  as_tibble()
sae <- sae %>%
  inner_join(aact$result_groups %>% select(-id) %>% rename(group_name = title)) %>%
  select(nct_id, group_name, ctgov_group_code, everything()) %>% 
  arrange(nct_id, ctgov_group_code) 

## Identify trials for which have results


no_sae <- aact$studies %>% 
  anti_join(sae) %>% 
  select(nct_id, phase, source) 
ids <- aact$id_information %>% 
  group_by(nct_id) %>% 
  summarise(other_id = paste(id_value, collapse = "|")) %>% 
  ungroup()
no_sae <- no_sae %>% 
  left_join(ids)
no_sae$nct_id %>% unique()

no_sae_intvn <- no_sae %>% 
  inner_join(intvn_lbl %>% 
               mutate(subjects_affected = "",
                      subjects_at_risk = "") %>% 
               select(nct_id, subjects_affected,
                      subjects_at_risk,
                      unq_name, type_label, name, ctg_str, ctg_str_original))
write_csv(no_sae_intvn, "Scratch_data/sae_complete.csv")
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


