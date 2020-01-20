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
  semi_join(trials)

intvn_lbl <- intvn_lbl %>% 
  semi_join(trials)

aact <- map(aact, function(x) {
  if("nct_id" %in% names(x)) x %>% filter(nct_id %in% trials$nct_id)
})

## Identify Jo trials
jo_trials <- read_tsv("Data/jo_trials.tsv")
names(jo_trials) <- names(jo_trials) %>% str_to_lower()
jo_trials <- jo_trials %>% 
  rename(nct_id =id)

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


trial_comparison_wide <- trial_comparison %>% 
  select(-intvn_txt_n, -intvn_txt_cmpnt, -intvn_txt_lbl) %>% 
  distinct() %>% 
  group_by(nct_id, arm_id, intvn_text_original) %>% 
  mutate(ordr = paste0("who_names", seq_along(nct_id))) %>% 
  ungroup() %>% 
  spread(ordr, who_atc_name) %>% 
  arrange(nct_id, arm_id, intervention_id)

trial_comparison_wide <- trial_comparison_wide %>% 
  mutate(in_jo = nct_id %in% jo_trials$nct_id) %>% 
  select(nct_id, in_jo, everything())

setdiff(jo_trials$nct_id, aact$studies$nct_id)

jo_classes <- trial_comparison_wide %>% 
  filter(in_jo)

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
  mutate(in_jo = nct_id %in% jo_trials$nct_id,
         baseline_res = nct_id %in% baseline$nct_id,
         sae_res = nct_id %in% sae$nct_id) %>%
  arrange(first_received_date, nct_id) %>% 
  left_join(csr_chk %>% select(nct_id, source, csr_status)) %>% 
  mutate(csr_status = case_when(in_jo ~ "extracted_jo",
                            baseline_res & sae_res ~ "Not reqd, in clinicaltrials.gov", 
                            TRUE ~ csr_status))

 
## Add text fields for review 
text_info <- aact$studies %>% 
  select(nct_id, brief_title, official_title) %>% 
  inner_join(aact$brief_summaries %>% 
        select(-id, brief_summary = description)) %>% 
  left_join(aact$detailed_descriptions %>% select(-id, detailed_summary = description)) %>% 
  left_join(aact$eligibilities %>% select(nct_id, eligibility_criteria = criteria))
text_info <- text_info %>% 
  mutate(official_title = if_else(brief_title == official_title, "", official_title))

text_info <- text_info %>% 
  separate(eligibility_criteria, sep = fixed("Exclusion Criteria:\n\n"), into = c("inclusion", "exclusion"), remove = FALSE) %>% 
  mutate(exclusion = paste0("Exclusion Criteria:\n\n", exclusion))

SecondaryFind <- function(incld_excld = text_info$inclusion){
  a_present <- str_detect(incld_excld, "(secondary hypertension)|(hypertension secondary)")
  a <- str_locate(incld_excld, "(secondary hypertension)|(hypertension secondary)")
  a <- as_tibble(a)
  secondary <- pmap(list(a$start, a$end, incld_excld), function(x, y, z) str_sub(z, start = x-50, end = y+50))
  secondary[!a_present] <- ""  
  map_chr(secondary, paste, collapse = "")
}


text_info$secondary_inclusion <- SecondaryFind(text_info$inclusion)
text_info$secondary_exclusion <- SecondaryFind(text_info$exclusion)

## Read in parsed eligibility criteria
elig <- readxl::read_excel("../ihw_eligibility/eligibility_exclusions_final_calc.xlsx")
names(elig) <- names(elig) %>% str_to_lower() %>% str_replace_all("\\s{1,}|/|â", "_")
names(elig) <- names(elig) %>% str_replace_all("_{1,}", "_")
setdiff(trials_over$nct_id, elig$trial_id)
elig <- elig %>% 
  rename(nct_id = trial_id) %>% 
  semi_join(trials_over)

## blood pressure criteria ----
elig_bp <- elig %>% 
  filter(str_detect(term %>% str_to_lower(), "\\b(hypertens|(blood press)|(sbp)|(dbp))")) %>% 
  distinct() %>% 
  arrange(nct_id, desc(inclusion_exclusion), negated, term) %>% 
  select(nct_id, inclusion_exclusion, negated, term, criteria_sentences, everything())
elig_no_bp <- elig %>% 
  anti_join(elig_bp %>% select(nct_id))

elig_bp_lvl <- elig_bp %>% 
  filter(!is.na(measurement_constraint)) %>% 
  rename(bp_lvls = measurement_constraint)
elig_bp_no_lvl <- elig_bp %>% 
  anti_join(elig_bp_lvl %>% select(nct_id)) %>% 
  select(-measurement_constraint)

write_csv(elig_bp_no_lvl %>% distinct(inclusion_exclusion, negated, criteria_sentences), "clipboard")
elig_bp_lvl_rv <- read_csv("Created_metadata/blood_pressure_inclusion.csv")
elig_bp_no_lvl <- elig_bp_no_lvl %>% 
  inner_join(elig_bp_lvl_rv)

elig_bp_lvl <- bind_rows(elig_bp_lvl, elig_bp_no_lvl)
elig_bp_lvl %>% distinct(nct_id)
elig_bp_lvl <- trials_over %>% 
  select(nct_id) %>% 
  left_join(elig_bp_lvl) %>% 
  arrange(inclusion_exclusion)

## Find secondary hypertension ----
elig %>% filter(str_detect(criteria_sentences, "secondary")) %>% distinct(criteria_sentences, term) %>% 
  write_csv("temp.csv")
secondary_txt <- read_csv("Created_metadata/secondary.csv")
elig_secondary <- elig %>% 
  semi_join(secondary_txt) %>% 
  filter( (inclusion_exclusion == "Exclusion" & negated == FALSE ) |
          (inclusion_exclusion == "Inclusion" & negated == TRUE))

essential <- elig %>% 
  anti_join(elig_secondary %>% select(nct_id)) %>% 
  filter(str_detect(criteria_sentences %>% str_to_lower(), "essential|idiopathic")) %>% 
  distinct(criteria_sentences) %>% 
  write_csv("clipboard")
elig_essential <-  elig %>% 
  semi_join(essential) %>% 
  filter( (inclusion_exclusion == "Exclusion" & negated == TRUE ) |
            (inclusion_exclusion == "Inclusion" & negated == FALSE))
elig_not_secondary <- bind_rows(essential_inclusion = elig_essential,
                                secondary_exclusion = elig_secondary,
                                .id = "essential_secondary")

elig_not_secondary %>% 
  distinct(nct_id)
elig_not_secondary <- trials_over %>% 
  select(nct_id) %>% 
  left_join(elig_not_secondary)

## Find drug usage ----
drug_elig <- elig %>% 
  filter(classification_result == "Drug Exposure") %>% 
  distinct(term) %>% as.data.frame() %>% 
  write_csv("clipboard")
drug_elig <- read_csv("Created_metadata/drugs_eligibility.csv")
drug_elig <- elig %>% 
  inner_join(drug_elig)
drug_elig %>% distinct(criteria_sentences, drug_class, inclusion_exclusion, negated) %>% write_csv("temp.csv")

## write files
write_csv(trials_over, "Data_extraction/01_trials_overview.csv")
write_csv(primary, "Data_extraction/02_primary_outcomes.csv")
write_csv(elig, "Data_extraction/03_age_sex_criteria.csv")
write_csv(trial_comparison_wide, "Data_extraction/04_drug_comparisons.csv")
write_csv(baseline, "Data_extraction/05_baseline_results.csv")
write_csv(sae, "Data_extraction/06_ae_sae_results.csv")
write_csv(sae_time_frame, "Data_extraction/07_ae_sae_results_time_frame.csv")
write_csv(text_info, "Data_extraction/08_text_summaries.csv")


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

## Create table with arms in it and ask for clarification and assignment
trial_comparison_slct <- trial_comparison_wide %>% 
  inner_join(trials_over %>% select(nct_id, csr_status)) %>% 
  filter(!csr_status %in% c("extracted_jo",
                               "Not reqd, in clinicaltrials.gov"))
write_csv(trial_comparison_slct, "clipboard")


trial_comparison_wide %>% 
  inner_join(trials_over %>% select(nct_id, csr_status)) %>% 
  filter(csr_status %in% c("extracted_jo",
                            "Not reqd, in clinicaltrials.gov")) %>% 
  distinct(nct_id)
