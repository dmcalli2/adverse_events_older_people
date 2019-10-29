#01_ctg_SAEs_total
library(tidyverse)

aact <- readRDS("../Trial_identify/clinical_trials_august_2017/ctg/Output_data/all_aact_tables.Rds")
load("../Trial_identify/clinical_trials_august_2017/ctg/Output_data/trials_conditions_lookup_atc_code_lookup.Rdata")

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

## Examine SAEs
sae <- aact$reported_events %>%
  filter(adverse_event_term == "Total, serious adverse events")
sae <- sae %>%
  select(id, nct_id, result_group_id, ctgov_group_code, time_frame,
         default_vocab, subjects_affected, subjects_at_risk) %>%
  as_tibble()
sae <- sae %>%
  inner_join(aact$result_groups %>% select(-id)) %>%
  select(-ctgov_group_code, -result_type, -result_group_id)

## Identify trials for which dont have SAEs and add in treatment arms to
## create lookup table
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


