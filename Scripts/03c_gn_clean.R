# 03c_gn_clean
library(tidyverse)

## review "totals"
source("Scripts/03b_gn_clean.R")
rm(list = ls())
source("Scripts/00_functions.R")

ae_all <- readRDS("Scratch_data/cleaned_guy_neave_combined.Rds")
ae <- ae_all$ae_sae_all
fup <- ae_all$follow_up
rm(ae_all)


## remove duplicates introduced in merging
ae <- ae %>% 
  distinct()
ae <- ae %>% 
  filter(!(nct_id == "NCT00220233" & is.na(total)))

## Extract totals and review ones without total so know which to aggregate
ae %>% 
  count(arm_name, sort = T)
total <- ae %>% 
  filter(arm_name == "Total")
not_total <- ae %>% 
  filter(!arm_name == "Total")

not_total %>%
distinct(nct_id, arm_name, total) %>%
  filter(nct_id == "NCT00260923") %>% 
write_tsv("clipboard")

## reviewed CSRs adn/or ctg to ensure that, where ambiguous from totals/names
# all groups included in totals were mutually exclusive
## Need to reexamine NCT00220233, confusing
tot_rvd <- read_tsv("Created_metadata/allocate_groups_totals.txt")
## Note that orderning is important here, need to cal;cualte SDs first, then means, the totals

## 55 trials with totals
total_calc <- not_total %>% 
  semi_join(tot_rvd %>% filter(in_grp ==1)) %>% 
  group_by(nct_id) %>% 
  summarise(arm_name = "Total",
            age_m = weighted.mean(age_m, total),
            age_s = PoolSD(age_s, total),
            bmi_m = weighted.mean(bmi_m, total),
            bmi_s = PoolSD(bmi_s, total),
            ae = sum(ae),
            sae = sum(sae),
            total = sum(total),
            female = sum(female),
            male = sum(male)) %>% 
  ungroup() 

total_calc %>% distinct(nct_id)

## combine calculated and reported, then write to CSV as final result
total_both <- bind_rows(calculated = total_calc,
                        reported = total, .id = "total_derived") 

## combine all into single dataset for review
rv <- ae %>% 
  filter(!arm_name == "Total") %>% 
  bind_rows(total_both) %>% 
  mutate(arrange_by = if_else(arm_name == "Total", "ZZZ_total", arm_name)) %>% 
  arrange(nct_id, arrange_by)%>% 
  select(-arrange_by)

rv2 <- rv %>% 
  left_join(fup %>% filter(!is.na(fu_days)) %>% distinct())
rv2 %>% 
  group_by(nct_id, arm_name) %>% 
  mutate(n_measures = length(arm_name)) %>% 
  ungroup() %>% 
  filter(n_measures >=2)
write_csv(rv2, "Data/final_nv_guy.csv")
write_csv(rv2, "../adverse_events_older_people_ng/sae_extraction_20200131_rearranged.csv", na = "")
