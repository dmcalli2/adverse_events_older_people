library(tidyverse)
library(readxl)
library(here)

## Include following to allow compiling of report
setwd(here())


ls_shts <- excel_sheets("Data/sae_extraction_20200120.xlsx")

ls_shts <- ls_shts[-1]
names(ls_shts) <- ls_shts
shts <- map(ls_shts, ~ read_excel(path = "Data/sae_extraction_20200120.xlsx", sheet = .x)) 

## rename "comment" columns so differs across sheets
shts <- map2(shts, names(shts), ~ .x %>% set_names(str_replace_all(names(.x), "comment", 
                                                                   paste(.y, "comment", sep = "_")
                                                                   ))
             )

## cols in each and names
map(shts, names)
map(shts, ncol)

## Rows in each
map(shts, nrow)

## Trials in each
map(shts, ~ sum(!duplicated(.x$nct_id)))

## Remove null reports
null_results <- shts$results_extract %>% 
  filter(measure %in% c("Not reported", "NO SAE data in CSR summary") |
           value %in% "Not reported in clinicaltrials.gov entry")

shts$results_extract %>% 
  filter(nct_id == "NCT00316095")

shts$results_extract <- shts$results_extract %>% 
  filter(!measure %in% c("Not reported", "NO SAE data in CSR summary"),
         !value %in% "Not reported in clinicaltrials.gov entry")

## Which were not able to get baseline for
no_base <- anti_join(shts$csr_info, shts$baseline_extract)

## Which were not able to get results for
no_res <-  anti_join(shts$csr_info, shts$results_extract)

setdiff(no_base$nct_id, no_res$nct_id)
setdiff(no_res$nct_id, no_base$nct_id)
intersect(no_base$nct_id, no_res$nct_id)

tbl_incomplete <- bind_rows(no_base %>% select(nct_id), no_res %>% select(nct_id)) %>% 
  distinct() %>% 
  mutate(base = !nct_id %in% no_base$nct_id,
         res = !nct_id %in% no_res$nct_id) %>% 
  arrange(desc(base), desc(res))
tbl_incomplete

## CHeck arms match
names_match <- bind_rows(shts$baseline_extract %>% select(nct_id, arm_name = arm_name_definitive),
                         shts$results_extract %>% select(nct_id, arm_name)) %>% 
  distinct() %>% 
  semi_join(shts$baseline_extract %>% select(nct_id, arm_name = arm_name_definitive)) %>% 
  semi_join(shts$results_extract %>% select(nct_id, arm_name)) %>% 
  filter(!is.na(nct_id), !nct_id == "OR")

setdiff(shts$results_extract$arm_name, shts$baseline_extract$arm_name_definitive)
res_name_diff <- shts$results_extract %>% 
  semi_join(shts$baseline_extract %>% select(nct_id)) %>% 
  anti_join(shts$baseline_extract %>% select(nct_id, arm_name = arm_name_definitive)) %>% 
  distinct(nct_id, arm_name) %>% 
  rename(arm_name_res = arm_name)

bas_name_diff <- shts$baseline_extract %>% 
  select(nct_id, arm_name = arm_name_definitive) %>% 
  semi_join(shts$results_extract %>% select(nct_id)) %>% 
  anti_join(shts$results_extract %>% select(nct_id, arm_name)) %>% 
  distinct(nct_id, arm_name) %>% 
  rename(arm_name_base = arm_name)

names_mismatch <- bind_rows(res_name_diff %>% select(nct_id),
                            bas_name_diff %>% select(nct_id)) %>% 
  distinct() %>% 
  left_join(res_name_diff) %>% 
  left_join(bas_name_diff) %>% 
  distinct() %>% 
  arrange(nct_id, arm_name_res, arm_name_res) %>% 
  left_join(shts$csr_info %>% select(nct_id, `Extractor (eg Naeve or Guy)`))

## Examine results terms
shts$results_extract %>% 
  count(measure, measure_type) %>% 
  spread(measure_type, n, fill = 0L)

## One with an NA as reported only as "very low", I think we treat this as missing
shts$results_extract %>% 
  filter(measure == "SAE", is.na(measure_type)) %>% 
  t()


## All Rseults comments are fine
## Show consistent approach to FU time across Guy and Naeve
res_cmnt <- shts$results_extract %>% 
  select(starts_with("results_ext")) %>% 
  distinct() 
res_cmnt <- res_cmnt %>% 
  filter(!(is.na(results_extract_comment1) & is.na(results_extract_comment2) & is.na(results_extract_comment1)))
res_cmnt_rv <- res_cmnt %>% 
  inner_join(shts$results_extract %>% select(nct_id, starts_with("results_extract"))) %>% 
  distinct()
res_cmnt_rv <- res_cmnt_rv %>% 
  inner_join(shts$results_extract)

## Baseline comments all concern fact that ITT population was slightly higher. Very small differences
shts$baseline_extract %>%
  select(nct_id, starts_with("baseline"), value) %>%
  filter(!(is.na(baseline_extract_comment1) & is.na(baseline_extract_comment2) & is.na(baseline_extract_comment3))) %>% 
  distinct()

## Baseline terms, all look fine, no further cleaning
shts$baseline_extract %>% 
  group_by(measure, measure_type) %>% 
  summarise(sum(!duplicated(nct_id))) %>% 
  ungroup()

