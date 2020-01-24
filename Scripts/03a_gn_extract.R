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

## Gather SAE and AE data
rslts <- shts$results_extract %>% 
  mutate(measure = case_when(
    measure %in% c("AE", "TEAE", "TEAEs") ~ "AE",
    measure %in% c("SAE", "Serious TEAE", "Serious TEAEs") ~ "SAE",
    measure %in% 
      c("Death", "Deaths",
      "No patient died during the single-blind or double-blind periods of the study. In the washout period and follow-up period, two deaths were reported.") ~ "Death",
    measure %in% c("Drug-related SAE") ~ "Drug-related SAE" ))

rslts %>% 
  distinct(measure_type) %>% 
  arrange(measure_type)

## Note the NA for measure type is where it only says "very low"
rslts %>% filter(is.na(measure_type))

## Separate results into measures
rslts_n <-  rslts %>% 
  filter(measure_type == "n")

rslts_nprcnt <- rslts %>% 
  filter(measure_type == "n (%)") %>% 
  separate(value, into = c("value_n", "value_pcnt"), sep = "\\(") %>% 
  mutate(value = parse_integer(value_n),
         measure_type = "n") %>% 
  select(-value_n, -value_pcnt)

rslts_prcnt <- rslts %>% 
  filter(measure_type == "%") 
rslts_prcnt %>% semi_join(rslts_n %>% select(nct_id, arm_name))
## All 4 trials where have baseline n
bline_main <- shts$baseline_extract %>% 
  semi_join(rslts_prcnt %>% select(nct_id, arm_name_definitive = arm_name)) %>% 
  filter(measure == "participants")  %>% 
  select(nct_id, measure, measure_type, arm_name = arm_name_definitive, value)

## One trial without arm level participant numbers
## For now (will aggregate later) will assume 1:1 randomisation
bline_mismatch <- shts$baseline_extract %>% 
  semi_join(rslts_prcnt %>% select(nct_id)) %>% 
  anti_join(rslts_prcnt %>% select(nct_id, arm_name_definitive = arm_name)) %>% 
  select(nct_id, measure, measure_type, arm_name = arm_name_definitive, value) 
bline_mismatch <- bline_mismatch %>% 
  filter(nct_id %in% "NCT02269176") %>% 
  select(-arm_name) %>% 
  left_join(rslts_prcnt %>% select(nct_id, arm_name)) %>% 
  mutate(value = as.integer(value)/2)
bline_main <- bind_rows(bline_main %>% mutate(value = as.integer(value)),
                        bline_mismatch)
bline_main <- bline_main %>% 
  select(nct_id, arm_name, n_participants = value)
rslts_prcnt %>% 
  anti_join(bline_main)

rslts_prcnt <- rslts_prcnt %>% 
  inner_join(bline_main) %>% 
  mutate(value = as.double(value),
         value = value * n_participants/100,
         measure_type = "n")

## 215 cleaned, as expected, as one is missing
rslts_clean <- bind_rows(rslts_n %>% 
                           select(nct_id, arm_name, measure,
                                            value) %>% 
                           mutate(value = as.integer(value)),
                         rslts_nprcnt %>% select(nct_id, arm_name, measure,
                                                value),
                         rslts_prcnt) %>% select(nct_id, arm_name, measure,
                                                 value)
ae <- rslts_clean %>% 
  filter(measure == "AE") %>% 
  select(nct_id, arm_name, ae = value) 
sae <- rslts_clean %>% 
  filter(measure == "SAE") %>% 
  select(nct_id, arm_name, sae = value) 

ae_sae <- bind_rows(ae %>% select(nct_id, arm_name),
                    sae %>% select(nct_id, arm_name)) %>% 
  distinct() %>% 
  left_join(ae) %>% 
  left_join(sae) %>% 
  arrange(nct_id, arm_name)

## Only 3 trials without sae data, but with AE data
no_sae <- ae_sae %>% 
  group_by(nct_id) %>% 
  mutate(allmis = all(is.na(sae))) %>% 
  ungroup() %>% 
  filter(allmis) %>% 
  distinct(nct_id, arm_name) %>% 
  group_by(nct_id) %>% 
  mutate(arm = seq_along(arm_name)) %>% 
  spread(arm, arm_name, fill = "")

## Note one trial collapses arms for AE analysis (eg Total valsartan rather valsartan at different doses)
ae_sae %>% filter(nct_id == "NCT00409760")

## Folow-up
## Two trials with missingess on follow-up, take from CSDR
rslts %>% filter(is.na(`follow-up`)) %>% distinct(nct_id)
fup <- rslts %>% 
  mutate(fup_units = str_extract(`follow-up`, "[a-z]{1,}"),
         fup_x = str_extract(`follow-up`, "[0-9]{1,}") %>% as.double()) %>% 
  distinct(nct_id, `follow-up`, fup_units, fup_x) %>% 
  mutate(fu_days = case_when(
    fup_units == "days" ~ fup_x,
    fup_units == "weeks" ~ fup_x * 7,
    fup_units == "months" ~ fup_x * (365.25/12),
    fup_units == "year" ~ fup_x * (365.25)
  )) %>% 
  select(nct_id, fu_days)
## All with results are have follow-up data
anti_join(rslts_clean, fup)

## Now clean baseline data ----
base <- shts$baseline_extract
## fix error w

base <- base %>% 
  filter(!measure == "continued on step 1 treatment")

## Drop ethnicity, only one trial with a measure
## Drop weight as no height so no use for bmi
## keep age median and range as is for now, maybe convert later
## similarly with BMI <25
base <- base %>% 
  filter(!measure %in% c("ethnicity", "weight kg")) 

## 
nmbr <- base %>% 
  filter(measure_type == "n") %>% 
  mutate(baseline_extract_comment1 = case_when(
    str_detect(value, "anticipated") ~ "anticipated", 
    str_detect(value, "Approximately") ~ "anticipated", 
    TRUE ~ baseline_extract_comment1),
         measure_type = "nmbr",
         value = parse_double(value %>% str_extract("[0-9]{1,}")),
    measure = case_when(measure == "male" ~ "male",
                        is.na(measure) ~ "participants",
                        TRUE ~ measure))

mn <- base %>% 
  filter(measure_type == "mean") %>% 
  mutate(value = parse_double(value))

s <- base %>% 
  filter(measure_type == "sd") %>% 
  mutate(measure_type = "s",
         value = parse_double(value))

mean_sd <- base %>% filter(measure_type %>% str_detect("mean"),
               measure_type %>% str_detect("sd")) %>% 
  separate(value, into = c("m", "s"), sep = "\\(") %>% 
  mutate(s = str_replace(s, "\\)", "") %>% str_trim) %>% 
  mutate_at(vars(m, s), parse_double) %>% 
  rename(mean = m, sd = s) %>% 
  select(-measure_type) %>% 
  gather("measure_type", "value", mean, sd)


n_prcnt <- base %>% filter(measure_type %>% str_detect("n"),
                           measure_type %>% str_detect("\\%")) %>% 
  separate(value, into = c("nmbr", "prcnt"), sep = "\\(") %>% 
  mutate(prcnt = str_replace(prcnt, "\\)", "") %>% str_trim) %>% 
  mutate_at(vars(nmbr, prcnt), parse_double) %>% 
  mutate(measure = case_when(
    str_detect(measure_type, "female") ~ "female",
    str_detect(measure, "female") ~ "female",
    str_detect(measure_type, "male") ~ "male",
    str_detect(measure, "male") ~ "male",
    TRUE ~ measure)) %>% 
  select(-measure_type) %>% 
  gather("measure_type", "value", nmbr, prcnt)

mf <-  base %>% 
  filter(measure_type %in% c("M:F", "M/F")) %>% 
  mutate(measure_type = if_else(str_detect(value, "\\%"), "prcnt", "nmbr"),
         value = value %>% 
           str_replace("M\\:", "") %>% 
           str_replace("F\\:", ":") %>% 
           str_replace_all("\\%", "")) %>% 
  separate(value, into = c("male", "female"), sep = "\\:|\\/") %>% 
  mutate_at(vars(male, female), parse_double) %>% 
  gather("measure", "value", male, female)

fm <-  base %>% 
  filter(measure_type %in% c("F:M")) %>% 
  mutate(measure_type = "nmbr") %>% 
  separate(value, into = c("female", "male"), sep = "\\:|\\/") %>% 
  mutate_at(vars(male, female), parse_double) %>% 
  gather("measure", "value", male, female)

fm_ratio <-  base %>% 
  filter(measure_type %in% c("F:M ratio")) %>% 
  mutate(measure_type = "prcnt") %>% 
  separate(value, into = c("female", "male"), sep = "\\:|\\/") %>% 
  mutate_at(vars(male, female), parse_double) %>% 
  mutate(value = male/(female + male),
         measure = "male") %>% 
  select(-female, -male)

mf_gen <- base %>% 
  filter(measure_type %in% c("males %", "males n", "male percentage", "male", "female")) %>% 
  mutate(measure = case_when(
    measure_type %in% c("males %", "males n", "male percentage", "male") ~ "male",
    measure_type %in% c("female") ~ "female",
    TRUE ~ NA_character_),
    measure_type = if_else(str_detect(measure_type, "\\%|percent") | str_detect(value, "\\%|percent"),
                           "prcnt",
                           "nmbr"),
    value = parse_double(value %>% str_replace_all("\\%|\\:", "")))

age_bands <- base %>% 
  filter(!is.na(value)) %>% 
  filter(measure_type %in% c("age bands",  "age group")) %>% 
  separate(value, into = c("band1", "band2"), sep = "\\,") %>% 
  gather("bands", "value", band1, band2) %>% 
  separate(value, into = c("band", "value"), sep = "\\:") %>% 
  mutate(measure_type = paste0("age ", str_trim(band)),
         value = parse_double(value)) %>% 
  select(nct_id, arm_name_definitive, measure, measure_type, value, population, starts_with("baseline"))

age_se <- base %>% 
  filter(measure == "age", measure_type == "se") %>% 
  left_join(nmbr %>% select (nct_id, arm_name_definitive, nmbr = value)) %>% 
  mutate(value = parse_number(value) * nmbr^0.5,
         measure_type = "s") %>% 
  select(-nmbr)

age_range <- base %>% 
  filter(measure == "age", measure_type %in% c("range")) %>% 
  separate(value, into = c("min_range", "max_range"), sep = "\\-") %>% 
  gather("measure_type", "value", min_range, max_range) %>% 
  mutate(value = str_extract(value, "[0-9]{1,}") %>% parse_double())

base <- base %>% 
  filter(!(measure_type %>% str_detect("mean") &
             measure_type %>% str_detect("sd")),
         !(measure_type %>% str_detect("n") &
             measure_type %>% str_detect("\\%")),
         !measure_type %in% c("n", "mean", "sd", "M:F", "M/F","F:M","F:M ratio",
                              "males %", "males n", "male percentage", "male", "female",
                              "age bands",  "age group"),
         !(measure == "age" & measure_type %in% c("se", "range")))


base %>% 
  count(measure, measure_type) %>% 
  arrange(measure, desc(n))%>% 
  as.data.frame()

## Drop ethnicity, only one trial with a measure
## Drop weight as no height so no use for bmi
## keep age median and range as is for now, maybe convert later
## similarly with BMI <25
base_remain <- base %>% 
  mutate(value = parse_double(value))


base_clean <- bind_rows(nmbr, mn, s, mean_sd, n_prcnt, mf, fm, fm_ratio, mf_gen, 
                        age_bands, age_se, age_range, base_remain, .id = "tlb_source")
base_clean <- bind_rows(base_clean) %>% 
  mutate(measure = str_to_lower(measure))

## As expected in terms of exclusions
xmn <- anti_join(shts$baseline_extract,
          base_clean %>% select(nct_id, arm_name_definitive, measure))

base_clean %>% 
  group_by(measure) %>% 
  summarise(trials = sum(!duplicated(nct_id)))

rm(nmbr, mn, s, mean_sd, n_prcnt, mf, fm, fm_ratio, mf_gen, 
   age_bands, age_se, age_range, base_remain, xmn)

saveRDS(list(base = base_clean, ae_sae = ae_sae, follow_up = fup), "Scratch_data/cleaned_guy_neave.Rds")
