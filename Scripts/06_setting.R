#6_setting

## put your username and password if want to run again
library(RPostgreSQL)
library(tidyverse)

drv <- dbDriver('PostgreSQL')
con <- dbConnect(drv, dbname="aact",host="aact-db.ctti-clinicaltrials.org", port=5432, user="dmcalli", password=mypass)
dbListTables(con)

nct_ids <- read_tsv("Processed_data/tots2.tsv")
nct_ids <- nct_ids$nct_id
nct_ids <- paste(nct_ids %>% unique(), collapse = "', '")
nct_ids <- paste0("('", nct_ids, "')")

countries <- dbGetQuery(con, paste0("SELECT * FROM countries WHERE nct_id IN ", nct_ids))
countries2 <- countries %>% 
  as_tibble() %>% 
  filter(is.na(removed) | removed == FALSE) %>% 
  select(-removed) %>% 
  distinct(nct_id, name)


nct_ids <- read_tsv("Processed_data/tots2.tsv") %>% 
  distinct(nct_id) %>% 
  left_join(countries2) %>% 
  mutate(name = if_else(is.na(name), "unknown", name)) %>% 
  rename(country = name)

write_tsv(nct_ids, "Processed_data/countries_each_trial.tsv")
