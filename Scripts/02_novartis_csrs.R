library(tidyverse)

rvd <- read_csv("Data/reviewED_sae_complete.csv")

rvd <- rvd %>% filter(source == "Novartis") %>% 
  distinct(nct_id, other_id)

a <- list.files("../CSRs/novartis_excel/")

res <- map(rvd$other_id, ~ str_detect(a , .x %>% str_replace("E1$", "")))
rvd$matching <- map(res, ~ a[.x])
     
rvd <- rvd %>% 
  mutate(nmbr_csr = map_int(matching, length))
## One trial, note could not find any CSR for this trial on Novartis website
rvd %>% 
  filter(nmbr_csr ==0)

copyfiles <- rvd$matching %>% unlist()
file.copy(from = paste0("../CSRs/novartis_excel/", copyfiles), to = paste0("Data/Novartis_csrs/", copyfiles))

setdiff(copyfiles, list.files("Data/Novartis_csrs/"))
