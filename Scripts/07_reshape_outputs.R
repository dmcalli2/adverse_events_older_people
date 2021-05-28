library(tidyverse)

mydf <- readRDS("Processed_data/check_variation_with_data_changes.Rds")

QuickFx <- function(a){
  names(a) <- c("main_res", "unadsr", "adjsr")
  b <- bind_rows(unad = a$unadsr,
                 adj = a$adjsr, .id = "model")
  b <- b  %>% 
    mutate_at(vars(-term, -model), ~ round(1/.x,2)) %>% 
    rename(uci2 = lci,
           lci2 = uci) %>% 
    select(model, term, est, lci = lci2, uci = uci2) %>% 
    mutate(res = paste0(est, " (", lci, "-", uci, ")")) %>% 
    select(model, term, res) %>% 
    spread(term, res) %>% 
    select(model, stand_sr = cept, oldr_sr = oldr, rat_sr = dffr)
  
  a$main_res <- a$main_res %>% 
    select(model, irr = rate, sir_not_inv = ratio)
 a$main_res %>% 
   left_join(b)
  
}
mydf$res2 <- map(toloop$res, QuickFx)

mydf_final <- mydf %>% 
  select(-res) %>% 
  unnest(res2)

write_csv(mydf_final, "Outputs/examine_impact_small_differences_in_data.csv")
