#05_analysis.R
library(tidyverse)
library(rstanarm)
library(broom.mixed)

dfs <- readRDS("Processed_data/age_sex_bmi_ae_sae2.Rds")
list2env(dfs, envir = .GlobalEnv)
rm(dfs)
expected <- readRDS("data/SAE_ratio_observed_expected.Rds") 

## Restrict all trials ----
trials <- trials %>% 
  filter(!type_comparison %in% c("same_class7", "diff_class7"),
         !phase == "Phase 2/Phase 3") 

tots <- tots %>% 
  filter(!is.na(sae)) %>% 
  inner_join(trials %>% select(nct_id, aliskiren, minimum_age, maximum_age, hard_outcome, type_comparison, fu_days, phase))
tots$subjects_at_risk[tots$nct_id == "NCT00553267"] <- 947
tots <- tots %>% 
  mutate(pt = fu_days * (subjects_at_risk-sae) + 0.5 * sae*fu_days,
         older = as.integer(minimum_age >=60),
         rate = 1000*sae/pt)
## note many missing expected because has no mean age or sd age
tots <- tots %>% 
  left_join(expected)

## run regression model, unadjusted no effect ----
library(rstanarm)
unad <- glm(sae ~ older + offset(log(pt)), data = tots, family = "poisson")
summary(unad)
unad_stan <- stan_glmer(sae ~ older + offset(log(pt)) + (1|nct_id), data = tots, family = "poisson")
summary(unad_stan)

tapply(1000*tots$sae/tots$pt, tots$older, identity)

## Adjusted model, same association as before in MPH analysis, 2.22 fold higher rate ----
mod2 <- update(unad, . ~ . + aliskiren + hard_outcome + type_comparison + phase)
summary(mod2)
mod2_stan <- update(unad_stan, . ~ . + aliskiren + hard_outcome + type_comparison + phase)
summary(mod2_stan)
mod2_stan_no_hard <- update(unad_stan, . ~ . + aliskiren + type_comparison + phase,
                            data = tots %>% filter(hard_outcome ==0))
mod3_stan <- update(unad_stan, . ~ . + aliskiren + type_comparison + phase)
posterior_interval(mod2_stan_no_hard, 0.95)
ExportRes <- function(modelname){
  # browser()
  a <- broom.mixed::tidy(modelname)
  b <- posterior_interval(modelname, 0.95) %>% 
    as_tibble(rownames = "term")
  a <- a %>% 
    inner_join(b)
  a <- a %>% 
    mutate_at(vars(estimate, `2.5%`, `97.5%`), function(x) x %>% 
                exp() %>% 
                formatC(digits = 2, format = "f", flag = "0")) %>% 
    mutate(res = paste0(estimate, " (", `2.5%`, "-",`97.5%`, ")"))
  a %>% 
    select(term, res)
}
unad <- ExportRes(unad_stan)
minad <- ExportRes(mod3_stan)
adj <- ExportRes(mod2_stan)

## Run model with expected count as offset to estimate ratio ----
tots <- tots %>% 
  mutate(expected_count = rate_mean * pt/1000,
         ssaer = sae/expected_count)

## plot expected counts on same plot as tots
tots_lng <- tots %>% 
  gather("obs_exp", "point", rate, rate_mean)

plot2 <- ggplot(tots, aes(x = age_m, y = rate, size = pt, colour = factor(older), 
                          shape = type_comparison,
                          alpha = factor(1- hard_outcome, labels = c("hard", "soft")))) + 
  geom_point() +
  geom_point(mapping= aes(y = rate_mean), colour = "black", size = 1) +
  # geom_linerange(mapping= aes(ymin = rate_lci, ymax = rate_uci), colour = "black", size = 1) +
  geom_smooth(mapping= aes(y = rate_mean, group = 1), colour = "grey", se = FALSE) +
  # facet_grid( ~phase) +
  scale_size(guide = FALSE)
plot2

mod0_ratio <- stan_glmer(sae ~ older + offset(log(expected_count)) + (1|nct_id), data = tots, family = "poisson")
summary(mod0_ratio)
mod1_ratio <- update(mod0_ratio, . ~ . + aliskiren + hard_outcome + type_comparison + phase, data = tots, family = "poisson")
summary(mod1_ratio)

unad_ratio <- ExportRes(mod0_ratio)
adj_ratio  <- ExportRes(mod1_ratio)

for_tbls <- bind_rows(rate_unad = unad,
                      rate_adj = adj,
                      ratio_unad = unad_ratio,
                      ratio_adj = adj_ratio, .id = "tbl_type")
for_tbls <- for_tbls %>% 
  separate(tbl_type, into = c("measure", "model"))
for_tbls %>% 
  spread(measure, res)  %>% 
  arrange(desc(model)) %>% 
  filter(term == "older")

## Run sensitivity analysis dropping each variable in turn
# sa1 <- map(tots$nct_id, ~ ExportRes(update(mod2_stan, data = tots %>% filter(!nct_id == .x))))
# saveRDS(sa1, "Scratch_data/sensitivity_analysis.Rds")
sa1 <- readRDS("Scratch_data/sensitivity_analysis.Rds")
names(sa1) <- tots$nct_id
# effect evident even when leave each out
sa1_smry <- bind_rows(sa1, .id = "nct_id") %>% 
  filter(term == "older") %>% 
  inner_join(tots %>% select(nct_id, older))

# sa2 <- map(tots$nct_id, ~ ExportRes(update(mod1_ratio, data = tots %>% filter(!nct_id == .x))))
# saveRDS(sa2, "Scratch_data/sensitivity_analysis_ratio.Rds")
sa2 <- readRDS("Scratch_data/sensitivity_analysis_ratio.Rds")
names(sa2) <- tots$nct_id
# effect evident even when leave each out
sa2_smry <- bind_rows(sa2, .id = "nct_id") %>% 
  filter(term == "older") %>% 
  inner_join(tots %>% select(nct_id, older))
