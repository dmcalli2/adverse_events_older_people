#05_analysis.R
library(tidyverse)

dfs <- readRDS("Processed_data/age_sex_bmi_ae_sae.Rds")
list2env(dfs, envir = .GlobalEnv)
rm(dfs)
expected <- readRDS("data/SAE_ratio_observed_expected.Rds") 
expected <- expected %>% 
  select(nct_id, expected_rate = expected_count, ratio, ratio_point, ratio_lower, ratio_upper)

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
## note many missing because has no mean age or sd age
tots <- tots %>% 
  left_join(expected)


## run regression model, unadjusted no effect ----
unad <- glm(sae ~ older + offset(log(pt)), data = tots, family = "poisson")
summary(unad)
tapply(1000*tots$sae/tots$pt, tots$older, identity)

## Adjusted model, same association as before in MPH analysis, 2.22 fold higher rate ----
mod2 <- update(unad, . ~ . + aliskiren + hard_outcome + type_comparison + phase)
summary(mod2)
a <- broom::tidy(mod2)
a$rr <- exp(a$estimate)
a$rr[1] <- NA_real_
## do plot
plot1 <- ggplot(tots, aes(x = age_m, y = rate, size = pt, colour = factor(older), shape = type_comparison)) + 
  geom_point(alpha = 0.8) +
  facet_grid(factor(hard_outcome, labels = c("soft", "hard")) ~phase)
plot1

## Adjust for age, simple imputation first, later consider more complex ----
tots <- tots %>% 
  mutate(age_imp = if_else(is.na(age_m), mean(age_m, na.rm = TRUE), age_m))

mod3 <- update(mod2, . ~ . + age_imp)
summary(mod3)

## mod exploratory analysis, it is the "hard outcome" adjustment which makes the difference ----
mod2_1 <- update(unad, . ~ . + aliskiren )
summary(mod2_1)
mod2_2 <- update(unad, . ~ . + hard_outcome )
summary(mod2_2)
mod2_3 <- update(unad, . ~ . + type_comparison  )
summary(mod2_3)
mod2_4 <- update(unad, . ~ . + phase  )
summary(mod2_4)

## do restricted model excluding hard outcome trials, similar associations ----
mod2_5 <- update(mod2, . ~ . - hard_outcome, data = tots %>% filter(!hard_outcome==1))
summary(mod2_5)

## model association with age within older and non older ----
model_age_again <- glm(sae ~ age_m + phase + offset(log(pt)), data = tots %>% filter(older == 0, hard_outcome == 0), family = "poisson")
summary(model_age_again)

## Run model with uncertainty in covariate ----
tots <- tots %>% 
  mutate(expected_count = expected_rate * pt/1000)

unad_e <- glm(sae ~ older + offset(log(expected_count)), data = tots, family = "poisson")
summary(unad_e)
mod2_e <- update(unad_e, . ~ . + aliskiren + hard_outcome + type_comparison + phase)
summary(mod2_e)
summary(mod2)

## get standard errors
tots <- tots %>% 
  mutate(se1 = ratio_point-ratio_lower,
         se2 = ratio_upper-ratio_point,
         se_other = log(ratio_upper/ratio_lower)/(2*1.96))



## do plot
plot1 <- ggplot(tots, aes(x = age_m, y = ratio, ymin = ratio_lower, ymax = ratio_upper, colour = factor(older), shape = type_comparison)) + 
  geom_point(alpha = 0.8) +
  geom_errorbar() +
  facet_grid(factor(hard_outcome, labels = c("soft", "hard")) ~phase) 
plot1


library(runjags)



modelstring <- "
"

