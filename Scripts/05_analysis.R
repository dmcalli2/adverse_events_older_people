#05_analysis.R
library(tidyverse)

dfs <- readRDS("Processed_data/age_sex_bmi_ae_sae.Rds")
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
## note many missing because has no mean age or sd age
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

## Run model with expected count as offset ----
tots <- tots %>% 
  mutate(expected_count = rate_mean * pt/1000,
         ssaer = sae/expected_count)

unad_e <- glm(sae ~ older + offset(log(expected_count)), data = tots, family = "poisson")
summary(unad_e)
mod2_e <- update(unad_e, . ~ . + aliskiren + hard_outcome + type_comparison + phase)
summary(mod2_e)
summary(mod2)


## do plot
plot1 <- ggplot(tots, aes(x = age_m, y = ssaer, colour = factor(older), shape = type_comparison)) + 
  geom_point(alpha = 0.8) +
  facet_grid(factor(hard_outcome, labels = c("soft", "hard")) ~phase) 
plot1

library(rstanarm)
mod0_ratio <- stan_glmer(sae ~ older + offset(log(expected_count)) + (1|nct_id), data = tots, family = "poisson")
summary(mod0_ratio)
mod1_ratio <- update(mod0_ratio, . ~ . + aliskiren + hard_outcome + type_comparison + phase, data = tots, family = "poisson")
summary(mod1_ratio)
mod2_ratio <- update(mod1_ratio, . ~ . + age_m, data = tots)
summary(mod2_ratio)
mod3_ratio <- update(mod2_ratio, . ~ . + age_m + I(age_m^2), data = tots)
summary(mod3_ratio)
# launch_shinystan(mod2_ratio)

tots_new <- tots %>% 
  filter(!is.na(expected_count))

## do plot
sims <- as.matrix(mod2_ratio)
sims <- sims[, 
             c('(Intercept)', 'aliskiren', 'hard_outcome', 'type_comparisondiff_class3', 'type_comparisondiff_class5', 'phasePhase 4', 'age_m', 'older')]
tots_new2 <- tots_new %>% 
  select(age_m, aliskiren, hard_outcome, type_comparison, phase, older) %>% 
  group_by(aliskiren, hard_outcome, type_comparison, phase, older) %>% 
  nest() %>% 
  ungroup()
tots_new2$data <- map(tots_new2$data, ~ tibble(age_m = min(.x$age_m):max(.x$age_m))) 
tots_new2 <- tots_new2 %>% 
  unnest(data)
sim_pred <- model.matrix(~ aliskiren + hard_outcome + type_comparison + phase + age_m + older, tots_new2)
sim_pred <- sim_pred[,colnames(sims)]
res <- sims %*% t(sim_pred)
res[] <- exp(res)
res_m <- colMeans(res)
res_lci <- apply(res, 2, quantile, probs = c(0.025))
res_uci <- apply(res, 2, quantile, probs = c(0.975))

tots_new2 <- tots_new2 %>% 
  mutate(ratio_point = res_m,
         ratio_lower = res_lci,
         ratio_upper = res_uci)

plot1 <- ggplot(tots_new2, aes(x = age_m, 
                               y = ratio_point, 
                               ymax = ratio_upper,
                               ymin = ratio_lower, 
                                colour = type_comparison, 
                                linetype = factor(aliskiren), 
                                shape = factor(aliskiren), 
                                fill = type_comparison,
                                group = interaction(type_comparison, older, aliskiren))) + 
  geom_ribbon(alpha = 0.1, colour = NA) +
  geom_line() +
  geom_point() +
  facet_grid(factor(hard_outcome, labels = c("soft", "hard")) ~phase) +
  scale_y_log10()
plot1

summary(mod1_ratio)
