library(tidyverse)
# library(msm)
library(truncnorm)

## data ----
coefs <- read_csv("SAIL_data/age_18_plus_fp_coef.csv")
vcov <- read_csv("SAIL_data/age_18_plus_fp_vcov.csv")
vcov <- subset(vcov, select = c(2:5))
names(vcov) <- coefs$X1
vcov <- as.matrix(vcov)


dfs <- readRDS("Processed_data/age_sex_bmi_ae_sae.Rds")
list2env(dfs, envir = .GlobalEnv)
rm(dfs)

## Functions ----
mean.tnorm<-function(mu,sd,lower,upper){
  ##return the expectation of a truncated normal distribution
  lower.std=(lower-mu)/sd
  upper.std=(upper-mu)/sd
  mean=mu+sd*(dnorm(lower.std)-dnorm(upper.std))/
    (pnorm(upper.std)-pnorm(lower.std))
  return(mean)
}
var.tnorm<-function(mu,sd,lower,upper){
  ##return the variance of a truncated normal distribution
  lower.std=(lower-mu)/sd
  upper.std=(upper-mu)/sd
  variance=sd^2*(1+(lower.std*dnorm(lower.std)-upper.std*dnorm(upper.std))/
                   (pnorm(upper.std)-pnorm(lower.std))-((dnorm(lower.std)-dnorm(upper.std))/
                                                          (pnorm(upper.std)-pnorm(lower.std)))^2)
  return(variance)
}

predict_fn <- function(coefs, vcov, age, GNDR_CD, output){
  ## calculate effect for a set of covariates
  cov <- c(1, I((age/100)^0.5), I((age/100)^2), GNDR_CD)
  lp <- coefs[1] + I((age/100)^0.5)*coefs[2]+ I((age/100)^2)*coefs[3] + GNDR_CD*coefs[4]
  
  se <- c(t(cov) %*% vcov %*% cov)
  
  
  pred <- exp(lp)*1000
  
  lci <- exp(lp-1.96*se)*1000
  uci <- exp(lp+1.96*se)*1000
  
  if(output=="point"){
    pred
  }
  else if(output=="lower"){
    lci
  }
  else if (output=="upper"){
    uci
  }
  else if (output == "verbose"){
    list(lp = lp, se = se, lci = lci, uci = uci)
  }
  else{
    NA
  }
}

get_parameters <- function(maxage, minage, meanage, sdage){
  lower <- minage
  upper <- maxage
  upper <- ifelse(is.na(upper), 100L, upper)
  trial_mean <- meanage
  trial_sd <- sdage
  trial_var <- sdage^2
  # 
  #   ## Create grid
  mu_x <- seq(lower, upper, 0.05)
  sd_x <- seq(1, upper-lower, 0.05)
  full_grid <- expand.grid(mu_x = mu_x, sd_x = sd_x)
  # 
  #   ## Calculate for all values of grid, is vectorised so is fast
  full_grid$mean_x <- mean.tnorm(full_grid$mu_x, full_grid$sd_x, lower, upper)
  full_grid$var_x <- var.tnorm(full_grid$mu_x, full_grid$sd_x, lower, upper)
  # print(nrow(full_grid))
  # 
  #   ## Identify closest values
  full_grid <- full_grid %>%
    as_tibble() %>%
    mutate(mean_diff = abs(trial_mean - mean_x),
           var_diff = abs(trial_var - var_x),
           total_diff = mean_diff + var_diff) %>%
    arrange(total_diff, mean_diff, var_diff)
  #   ## Append original parameters
  estimate <- full_grid %>%
    slice(1:10) %>%
    mutate(trial_mean = trial_mean,
           trial_var = trial_var,
           trial_lower = lower,
           trial_upper = upper,
           trial_sd = trial_sd) %>%
    select(trial_mean, mean_x, trial_var, var_x, mu_x, sd_x, trial_sd, everything())
  #   # Produce estimates for distribution
  freqs <- msm::ptnorm(seq(0, 100, 1),
                       estimate$mu_x[1], estimate$sd_x[1],
                       estimate$trial_lower[1], estimate$trial_upper[1])
  tibble(age = 0:100, age_p = freqs - lag(freqs, 1, default = 0))
}

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

names(tots)

tots <- tots %>%
  filter(!is.na(male), 
         !is.na(female), 
         !is.na(age_m), 
         !is.na(age_s))
         
# replace missing maximum ages with infinity. All have minimum ages
 tots <- tots %>% 
   mutate(maximum_age = if_else(is.na(maximum_age), Inf, as.double(maximum_age)))

# # get parameters of truncated normal distribution, quite slow so run only when have to ----
# tots$age_dist <- pmap(list(tots$age_m, tots$age_s, tots$minimum_age, tots$maximum_age),
#                     function(mean, sd, min, max) get_parameters(maxage = max, minage = min, meanage = mean, sdage = sd))
# saveRDS(tots, "Scratch_data/age_distributions.Rds")
tots <- readRDS("Scratch_data/age_distributions.Rds")

## limit data to minimum and maximum ages across all datasets then estimate rates for each of these ----
min_age <- map_dbl(tots$age_dist, ~ .x %>% filter(.x$age_p >=0.01) %>% pull(age) %>% min()) %>% min()
max_age <- map_dbl(tots$age_dist, ~ .x %>% filter(.x$age_p >0.01) %>% pull(age) %>% max()) %>% max()

mycovs <- tibble(age = min_age:max_age) %>% 
  mutate(age1 = (age/100)^0.5,
         age2 = (age/100)^2)

smpl_coefs <- mvtnorm::rmvnorm(10000, coefs$model3.coefficients, sigma = vcov)
colnames(smpl_coefs) <- coefs$X1

# estimate linear predictors for men and women
mycovs$men <- map2(mycovs$age1, mycovs$age2, function(age1, age2){
                     smpl_coefs[, "(Intercept)"] +
                     smpl_coefs[, "I((age/100)^0.5)"] *age1 +
                     smpl_coefs[, "I((age/100)^2)"] *age2 +
                     0})
mycovs$women <- map2(mycovs$age1, mycovs$age2, function(age1, age2){
  smpl_coefs[, "(Intercept)"] +
    smpl_coefs[, "I((age/100)^0.5)"] *age1 +
    smpl_coefs[, "I((age/100)^2)"] *age2 +
    smpl_coefs[, "GNDR_CD"]})
mycovs$both <- map2(mycovs$men, mycovs$women, ~ tibble(iter = 1:10000, men = .x, women = .y))
mycovs <- mycovs %>% 
  select(-men, -women)
mycovs <- mycovs %>% 
  unnest(cols = both)
mycovs <- mycovs %>% 
  arrange(iter)
# convert to rate per person-year
mycovs <- mycovs %>% 
  mutate_at(vars(men, women), function(x) exp(x)) 
mycovs <- mycovs %>% 
  select(age, iter, men, women)

tots$smpls <- map(tots$age_dist, ~ inner_join(.x, mycovs))
tots$smpls <- map(tots$smpls, ~ .x %>% 
                    group_by(iter) %>% 
                    summarise(men =  weighted.mean(men, age_p),
                              women = weighted.mean(women, age_p)) %>% 
                    ungroup())
tots$age_dist <- NULL
tots <- tots %>% 
  mutate(male_p = male/(male + female))
tots$smpls <- map2(tots$smpls, tots$male_p, ~ .x %>% 
                    mutate(both = .y*men + (1-.y) * women))

## sample from beta distribution and multiple by number at risk to get CI for obserVED SAE
tots$sae_smpls <- map2(tots$sae, tots$subjects_at_risk, ~ rbeta(10000, .x + 0.005, .y-.x+0.005) * .y)
## or sample from Poisson
tots$sae_smpls_pois <- map(tots$sae, ~ rpois(10000, .x))
# compare beta and Poisson
pois <- map(tots$sae_smpls_pois, ~ tibble(m = mean(.x), s = sd(.x))) %>% 
  bind_rows(.id = "trial")
beta <- map(tots$sae_smpls, ~ tibble(m_beta = mean(.x), s_beta = sd(.x)))%>% 
  bind_rows(.id = "trial")
cmpr <- bind_cols(pois, beta %>% select(-trial))
plot(cmpr$m, cmpr$m_beta)
abline(a = 0, b = 1)
plot(cmpr$s, cmpr$s_beta)
abline(a = 0, b = 1)
## Divide by person time to get rate
tots$sae_smpls <- map2(tots$sae_smpls, tots$pt, ~ 1000 * .x/.y)

## Produce ratios ----
tots$smpls <- map(tots$smpls, ~ .x$both*1000)
tots$ratios <- map2(tots$sae_smpls, tots$smpls, ~ .y/.x)

## Summarise rates

tots$rate_mean <- map_dbl(tots$smpls, mean)
tots$rate_se <- map_dbl(tots$smpls, sd)
tots$rate_lci <- map_dbl(tots$smpls, ~ quantile(.x, 0.025))
tots$rate_uci <- map_dbl(tots$smpls, ~ quantile(.x, 0.975))

## Summarise ratios 
tots$ratio_mean <- map_dbl(tots$ratios, mean)
tots$ratio_se <- map_dbl(tots$ratios, sd)
tots$ratio_lci <- map_dbl(tots$ratios, ~ quantile(.x, 0.025))
tots$ratio_uci <- map_dbl(tots$ratios, ~ quantile(.x, 0.975))

tots <- tots %>% 
  select(nct_id, subjects_at_risk, sae, fu_days, pt, rate, 
         rate_mean, rate_se, rate_lci, rate_uci, 
         ratio_mean, ratio_se, ratio_lci, ratio_uci, 
         older, phase)

tots$v_line <- 1

tots$label <- ifelse(tots$older==1, "Older-people trials", "Standard trials")

a <- ggplot(tots, aes(x = nct_id, y = ratio_mean,  colour = phase))+
  geom_point()+
  geom_pointrange(aes(ymin = ratio_lci, ymax = ratio_uci))+
  facet_grid(label~., scales = "free", space = "free")+
  geom_hline(aes(yintercept = v_line))+
  coord_flip()+
  ggtitle("Ratio of observed to expected counts")+
  ylab("Observed count / Expected count")+
  xlab("Trial ID")
a
a <- a +
  coord_flip(ylim = c(0.9, 15))
a
saveRDS(tots, "Data/SAE_ratio_observed_expected.Rds")

