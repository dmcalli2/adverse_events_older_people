library(tidyverse)
library(msm)
library(truncnorm)

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

coefs <- read.csv("SAIL_data/age_18_plus_fp_coef.csv")
vcov <- read.csv("SAIL_data/age_18_plus_fp_vcov.csv")

vcov <- subset(vcov, select = c(2:5))
names(vcov) <- coefs$X
vcov <- as.matrix(vcov)

predict_fn <- function(coefs, vcov, age, GNDR_CD, output){
  
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
  else{
    NA
  }
}



dfs <- readRDS("Processed_data/age_sex_bmi_ae_sae.Rds")
list2env(dfs, envir = .GlobalEnv)
rm(dfs)

## Restrict all trials
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

tots2 <- tots%>%
  filter(!is.na(male), 
         !is.na(female), 
         !is.na(age_m), 
         !is.na(age_s))
         

get_parameters <- function(x){
  lower <- tots2$minimum_age[x]
  upper <- tots2$maximum_age[x]
  upper <- ifelse(is.na(upper), 100L, upper)
  trial_mean <- tots2$age_m[x]
  trial_sd <- tots2$age_s[x]
  trial_var <- tots2$age_s[x]^2
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
  list(estimate = estimate, freqs = freqs)
  estimate%>%
    slice(1L)
  estimate
}


test_fn <- function(x){
  
  
  param <- get_parameters(x)
  min_age_trial <- ifelse(is.na(tots2$minimum_age[x]), -Inf, tots2$minimum_age[x])
  max_age_trial <- ifelse(is.na(tots2$maximum_age[x]), Inf, tots2$maximum_age[x])
  
  
  males <- data.frame( age = rtruncnorm(((tots2$male[x] / tots2$subjects_at_risk[x])*10000), 
                                        a = min_age_trial, b = max_age_trial,  mean = param$mu_x[1], sd = param$sd_x[1]),
                       GNDR_CD = 1)
  females <- data.frame( age = rtruncnorm(((tots2$female[x] / tots2$subjects_at_risk[x])*10000), 
                                          a = min_age_trial, b = max_age_trial, mean = param$mu_x[1], sd = param$sd_x[1]), 
                         GNDR_CD = 2)
  sample <- males%>%
    full_join(females)
  
  p <- map2(sample$age, sample$GNDR_CD, function(x, y) {predict_fn(coefs = coefs$model3.coefficients, 
                                                                   vcov = vcov, 
                                                                   age = x, 
                                                                   GNDR_CD = y, 
                                                                   output = "point")})
  sample$estimate <- unlist(p)
  sample2 <-  sum(sample$estimate)/nrow(sample)
  sample2 
}

test_fn(1)
expected_count <- map(1:nrow(tots2), test_fn)

tots2$expected_count <- unlist(expected_count)
tots2$ratio <- tots2$rate/tots2$expected_count


uncertainty_fn <- function(x, output){
  
  
  param <- get_parameters(x)
  min_age_trial <- ifelse(is.na(tots2$minimum_age[x]), -Inf, tots2$minimum_age[x])
  max_age_trial <- ifelse(is.na(tots2$maximum_age[x]), Inf, tots2$maximum_age[x])
  
  
  males <- data.frame( age = rtruncnorm(((tots2$male[x] / tots2$subjects_at_risk[x])*10000), 
                                        a = min_age_trial, b = max_age_trial,  mean = param$mu_x[1], sd = param$sd_x[1]),
                       GNDR_CD = 1)
  females <- data.frame( age = rtruncnorm(((tots2$female[x] / tots2$subjects_at_risk[x])*10000), 
                                          a = min_age_trial, b = max_age_trial, mean = param$mu_x[1], sd = param$sd_x[1]), 
                         GNDR_CD = 2)
  sample <- males%>%
    full_join(females)
  
  p <- map2(sample$age, sample$GNDR_CD, function(x, y) {predict_fn(coefs = coefs$model3.coefficients, 
                                                                   vcov = vcov, 
                                                                   age = x, 
                                                                   GNDR_CD = y, 
                                                                   output = "point")})
  sample$estimate <- unlist(p)
  
  sample$observed <- rpois(nrow(sample), tots2$sae[x])
  sample$observed <- sample$observed/tots2$pt[x]*1000
  
  sample2 <-  sum(sample$estimate)/nrow(sample)
  
  ratio <- sample$observed/sample$estimate
  
  if(output=="point"){
    quantile(ratio, probs = 0.5)
  }
  else if(output=="lower"){
    quantile(ratio, probs = 0.025)
    
  }
  else if (output=="upper"){
    quantile(ratio, probs = 0.975)
    
  }
  else{
    NA
  }
}

uncertainty_fn(1:2, output = "point")

ratio_point <- map2(1:nrow(tots2),"point", uncertainty_fn)
 
ratio_lower <- map2(1:nrow(tots2),"lower", uncertainty_fn)
ratio_upper <- map2(1:nrow(tots2),"upper", uncertainty_fn)

tots2$ratio_point <- unlist(ratio_point)
tots2$ratio_lower <- unlist(ratio_lower)
tots2$ratio_upper <- unlist(ratio_upper)
rm(ratio_point)
rm(ratio_lower)
rm(ratio_upper)
class(tots2$ratio_upper)

tots2$v_line <- 1

tots2$label <- ifelse(tots2$older==1, "Older-people trials", "Standard trials")

ggplot(tots2, aes(x = nct_id, y = ratio_point,  colour = phase))+
  geom_point()+
  geom_pointrange(aes(ymin = ratio_lower, ymax = ratio_upper))+
  facet_grid(label~., scales = "free", space = "free")+
  geom_hline(aes(yintercept = v_line))+
  coord_flip()+
  ggtitle("Ratio of observed to expected counts")+
  ylab("Observed count / Expected count")+
  xlab("Trial ID")

saveRDS(tots2, "Data/SAE_ratio_observed_expected.csv")

# get_parameters <- function(x){
#   lower <- tots$minimum_age[x]
#   upper <- tots$maximum_age[x]
#   upper <- ifelse(is.na(upper), 100L, upper)
#   trial_mean <- tots$age_m[x]
#   trial_sd <- tots$age_s[x]
#   trial_var <- tots$age_s[x]^2
#   # 
#   #   ## Create grid
#   mu_x <- seq(lower, upper, 0.05)
#   sd_x <- seq(1, upper-lower, 0.05)
#   full_grid <- expand.grid(mu_x = mu_x, sd_x = sd_x)
#   # 
#   #   ## Calculate for all values of grid, is vectorised so is fast
#   full_grid$mean_x <- mean.tnorm(full_grid$mu_x, full_grid$sd_x, lower, upper)
#   full_grid$var_x <- var.tnorm(full_grid$mu_x, full_grid$sd_x, lower, upper)
#   # print(nrow(full_grid))
#   # 
#   #   ## Identify closest values
#   full_grid <- full_grid %>%
#     as_tibble() %>%
#     mutate(mean_diff = abs(trial_mean - mean_x),
#            var_diff = abs(trial_var - var_x),
#            total_diff = mean_diff + var_diff) %>%
#     arrange(total_diff, mean_diff, var_diff)
#   #   ## Append original parameters
#   estimate <- full_grid %>%
#     slice(1:10) %>%
#     mutate(trial_mean = trial_mean,
#            trial_var = trial_var,
#            trial_lower = lower,
#            trial_upper = upper,
#            trial_sd = trial_sd) %>%
#     select(trial_mean, mean_x, trial_var, var_x, mu_x, sd_x, trial_sd, everything())
#   #   # Produce estimates for distribution
#   freqs <- msm::ptnorm(seq(0, 100, 1),
#                        estimate$mu_x[1], estimate$sd_x[1],
#                        estimate$trial_lower[1], estimate$trial_upper[1])
#   list(estimate = estimate, freqs = freqs)
#   estimate%>%
#     slice(1L)
#   estimate
# }
# 
# param <- get_parameters(1)
# 
# 
# 
# 
# males <- data.frame( age = rtruncnorm(((tots$male[1] / tots$subjects_at_risk[1])*10000), 
#                                       a = 45, mean = param$mu_x[1], sd = param$sd_x[1]),
#                GNDR_CD = 1)
# females <- data.frame( age = rtruncnorm(((tots$female[1] / tots$subjects_at_risk[1])*10000), 
#                                    a = 45, mean = param$mu_x[1], sd = param$sd_x[1]), 
#                GNDR_CD = 1)
# sample <- males%>%
#   full_join(females)
# 
# 
# 
# 
# p <- map2(sample$age, sample$GNDR_CD, function(x, y) {predict_fn(coefs = coefs$model3.coefficients, 
#                                                                 vcov = vcov, 
#                                                                 age = x, 
#                                                                 GNDR_CD = y, 
#                                                                 output = "point")})
# l <- map2(sample$age, sample$GNDR_CD, function(x, y) {predict_fn(coefs = coefs$model3.coefficients, 
#                                                                  vcov = vcov, 
#                                                                  age = x, 
#                                                                  GNDR_CD = y, 
#                                                                 output = "lower")})
# u <- map2(sample$age, sample$GNDR_CD, function(x, y) {predict_fn(coefs = coefs$model3.coefficients, 
#                                                                  vcov = vcov, 
#                                                                  age = x, 
#                                                                  GNDR_CD = y, 
#                                                                 output = "upper")})
# sample$estimate <- unlist(p)
# sample$lci <- unlist(l)
# sample$uci <- unlist(u)
# 
# sample2 <- sample%>%
#   summarise(estimate = sum(estimate)/nrow(sample), 
#             lci = sum(lci)/nrow(sample),
#             uci = sum(uci)/nrow(sample))
# 
# ((tots$sae/tots$subjects_at_risk)/ (tots$fu_days/365))
# tots$sae[1]/tots$pt[1]*1000
# tots$subjects_at_risk[1]*tots$fu_days[1]
