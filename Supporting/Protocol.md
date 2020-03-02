# Analysis plan

# Original trial selection

As part of an MPH project, from a denominator set of trials taken from clinicaltrials.gov which we previously described (https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-019-1427-1) we reviewed all phase 3 and 4 clinical trials whose eligibility criteria that specified a minimum age for inclusion of 60 or above. We found that RAAS drugs were commonly studied in such clinic trials carried out in older people. Hypertension is also an important public health issue and therefore RAAS drugs were selected as the exemplar for comparing clinical trials in older and younger people. Of the trials that this yielded we manually searched titles and brief summaries to identify those which related to drugs acting on the RAAS to treat hypertension, either as a single agent or in combination with another drug. In all cases the medication in the treatment arm included a drug which acted on the RAAS and in some cases the comparator drug was also of this class. Drugs which were considered to be included were ACE inhibitors, angiotensin receptor blockers and renin inhibitors.

For 13 of these older people trials, the following drugs had been studied aliskiren, irbesartan, olmesartan, telmisartan, and valsartan. From the denominator set of trials we therefore searched for any trial including these drugs where the indication was also stated as hypertension. We excluded any trials with complex indications (eg hypertension and diabetes, hypertension in heart failure).

This yielded 144 trials including the initial 13 trials. For the initial MPH project we matched 13 older people trials (up to 1:3 matching) to 22 non-older people trials from this set, and extracted adverse event and serious adverse event data. We found that the serious adverse event rate was around 2-fold higher in the older people trials than in the matched "general population" trials.

We subsequently decided to (a) obtain and compare serious adverse events data for the remaining set of trials to simplify the analytically complexity (by avoiding matching) and (b) obtain a comparison dataset from primary care to determine whether the ratio between SAEs in older people versus general population trials was larger, similar or smaller than that which would be expected based on the rates of hospitalisation and death in the community.

This document is designed to allow us to pre-specify these subsequent analyses. It was completed after the analysis completed for the MPH project, and after additional data was extracted from clincialtrials.gov, clinical study reports and papers and harmonised into a single dataset. It was however completed prior to any further plotting or analyses.

# Planned modelling

## Restriction
First we will restrict the non-older people trials on the following variables:- 

- comparison type - with placebo arm, different WHO ATC class at different levels of hierarchy (3-character, 5-character or 7-character) or same drug. 4 trials are against the same drug(s) and 7 are against drugs in the same class (to the 5th character) in the non-older people set but none of these are in the older people set, so exclude
- sex eligibility criteria - all trials allow both sexes
- hard outcome - both older and non-older people trials include heard outcomes
- phase - both include phase 3 and phase 4, the older trials do not include phase 2/3, so exclude this one trial

This leaves 111 trials, of which 11 are older people trials.

## Trials modelling
Having excluded these trials, in Poisson regression models, we will model serious adverse events on older people trial status adjusted for the following:-

1. Unadjusted, offset by estimated person time (calculated as  follow-up x (number of participants - 0.5 x number of SAEs))
2. 1 + direct renin inhibitors (y/n),  comparison type (placebo, different to 3-character, different to 5-character), phase (3 or 4), outcome type (hard or not)
3. 2 + age
4. 2 + sex
5. 2 + age + sex

All models will be fit in the Bayesian software package JAGGS, as this will allow us to include the standard deviation as well as mean for age as a covariate.

## Population modelling

We will estimate the age-sex specific incidence rate for first hospitalisation or death in the SAIL platform. From among our existing cohort of people with hypertension identified via primary care Read codes described previously (https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-019-1427-1) we will exclude any patient who has not received a treatment with any drug in the RAAS class, and those with recent myocardial infarction, stroke or TIA (to attempt to exclude patients started on these drugs for these indications rather than hypertension). The We will then model the event rate in Poisson regression models. In sensitivity analyses we will repeat the modelling without these exclusions.

## Estimating expected SAE

Since, in adults, hospitalisation and death is the major driver of serious adverse events, we will use these rates, age-sex distribution of each trials to calculate the expected SAE rate for each trial. We will then compare the ratio for the expected SAE rate for older versus younger people trials. Next, using the as well as the number of participants and follow-up time for each trial we will repeat the trials modelling, offsetting for this community-derived expected SAE rate, rather than the trial follow-up. This will allow us to determine whether any difference in the odler people versus younger people trials persists after accounting for the expected serious adverse events according to age and sex.

For the puposes of this analysis we will treat both the age-sex distribution of the trials and population incidence rates as fixed. From a separate project, we already know the sex-specific age distribution of each trial. We will assume that age has a truncated normal distribution (truncating at the age eligibility limits) and will estimate the paramaters of the truncated normal distribution for men and women. We will do so using the mean and standard deviation as well as the truncation limits for each trial, as well as the mean difference in age between men and women for hypertension trials (which we have for 5 trials from this set).

```
for(i in 1:n){
age_est[i] ~ dtnorm(m1, s1)T(min_age, max_age)
}
for(j in 1:n_women){
women[j] ~ dtnorm(m1 + m_diff, s2)T(min_age, max_age)
m_diff ~ dnorm(diff_mean, diff_se)
}
women_tot <- sum(women)
men_tot <- sum(men)/n_men
both_mean <- (women_tot + men_tot)/(n_men + n_women)


```



