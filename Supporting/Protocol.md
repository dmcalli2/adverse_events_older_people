Analysis plan for study of adverse event reporting in older people
trials
================

# Original trial selection

As part of an MPH project, from a denominator set of trials taken from
clinicaltrials.gov which we previously described
(<https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-019-1427-1>)
we reviewed all phase 3 and 4 clinical trials where the eligibility
criteria specified a minimum age for inclusion of 60 or above. We found
that renin-angiontensin-aldosterone system (RAAS) drugs were one of the
commoner drug classes in this set of trials. Since hypertension is also
an important public health issue we selected RAAS drugs as an exemplar
for comparing clinical trials in older and younger people.

Thriteen of these “older” trials were identified. These trialled
aliskiren, irbesartan, olmesartan, telmisartan, and/or valsartan. To
identify a comparator set of trials (hereafter termed non-older), we
searched in the denominator set of trials for any trial where one or
more these drugs was included among the interventions and where one or
more of the stated indications was hypertension. We identified 155
trials.

We then excluded any trials with additional indications (eg hypertension
and diabetes, hypertension and heart failure) excluding 11 trials,
leaving 144 trials including the initial 13 older people trials. For the
MPH project we matched these 13 older people trials (in one to many
matching) to 22 non-older people trials from this dataset and extracted
adverse event and serious adverse event (SAE) data. We found that the
adverse event rate was similar in both, but that serious adverse event
rate was around 2-fold higher in the older people trials than in the
matched “general population” trials.

This suggests that such older trials, at least in part, do acheive the
stated aim of improving trial representativeness (by including sicker
patients at higher risk of adverse events).

We subsequently decided to extract and compare serious adverse events
data for the remaining set of trials in order to simplify the analytical
complexity (by avoiding matching) and increase the precision of our
estimates. We also opted to obtain a comparison dataset from primary
care to determine whether the ratio between SAEs in older people versus
general population trials was larger, similar or smaller than that which
would be expected based on the rates of hospitalisation and death in the
community. This comparison was used because hospitalisation and death
are a major part of the definition of serious adverse events.

This document pre-specifies these two subsequent analyses. It was
written after completion of the MPH project, and after the additional
trial data was extracted and harmonised into a single dataset from
clincialtrials.gov, clinical study reports and published papers. It was
however completed prior to plotting the new data or conducting any
statistical analyses.

# Analyses

## Restriction of trials

Of the 144 trials, 22 had no adverse event reporting in
clinicaltrials.gov, published papers or accessible clinical study
reports. 11 of these were older trials and 133 were non-older trials. We
restricted the latter on the following variables to improve
comparability to the older trials:-

  - comparison type - trials which included a placebo arm, different WHO
    ATC class at different levels of hierarchy (3-character, 5-character
    or 7-character) or only a comparison against the same drug (eg
    different dosages). Among non-older trials, 4 compared the same drug
    (at different doses/regimes) and 7 compared drugs in the same
    5-character clas. There were no trials of this type in the older
    trials set, so these were excluded
  - sex eligibility criteria - all trials allow both sexes
  - hard outcome - both older and non-older trials include hard outcomes
  - phase - both older and non-older trials included phase 3 and phase
    4, the older trials do not include “phase 2/3”, so we excluded this
    one trial

This leaves 111 trials, of which 11 are older people trials.

## Regression models in trials

In Poisson regression models, we will model serious adverse events on
older people trial status with the following covariates:-

1.  Unadjusted, offset by estimated person time (calculated as follow-up
    x (number of participants - 0.5 x number of SAEs))
2.  1 + direct renin inhibitors (y/n), comparison type (placebo,
    different to 3-character, different to 5-character), phase (3 or 4),
    outcome type (hard or not)
3.  2 + age
4.  2 + sex
5.  2 + age + sex

The primary analysis will be the older people trial status variable from
model 2.

Although the time since a product has been first trialled may influence
the adverse event rate (eg because clincians are more confident in
recruiting higher risk patients when drugs are more mature), we do not
intend to adjust for this in the main analysis. In our view this would
be be conditioning on a mediator. If older trials have higher rates
because they are done at a more mature stage in the evaluation of a
drug, this does not reduce their applicability from the point of view of
the decision-maker.

All models will be fit using JAGS (a Bayesian software package), as this
will allow us to include the age distribution (from the mean and
standard deviation) as a covariate and not merely the mean age.

## Community comparison

First, we will compare the age-sex dsitribution of older people in the
community with older people in trials. We will do so graphcially, and
using summary statistics.

### Regression modelling in community

We will then estimate the age-sex specific incidence rate for first
hospitalisation or death in the SAIL platform, which contains community
data for people in Wales. We will identify first recorded prescription
of a RAAS drug for any participant with an previous diagnostic Read for
hypertension (using a code list described previously
(<https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-019-1427-1>)).
We will exclude participants who had been registered with a practice for
less than 6 months, to avoid including people who were already
established on a RAAS drug but who may have moved practice. We will
exclude those with myocardial infarction, stroke or TIA in the 6 months
prior to starting the RAAS drug (to attempt to exclude patients started
on these drugs for these indications rather than hypertension). We will
then model the event rate in Poisson regression models. In sensitivity
analyses we will repeat the modelling without excluding MI, stroke or
TIA.

### Estimating expected SAE ratios

Since, in adults, hospitalisation and death is the major driver of
serious adverse events, we will use these rates as well as the age-sex
distribution of each trial to calculate the expected SAE rate for each
trial. We will then calculate the ratio for the **expected** SAE rate
for older versus younger people trials and compare this to the
**estimated rate ratio** from the regression modelling using the trials
data.

We will obtain uncertainty intervals using parametric bootstrapping. We
will sample from a Poisson distribution taking as a parameter the
(exponentiated) age-sex specific linear predictors of the regression
equations derived from the community sample above, summing over the
age-sex distribution of the trials to obtain estimates for each trial.
We will then compare the rates for older versus non-older trials. We
will also obtain samples for the coefficient for the older-people status
variable from the trials regression models (which is straightforward as
these are fit using MCMC). We compare these two sets of samples to
obtain the ratio of ratios (estimated older/non-older versus expected
older/non-older). We hypothesise that expected ratio (based on community
rates) will be larger than the estimated ratio (based on the trials data
lone) - ie that older trials will have higher rates than non-older
trials, but not as much higher as expected based on community
hospitalisation and death rates.

For the puposes of this analysis we will treat the age-sex distribution
of the trials as fixed. From a separate project, we already know the
sex-specific age distribution for five trials from this dataset. We will
use this information to calculate the sex-specific mean and standard
deviation for age for the remainder of the trials by re-writing the
formula from Table 6.5a in the Cochrane Handbook
(<https://training.cochrane.org/handbook/current/chapter-06>). Then,
assuming that these distributions follow a truncated normal distribution
with known mean, sd, upper and lower limits, we will iterate over a
plausible possible range for the central tendency and dispersion in
order to calcualte these values. We will then use the cumulative
distribution function of the truncated normal distribution to estimate
the proportion of participants at each age and apply the community
linear predictors to these in order to estimate the trial-level expected
SAE rates.

    # https://www.r-bloggers.com/truncated-normal-distribution/
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
