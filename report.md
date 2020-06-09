Evaluating approaches to backcalculating cases counts by date of
infection from cases counts by date of report
================

**Authors:** EpiForecasts, CMMID Covid working group, Sebastian Funk

## Summary

  - The current methods being used to backcalculate case counts by date
    of infection assume either a constant shift (typically by the
    mean/median of the delay from onset to report and the mean/median of
    the incubation period) or assume independence between cases and
    sample both the delay from onset to report and the incubation
    period.

  - Parametric approaches are problematic in that they make assumptions
    about the data generating process that may not be correct.

  - Here we evaluate a non-parametric approach that is based on
    convoluting the known delays in comparison to our currently deployed
    independent sampling method.

  - We find that the non-parametric approach is better able to reproduce
    simulated data and returns plausible results when used on reported
    Covid-19 cases from 4 countries.

  - We also implement an in model approach to estimating the
    time-varying reproduction number but this work is at an early stage.

  - These findings are preliminary and more methodological work is
    required along with comparisons to other approaches.

## Introduction

## Dependencies

## Methods

### Model

  - Overview: non-parameteric backcalculation assuming a median shift
    with independent gaussian noise.

  - Orginal implementation here:
    <https://github.com/jhellewell14/rtconfirm/blob/master/R/backcalc.R>

  - This implementation uses median shifted reported cases (smoothed
    using a rolling average over the width of the median generation
    interval) as a prior and then fits independent gaussian noise on top
    of this. For future cases (i.e with no data to shift to into the
    last reported case count is used).

  - Daily reporting effects are included and constrained by weekly
    cases.

  - Uncertainty in reporting delays and incubation periods is enforced
    by upweighting the prior for these parameters.

  - We implement an in model time-varying reproduction number estimation
    step based on the renewal equation and the imputed infections.

<!-- end list -->

``` r
model
#> S4 class stanmodel 'nowcast' coded as follows:
#> functions {
#>   // convolve a pdf and case vector using matrix multiplication
#>   vector convolve(vector cases, vector pdf, int direction) {
#>     int t = num_elements(cases);
#>     matrix[t, t] delay_mat = rep_matrix(0, t, t);
#>     int max_pdf = num_elements(pdf);
#>     row_vector[max_pdf] row_pdf = to_row_vector(pdf);
#>     vector[t] convolved_cases;
#> 
#>     for (s in 1:t) {
#>       if (direction) {
#>         int max_length = min(s, max_pdf);
#>         delay_mat[s, (s - max_length + 1):s] = row_pdf[(max_pdf - max_length + 1):max_pdf];
#>       }else{
#>         int max_length = min(t - s, max_pdf - 1);
#>         delay_mat[s, s:(s + max_length)] = row_pdf[1:(max_length + 1)];
#>       }
#>     }
#>     
#>    convolved_cases = delay_mat * to_vector(cases);
#> 
#>    return(convolved_cases);
#>   }
#> 
#>   // apply backsampling and upscaling based on a pdf
#>   vector backsample(vector cases, vector pdf) {
#>     int t = num_elements(cases);
#>     vector[t] backsampled_cases;
#>     int max_upscale = min(t, num_elements(pdf));
#>     int pdf_length = num_elements(pdf);
#>     vector[pdf_length] cdf;
#>     
#>     backsampled_cases = convolve(cases, pdf, 0);
#>     
#>     // apply upscaling
#>     cdf = cumulative_sum(pdf);
#>     
#>     for (i in  1:max_upscale) {
#>       backsampled_cases[(t - i + 1)] = (backsampled_cases[(t - i + 1)] + 1) / cdf[i];
#>     }
#>     
#>     //bound last day to equal day before
#>     backsampled_cases[t] = backsampled_cases[t -1];
#>     
#>     return(backsampled_cases);
#>   }
#>   
#>   // discretised lognormal pmf
#>   real discretised_lognormal_pmf(int y, real mu, real sigma) {
#>     return((lognormal_cdf(y, mu, sigma) - lognormal_cdf(y - 1, mu, sigma)));
#>   }
#>   
#>   // discretised gamma pmf
#>   real discretised_gamma_pmf(int y, real mu, real sigma) {
#>     // calculate alpha and beta for gamma distribution
#>     real alpha = (mu / sigma)^2;
#>     real beta = (sigma^2) / mu;
#>     return((gamma_cdf(y, alpha, beta) - gamma_cdf(y - 1, alpha, beta)));
#>   }
#> }
#> 
#> 
#> data {
#>   int t;                             // number of time steps
#>   int day_of_week[t];                // day of the week indicator (1 - 7)
#>   int <lower = 0> cases[t];          // observed cases
#>   vector<lower = 0>[t] shifted_cases;// median shifted smoothed cases
#>   real inc_mean_sd;                  // prior sd of mean incubation period
#>   real inc_mean_mean;                // prior mean of mean incubation period
#>   real inc_sd_mean;                  // prior sd of sd of incubation period
#>   real inc_sd_sd;                    // prior sd of sd of incubation period
#>   int max_inc;                       // maximum incubation period
#>   real rep_mean_mean;                // prior mean of mean reporting delay
#>   real rep_mean_sd;                  // prior sd of mean reporting delay
#>   real rep_sd_mean;                  // prior mean of sd of reporting delay
#>   real rep_sd_sd;                    // prior sd of sd of reporting delay
#>   int max_rep;                       // maximum report delay
#>   real <lower = 0> r_mean;           // prior mean of reproduction number
#>   real <lower = 0> r_sd;             // prior standard deviation of reproduction number
#>   real gt_mean_sd;                   // prior sd of mean generation time
#>   real gt_mean_mean;                 // prior mean of mean generation time
#>   real gt_sd_mean;                   // prior sd of sd of generation time
#>   real gt_sd_sd;                     // prior sd of sd of generation time
#>   int max_gt;                        // maximum generation time
#>   int model_type;                    // type of model: 1 = poisson otherwise negative binomial
#>   int estimate_r;                    // should the reproduction no be estimated (1 = yes)
#> }
#> 
#> transformed data{
#>   real r_alpha;                      // alpha parameter of the R gamma prior
#>   real r_beta;                       // beta parameter of the R gamma prior
#>   real initial_noise_mean;           // initial noise estimate
#>   real initial_R_mean;               // initial R estimate
#> 
#>   // calculate alpha and beta for gamma distribution
#>   r_alpha = (r_mean / r_sd)^2;
#>   r_beta = (r_sd^2) / r_mean;
#>   
#>   // precalculate parameters for noise
#>   initial_noise_mean = log(1);
#>   
#>   // precalculate parameters for R
#>   initial_R_mean = log(1);
#> }
#> parameters{
#>   real<lower = 0> initial_noise;                      // initial noise on the mean shifted observed cases
#>   vector<lower = 0>[t -1] noise_diff;                 // temporal scaling on noise
#>   simplex[7] day_of_week_eff_raw;                     // day of week reporting effect + control parameters
#>   real <lower = 0> inc_mean;                          // mean of incubation period
#>   real <lower = 0> inc_sd;                            // sd of incubation period
#>   real <lower = 0> rep_mean;                          // mean of reporting delay
#>   real <lower = 0> rep_sd;                            // sd of incubation period
#>   real<lower = 0> rep_phi;                            // overdispersion of the reporting process
#>   vector<lower = 0>[estimate_r] initial_R;            // initial reproduction number
#>   vector<lower = 0>[estimate_r > 0 ? t -1 : 0] R_diff;// scaling of reproduction number over time
#>   real<lower = 0> gt_mean[estimate_r];                // mean of generation time
#>   real <lower = 0> gt_sd[estimate_r];                 // sd of generation time
#>   real<lower = 0> inf_phi[estimate_r];                // overdispersion of the infection process
#> }
#> 
#> transformed parameters {
#>   // stored transformed parameters
#>   vector<lower = 0>[t] noise;                             // noise on the mean shifted observed cases
#>   vector<lower = 0>[t] infections;                        // infections over time
#>   vector<lower = 0>[t] reports;                           // reports over time
#>   vector[7] day_of_week_eff;                              // day of the week effect
#>   vector<lower = 0>[estimate_r > 0 ? t : 0] R;            // reproduction number over time
#>   vector[estimate_r > 0 ? t : 0] branch_reports;          // reports generated by the branching process
#>  {
#>   // temporary transformed parameters
#>   vector[max_rep] rev_delay;                              // reversed report delay pdf
#>   vector[max_inc] rev_incubation;                         // reversed incubation period pdf
#>   vector[t] onsets;                                       // onsets over time
#>   vector[estimate_r > 0 ? max_gt : 0] rev_generation_time;// reversed generation time pdf
#>   vector[estimate_r > 0 ? t : 0] infectiousness;          // infections over time
#>   vector[estimate_r > 0 ? t : 0] branch_onsets;           // onsets generated by the branching process
#>   vector[estimate_r > 0 ? t : 0] branch_infections;       // infections generated by the branching process
#>   
#>   // reverse the distributions to allow vectorised access
#>   for (j in 1:max_inc) {
#>     rev_incubation[j] =
#>         discretised_lognormal_pmf(max_inc - j + 1, inc_mean, inc_sd);
#>   }
#>   
#>   for (j in 1:max_rep) {
#>     rev_delay[j] =
#>         discretised_lognormal_pmf(max_rep - j + 1, rep_mean, rep_sd);
#>   }
#>     
#>   // define day of the week effect
#>   day_of_week_eff = 7 * day_of_week_eff_raw;
#> 
#>   // walk over noise terms
#>   noise[1] = initial_noise;
#>   for (s in 2:t) {
#>     noise[s] = noise[s - 1] * noise_diff[s - 1];
#>   }
#>   
#>   // generate infections from backcalculated and non-parameteric noise
#>   for (s in 1:t) {
#>     infections[s] = shifted_cases[s] * noise[s];
#>   }
#>   
#>   // onsets from infections
#>   onsets = convolve(infections, rev_incubation, 1);
#> 
#>   // reports from onsets
#>   reports = convolve(onsets, rev_delay, 1);
#> 
#>   for (s in 1:t) {
#>     // add reporting effects (adjust for simplex scale)
#>     reports[s] *= day_of_week_eff[day_of_week[s]];
#>   }
#>     
#>   ////////////////////////////////////////////////////// 
#>   // estimate reproduction no from a branching process
#>   //////////////////////////////////////////////////////
#>   if (estimate_r) {
#>     // calculate pdf of generation time from distribution
#>     for (j in 1:(max_gt)) {
#>        rev_generation_time[j] =
#>            discretised_gamma_pmf(max_gt - j + 1, gt_mean[estimate_r], gt_sd[estimate_r]);
#>      }
#>      
#>      // infectiousness from infections
#>      infectiousness = convolve(infections, rev_generation_time, 1);
#>   
#>      // reproduction number is a scaling of previous timestep
#>      R[1] = initial_R[estimate_r];
#>      for (s in 2:t) {
#>        R[s] = R[s - 1] * R_diff[s - 1];
#>      }
#>      
#>      // Estimate infections using branching process
#>      branch_infections = R .* infectiousness;
#>      
#>      // onsets from infections
#>      branch_onsets = convolve(branch_infections, rev_incubation, 1);
#> 
#>      // reports from onsets
#>      branch_reports = convolve(branch_onsets, rev_delay, 1);
#> 
#>      for (s in 1:t) {
#>       // add reporting effects (adjust for simplex scale)
#>       branch_reports[s] *= day_of_week_eff[day_of_week[s]];
#>     }
#>   }
#>  }
#> }
#> 
#> model {
#>   // reporting overdispersion
#>   rep_phi ~ exponential(1);
#> 
#>   // noise on median shift
#>   initial_noise ~ lognormal(initial_noise_mean, 0.1);
#>   noise_diff ~ lognormal(initial_noise_mean, 0.1);
#> 
#>   // daily cases given reports
#>   if (model_type == 1) {
#>     target += poisson_lpmf(cases | reports);
#>   }else{
#>     target += neg_binomial_2_lpmf(cases | reports, rep_phi);
#>   }
#> 
#>   // penalised priors for incubation period, and report delay
#>   target += normal_lpdf(inc_mean | inc_mean_mean, inc_mean_sd) * t;
#>   target += normal_lpdf(inc_sd | inc_sd_mean, inc_sd_sd) * t;
#>   target += normal_lpdf(rep_mean | rep_mean_mean, rep_mean_sd) * t;
#>   target += normal_lpdf(rep_sd | rep_sd_mean, rep_sd_sd) * t;
#>   
#>   ////////////////////////////////////////////////////// 
#>   // estimate reproduction no from a branching process
#>   //////////////////////////////////////////////////////
#>   if (estimate_r) {
#>     
#>     // infection overdispersion
#>     inf_phi ~ exponential(1);
#> 
#>     // initial prior on R
#>     initial_R ~ gamma(r_alpha, r_beta);
#>     R_diff ~ lognormal(initial_R_mean, 1);
#>      
#>     // penalised_prior on generation interval
#>     target += normal_lpdf(gt_mean | gt_mean_mean, gt_mean_sd) * t;
#>     target += normal_lpdf(gt_sd | gt_sd_mean, gt_sd_sd) * t;
#>   
#>     // Likelihood of Rt given infections
#>     if (model_type == 1) {
#>       target += poisson_lpmf(cases | branch_reports);
#>     }else{
#>       target += neg_binomial_2_lpmf(cases | branch_reports, inf_phi[estimate_r]);
#>     }
#>   }
#> }
#>   
#> generated quantities {
#>   int imputed_infections[t];
#>   
#>   // simulated infections - assume poisson (with negative binomial reporting)
#>   imputed_infections = poisson_rng(infections);
#> }
```

## Analysis

### Simulate data

  - Define a realistic basic reproduction number estimate that starts at
    2, decreases linearly to 0.5, remains constant, increases linearly
    to 1.2, decreases sharply to 1, again remains constant and finally
    increases linearly. Noise is added (2% standard deviation) for all
    values.

<!-- end list -->

``` r
## Define an initial rt vector 
rts <- c(rep(2, 20), (2 - 1:15 * 0.1), rep(0.5, 10), (0.5 + 1:7 * 0.1), rep(1.2, 20),
         rep(1, 10), (1 + 1:10 * 0.05))
## Add noise
##rts <- rts * rnorm(length(rts), 1, 0.02)

rts
#>  [1] 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00 2.00
#> [16] 2.00 2.00 2.00 2.00 2.00 1.90 1.80 1.70 1.60 1.50 1.40 1.30 1.20 1.10 1.00
#> [31] 0.90 0.80 0.70 0.60 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50
#> [46] 0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20
#> [61] 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.20 1.00 1.00 1.00
#> [76] 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40
#> [91] 1.45 1.50
```

  - In order to simulate cases by date of infection from a reproduction
    number trace an estimate of the generation interval is required. We
    use the current default estimate for Covid-19 from `EpiNow`.

<!-- end list -->

``` r
## Use the mean default generation interval for covid
generation_time_pdf <- rowMeans(EpiNow::covid_generation_times)

generation_time_pdf
#>  [1] 0.000000e+00 2.063055e-01 1.887694e-01 1.602335e-01 1.240479e-01
#>  [6] 9.071765e-02 6.487144e-02 4.603137e-02 3.269885e-02 2.340093e-02
#> [11] 1.673225e-02 1.204701e-02 8.778529e-03 6.435373e-03 4.739174e-03
#> [16] 3.587893e-03 2.700888e-03 2.063697e-03 1.558217e-03 1.195100e-03
#> [21] 9.091163e-04 6.717316e-04 5.056617e-04 3.622574e-04 2.497779e-04
#> [26] 1.576857e-04 1.031422e-04 6.222018e-05 3.392786e-05 1.669864e-05
#> [31] 7.260972e-06 3.493567e-06 1.600534e-06 6.507847e-07 1.398601e-07
```

``` r
incubation_period <- list(mean = EpiNow::covid_incubation_period[1, ]$mean,
                           mean_sd = EpiNow::covid_incubation_period[1, ]$mean_sd,
                           sd = EpiNow::covid_incubation_period[1, ]$sd,
                           sd_sd = EpiNow::covid_incubation_period[1, ]$sd_sd,
                           max = 30)
                    
reporting_delay <- list(mean = log(5),
                         mean_sd = log(1.1),
                         sd = log(2),
                         sd_sd = log(1.2),
                         max = 30)

generation_time <- list(mean = EpiNow::covid_generation_times_summary[1, ]$mean,
                        mean_sd = EpiNow::covid_generation_times_summary[1, ]$mean_sd,
                        sd = EpiNow::covid_generation_times_summary[1, ]$sd,
                        sd_sd = EpiNow::covid_generation_times_summary[1, ]$sd_sd,
                        max = 30)

## Sample a report delay as a lognormal - take 10 samples
delay_defs <- EpiNow::lognorm_dist_def(mean = reporting_delay$mean, mean_sd = reporting_delay$mean_sd,
                                       sd = reporting_delay$sd, sd_sd = reporting_delay$sd_sd,
                                       max_value = reporting_delay$max, samples = 1000)


## Sample a incubation period (again using the default for covid) - take 10 samples
incubation_defs <- EpiNow::lognorm_dist_def(mean = incubation_period$mean,
                                           mean_sd = incubation_period$mean_sd,
                                           sd = incubation_period$sd,
                                           sd_sd = incubation_period$sd_sd,
                                           max_value = 30, samples = 1000)

## Simulate cases with a decrease in reporting at weekends and an incease on Monday
## using a single sample of both distributions                                    
simulated_cases <- EpiNow::simulate_cases(rts, initial_cases = 10 , initial_date = as.Date("2020-03-01"),
                                          generation_interval = generation_time_pdf,
                                          delay_def = delay_defs[10, ],
                                          incubation_def = incubation_defs[10, ],
                                          reporting_effect = c(1.6, 1, 1, 1, 0.8, 0.4, 1))
simulated_cases
#>            date cases reference
#>   1: 2020-03-02     7 infection
#>   2: 2020-03-03     6 infection
#>   3: 2020-03-04     6 infection
#>   4: 2020-03-05    10 infection
#>   5: 2020-03-06    15 infection
#>  ---                           
#> 267: 2020-05-27   192    report
#> 268: 2020-05-28   146    report
#> 269: 2020-05-29   144    report
#> 270: 2020-05-30    77    report
#> 271: 2020-05-31   214    report
```

### Compare approaches on simulated data

``` r
## Extract simulated infections
simulated_reports <- simulated_cases[reference == "report"][, confirm := cases][, cases := NULL]

## Rt prior
rt_prior <- list(mean = 2.6, sd = 2)

## Reconstruction via backwards sampling
sampling_cases <- nowcast_pipeline(reported_cases = simulated_reports[, import_status := "local"], 
                                   target_date = max(simulated_reports$date),
                                   delay_defs = delay_defs, 
                                   incubation_defs = incubation_defs,
                                   nowcast_lag = 0, approx_delay = TRUE)
#> Warning in stats::rnbinom(n, x + 1, prob): NAs produced

## Non-parameteric reconstruction
non_parametric_cases <- nowcast(simulated_reports,
                                family = "negbin",
                                incubation_period = incubation_period,
                                reporting_delay = reporting_delay,
                                generation_time = generation_time, 
                                estimate_rt = FALSE,
                                rt_prior = rt_prior,
                                cores = 4, chains = 4,
                                samples = 1000, warmup = 1000,
                                return_all = TRUE, model = model,
                                verbose = TRUE)
#> Running for 2000 samples and 99 time steps
```

### Compare approaches on reported Covid-19 cases in Austria, the United Kingdom, United States of America and Russia

  - Get
data

<!-- end list -->

``` r
reported_cases <- NCoVUtils::get_ecdc_cases(countries = c("Austria", "United_Kingdom",
                                                          "United_States_of_America", "Russia"))
reported_cases <- NCoVUtils::format_ecdc_data(reported_cases)
reported_cases <- data.table::as.data.table(reported_cases)[, confirm := cases][, cases := NULL]
reported_cases <- reported_cases[date >= "2020-02-01"]
```

  - Run backcalculation on each country in
turn

<!-- end list -->

``` r
countries <- c("Austria", "United Kingdom", "United States of America", "Russia")

results <- lapply(countries,
                  function(country) {
        message("Nowcasting using sampling for: ", country)                
        cases <- data.table::copy(reported_cases)[region %in% country]  
                                        
        ## Reconstruction via backwards sampling
        sampling_cases <- nowcast_pipeline(reported_cases = cases[, import_status := "local"], 
                                           target_date = max(cases$date),
                                           delay_defs = delay_defs,
                                           incubation_defs = incubation_defs,
                                           nowcast_lag = 0, approx_delay = TRUE)
        
        message("Non-parametric nowcasting for: ", country)
        ## Non-parametric reconstruction
        non_parametric_cases <- nowcast(cases,
                                        family = "negbin",
                                        incubation_period = incubation_period,
                                        reporting_delay = reporting_delay,
                                        generation_time = generation_time, 
                                        estimate_rt = FALSE,
                                        samples = 2000, warmup = 1000,
                                        cores = 4, chains = 4,
                                        return_all = TRUE, model = model)
        
        return(list(sampling_cases, non_parametric_cases))
                                       })
#> Warning in stats::rnbinom(n, x + 1, prob): NAs produced

#> Warning in stats::rnbinom(n, x + 1, prob): NAs produced
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
#> Warning in stats::rnbinom(n, x + 1, prob): NAs produced

#> Warning in stats::rnbinom(n, x + 1, prob): NAs produced
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess

names(results) <- countries
```

## Results

### Simulated data

  - Reporting effects.

<!-- end list -->

``` r
non_parametric_cases$day_of_week[, as.list(summary(value)), by = "wday"]
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#> 1:    Monday 1.3321460 1.5882585 1.6496200 1.6516356 1.7106489 1.9419309
#> 2:   Tuesday 0.7691495 0.9483523 0.9913809 0.9938828 1.0343125 1.2626165
#> 3: Wednesday 0.8299799 0.9997393 1.0424612 1.0439416 1.0860256 1.3065855
#> 4:  Thursday 0.7630133 0.9711470 1.0132362 1.0146981 1.0553339 1.2301869
#> 5:    Friday 0.6734807 0.8026501 0.8400650 0.8408368 0.8777410 1.0472130
#> 6:  Saturday 0.3083335 0.3791148 0.3992795 0.4001306 0.4199423 0.5164145
#> 7:    Sunday 0.8195839 1.0119234 1.0513874 1.0548746 1.0964231 1.3402281
```

  - Recover reporting delays

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$rep_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$rep_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean          sd
#> 1:      mean 1.6085744 0.009525181
#> 2:        sd 0.6916494 0.017710901
```

  - Recover incubation period

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$inc_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$inc_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean          sd
#> 1:      mean 1.6207186 0.006204358
#> 2:        sd 0.4177093 0.006923570
```

  - Prepare data for
plotting

<!-- end list -->

``` r
simulated_cases <- simulated_cases[reference %in% c("infection", "report")][, median := cases][,
                                   type := ifelse(reference == "infection", 
                                                  "Simulated infections", 
                                                  "Simulated reported cases")][,
                                   `:=`(cases = NULL, reference = NULL)]

summarise_nowcasting_approaches <- function(sampling_cases, non_parametric_cases) {
  sampling_cases <- sampling_cases[type %in% "infection_upscaled"][,
               .(median = median(cases), bottom = quantile(cases, 0.025),
                 lower = quantile(cases, 0.25), upper = quantile(cases, 0.75),
                 top = quantile(cases, 0.975)), by = c("date", "type")][,
                 type := "Sampled"]

non_parametric_infections <- non_parametric_cases$infections[,
               .(median = median(value), bottom = quantile(value, 0.025),
                 lower = quantile(value, 0.25), upper = quantile(value, 0.75),
                 top = quantile(value, 0.975)), by = c("date")][,
                 type := "Non-parametric"]

non_parametric_prior_infections <- non_parametric_cases$prior_infections[sample == 1,
               .(median = value), by = c("date")][,
                 type := "Infections prior"]
out <- data.table::rbindlist(list(sampling_cases, non_parametric_infections,
                                  non_parametric_prior_infections), fill = TRUE)
return(out)
}

summarised_nowcasting_approaches <- summarise_nowcasting_approaches(sampling_cases, non_parametric_cases)

simulated_cases <- data.table::rbindlist(list(simulated_cases,
                                              summarised_nowcasting_approaches),
                                         fill = TRUE)
```

  - Plot date. *Note: Here we have cut-off cases prior to the start of
    March. This truncates the long tale observed in the sampling
    approach.*

<!-- end list -->

``` r
plot_data <- simulated_cases[date >= as.Date("2020-03-01")]

plot_obs <- ggplot2::ggplot(plot_data[type %in% "Simulated reported cases"], 
                           ggplot2::aes(x = date, col = type, fill = type)) +
  ggplot2::geom_col(ggplot2::aes(y = median), fill = "grey", col = "white",
                    show.legend = FALSE) +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(y = "Reported cases", x = "")


plot_sims_data <- plot_data[!type %in% "Simulated reported cases"][,
                             type := factor(type, levels = c("Simulated infections",
                                                             "Sampled",
                                                             "Infections prior",
                                                             "Non-parametric"))]
plot <- ggplot2::ggplot(plot_sims_data, ggplot2::aes(x = date, col = type, fill = type)) +
 ggplot2::geom_linerange(data = plot_sims_data[!type %in% "Simulated infections"],
                          ggplot2::aes(ymin = bottom, ymax = top), 
                         alpha = 0.2, size = 1.5) +
 ggplot2::geom_linerange(data = plot_sims_data[!type %in% "Simulated infections"],
                         ggplot2::aes(ymin = lower, ymax = upper), 
                         alpha = 0.4, size = 1.5) +
 ggplot2::geom_line(data = plot_sims_data[type %in% c("Simulated infections", "Infections prior")],
                    ggplot2::aes(y = median), size = 1.1, alpha = 0.6) +
 ggplot2::geom_point(data = plot_sims_data[type %in%  c("Simulated infections", "Infections prior")],
                     ggplot2::aes(y = median), size = 1.1, alpha = 0.8, show.legend = FALSE) +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(y = "Cases by date of infection", x = "Date", col = "Type")

plot_obs + plot + patchwork::plot_layout(ncol = 1)
```

<img src="figures/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

### Reported Covid-19 cases in the United Kingdom, United States of America and South Korea

  - Explore reporting effects by
country

<!-- end list -->

``` r
purrr::map(results, ~ .[[2]]$day_of_week[, as.list(summary(value)), by = "wday"])
#> $Austria
#>         wday      Min.   1st Qu.    Median      Mean  3rd Qu.     Max.
#> 1:    Monday 0.4553486 0.6597074 0.7133279 0.7191980 0.773323 1.181415
#> 2:   Tuesday 0.6435543 0.8928487 0.9620136 0.9676083 1.037469 1.487350
#> 3: Wednesday 0.6217315 0.9533216 1.0260229 1.0319547 1.102476 1.551721
#> 4:  Thursday 0.6351770 0.8889191 0.9575475 0.9639907 1.033012 1.434612
#> 5:    Friday 0.7145081 1.0078820 1.0844727 1.0910157 1.166702 1.736259
#> 6:  Saturday 0.8535405 1.1456865 1.2271969 1.2338026 1.316128 1.818715
#> 7:    Sunday 0.6418429 0.9157637 0.9863031 0.9924299 1.060449 1.463494
#> 
#> $`United Kingdom`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.6234883 0.9875295 1.0647041 1.0714275 1.1473135 1.795143
#> 2:   Tuesday 0.5369750 0.7564222 0.8184997 0.8237987 0.8834020 1.288060
#> 3: Wednesday 0.6056202 0.8462580 0.9121482 0.9172414 0.9814343 1.414607
#> 4:  Thursday 0.6096313 0.8767529 0.9490498 0.9557929 1.0271495 1.476008
#> 5:    Friday 0.6911887 0.9549942 1.0297273 1.0357731 1.1112875 1.563819
#> 6:  Saturday 0.7023413 0.9896022 1.0652523 1.0720039 1.1476437 1.560209
#> 7:    Sunday 0.7346365 1.0391061 1.1184163 1.1239625 1.2031509 1.705964
#> 
#> $`United States of America`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7090682 0.8663433 0.9029952 0.9045882 0.9408033 1.199528
#> 2:   Tuesday 0.7588929 0.9038461 0.9415399 0.9439228 0.9813710 1.217412
#> 3: Wednesday 0.7501136 0.8849629 0.9213808 0.9232633 0.9586499 1.190471
#> 4:  Thursday 0.7354623 0.9123161 0.9508406 0.9533327 0.9924186 1.220941
#> 5:    Friday 0.8049083 1.0224692 1.0645537 1.0661633 1.1085772 1.397414
#> 6:  Saturday 0.9236516 1.1016174 1.1457626 1.1489255 1.1921687 1.414941
#> 7:    Sunday 0.8111929 1.0179732 1.0585695 1.0598043 1.0998705 1.362117
#> 
#> $Russia
#>         wday      Min.   1st Qu.    Median      Mean  3rd Qu.     Max.
#> 1:    Monday 0.6945279 0.9380559 0.9992149 1.0034149 1.063491 1.510367
#> 2:   Tuesday 0.6705700 0.8975525 0.9543772 0.9586771 1.015618 1.342484
#> 3: Wednesday 0.6615815 0.8873466 0.9431653 0.9475883 1.003029 1.358241
#> 4:  Thursday 0.6502685 0.9620644 1.0264623 1.0310391 1.093293 1.479560
#> 5:    Friday 0.7254960 0.9480915 1.0096342 1.0140335 1.074367 1.412991
#> 6:  Saturday 0.7355654 1.0195432 1.0863036 1.0929837 1.159355 1.594220
#> 7:    Sunday 0.6639357 0.8895024 0.9481209 0.9522635 1.010432 1.389664
```

  - Prepare data for
plotting

<!-- end list -->

``` r
summarised_cases <- data.table::copy(reported_cases)[, median := confirm][, type := "Reported Cases"]


summarised_results <- purrr::map2(results, names(results),
                                  ~ summarise_nowcasting_approaches(.x[[1]], .x[[2]])[,region := .y])

summarised_results <- data.table::rbindlist(summarised_results)
all_country_data <- data.table::rbindlist(list(summarised_cases, summarised_results), fill = TRUE)
all_country_data <- all_country_data[date >= "2020-03-01"]
```

  - Plot
data

<!-- end list -->

``` r
plot <- ggplot2::ggplot(all_country_data, ggplot2::aes(x = date, col = type, fill = type)) +
  ggplot2::geom_col(data = all_country_data[type %in% "Reported Cases"],
                    ggplot2::aes(y = median), fill = "grey", col = "white") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = bottom, ymax = top), 
                         alpha = 0.2, size = 1) +
 ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), 
                         alpha = 0.4, size = 1) +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(y = "Cases", x = "Date", col = "Type") + 
  ggplot2::facet_wrap(~region, scales = "free_y")

plot
```

<img src="figures/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

## Discussion

## References
