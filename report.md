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

  - Overview: non-parameteric backcalculation assuming a mean shift with
    gaussian noise.

  - This implementation uses mean shifted reported cases (smoothed using
    a 7 da, backwards looking, rolling average) as a prior combined with
    a guassian process. For future cases (i.e with no data to shift to
    into the last reported case count is used).

  - Fixed daily reporting effects are adjusted for.

  - Uncertainty in reporting delays and incubation periods is enforced
    by upweighting the prior for these parameters.

  - We implement an in model time-varying reproduction number estimation
    step based on the renewal equation, the imputed infections and a
    gaussian process.

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
#>     return(normal_cdf((log(y + 1) - mu) / sigma, 0.0, 1.0) - normal_cdf((log(y) - mu) / sigma, 0.0, 1.0));
#>   }
#>   
#>   // discretised gamma pmf
#>   real discretised_gamma_pmf(int y, real mu, real sigma) {
#>     // calculate alpha and beta for gamma distribution
#>     real alpha = ((mu)/ sigma)^2;
#>     real beta = (mu) / (sigma^2);
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
#>   int model_type;                    // type of model: 0 = poisson otherwise negative binomial
#>   int estimate_r;                    // should the reproduction no be estimated (1 = yes)
#> }
#> 
#> transformed data{
#>   real r_alpha;                      // alpha parameter of the R gamma prior
#>   real r_beta;                       // beta parameter of the R gamma prior
#>   real delta;                        // modifier to make GP + definite
#>   real time[t];                      // time vector
#>   
#>   // calculate alpha and beta for gamma distribution
#>   r_alpha = (r_mean / r_sd)^2;
#>   r_beta = r_mean / (r_sd^2);
#>   
#>   // Assign + definite term
#>    delta = 1e-9;
#>   
#>   // Make time vector
#>    for (s in 1:t) {
#>      time[s] = s;
#>    }
#> }
#> parameters{
#>   simplex[7] day_of_week_eff_raw;                     // day of week reporting effect + control parameters
#>   real <lower = 0> inc_mean;                          // mean of incubation period
#>   real <lower = 0> inc_sd;                            // sd of incubation period
#>   real <lower = 0> rep_mean;                          // mean of reporting delay
#>   real <lower = 0> rep_sd;                            // sd of incubation period
#>   real<lower = 0> rep_phi[model_type];                // overdispersion of the reporting process
#>   real<lower=0> rho;                                  // length scale of noise GP
#>   real<lower=0> alpha;                                // scale of of noise GP
#>   vector[t] eta;                                      // unconstrained noise
#>   vector<lower = 0>[estimate_r] initial_R;            // baseline reproduction number estimate
#>   real<lower = 1> gt_mean[estimate_r];                // mean of generation time
#>   real <lower = 0> gt_sd[estimate_r];                 // sd of generation time
#>   real<lower=0> inf_phi[model_type*estimate_r];       // overdispersion of the infection process
#>   real<lower=0> R_rho[estimate_r];                    // length scale of R GP
#>   real<lower=0> R_alpha[estimate_r];                  // scale of R GP
#>   vector[estimate_r > 0 ? t : 0] R_eta;               // unconstrained R noise
#> }
#> 
#> transformed parameters {
#>   // stored transformed parameters
#>   vector<lower = 0>[t] noise;                             // noise on the mean shifted observed cases
#>   vector<lower = 0>[t] infections;                        // infections over time
#>   vector<lower = 0>[t] reports;                           // reports over time
#>   vector[7] day_of_week_eff;                              // day of the week effect
#>   vector[estimate_r > 0 ? t : 0] branch_reports;          // reports generated by the branching process
#>   vector[estimate_r > 0 ? t : 0] R;                       // reproduction number over time
#>   
#>  {
#>   // temporary transformed parameters
#>   vector[max_rep] rev_delay;                              // reversed report delay pdf
#>   vector[max_inc] rev_incubation;                         // reversed incubation period pdf
#>   vector[t] onsets;                                       // onsets over time
#>   vector[estimate_r > 0 ? max_gt : 0] rev_generation_time;// reversed generation time pdf
#>   vector[estimate_r > 0 ? t : 0] infectiousness;          // infections over time
#>   vector[estimate_r > 0 ? t : 0] branch_onsets;           // onsets generated by the branching process
#>   vector[estimate_r > 0 ? t : 0] branch_infections;       // infections generated by the branching process
#>   matrix[t, t] K;                                         // covariance matrix 
#>   matrix[t, t] L_K;                                       // cholesky docomposed covariance matrix
#>   vector[estimate_r > 0 ? t : 0] R_noise;                 // noise on R centred on 1
#>   
#>   // reverse the distributions to allow vectorised access
#>   for (j in 1:max_inc) {
#>     rev_incubation[j] =
#>         discretised_lognormal_pmf(max_inc - j + 1, inc_mean, inc_sd);
#>   }
#>   
#>   for (j in 1:(max_rep)) {
#>     rev_delay[j] =
#>         discretised_lognormal_pmf(max_rep - j + 1, rep_mean, rep_sd);
#>   }
#>     
#>   // define day of the week effect
#>   day_of_week_eff = 7 * day_of_week_eff_raw;
#> 
#>   // GP in noise
#>   K = cov_exp_quad(time, alpha, rho);
#>   // diagonal elements with offset to make + definite
#>   for (n in 1:t) {
#>       K[n, n] = K[n, n] + delta;
#>   }
#>   
#>   L_K = cholesky_decompose(K);
#>   noise = exp(L_K * eta);
#>   
#>   for (s in 1:t) {
#>     if(noise[s] == 0) {
#>       noise[s] = 0.0001;
#>     }
#>   }
#> 
#>   // generate infections from prior infections and non-parameteric noise
#>   infections = shifted_cases .* noise;
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
#>     for (j in 1:(max_gt - 1)) {
#>        rev_generation_time[j] =
#>            discretised_gamma_pmf(max_gt - j, gt_mean[estimate_r], gt_sd[estimate_r]);
#>      }
#>      
#>      // set same day to be 0
#>      rev_generation_time[max_gt] = 0;
#>      
#>      infectiousness = convolve(infections, rev_generation_time, 1);
#> 
#>      // Construct R over time
#>       K = cov_exp_quad(time, R_alpha[estimate_r], R_rho[estimate_r]);
#>      // diagonal elements with offset to make + definite
#>       for (n in 1:t) {
#>         K[n, n] = K[n, n] + delta;
#>         }
#>       L_K = cholesky_decompose(K);
#>       R_noise = exp(L_K * R_eta);
#>   
#>      // Estimate infections using branching process
#>       for (s in 1:t) {
#>          R[s] = initial_R[estimate_r] * R_noise[s];
#>          branch_infections[s] = R[s] * infectiousness[s];
#>          
#>          // Make sure all dates have a non-zero value
#>          if (branch_infections[s] == 0){
#>             branch_infections[s] = 0.0001; 
#>          }
#>        }
#>      // onsets from infections
#>      branch_onsets = convolve(branch_infections, rev_incubation, 1);
#> 
#>      // reports from onsets
#>      branch_reports = convolve(branch_onsets, rev_delay, 1);
#> 
#>      for (s in 1:t) {
#>       // add reporting effects (adjust for simplex scale)
#>       branch_reports[s] *= (day_of_week_eff[day_of_week[s]]);
#>     }
#>   }
#>  }
#> }
#> 
#> model {
#>   
#>   // priors for noise GP
#>   rho ~ lognormal(1.098612, 0.5); //log(3)
#>   alpha ~ std_normal();
#>   eta ~ std_normal();
#> 
#>   // reporting overdispersion
#>   if (model_type) {
#>     rep_phi[model_type] ~ exponential(1);
#>   }
#>   
#>   // daily cases given reports
#>   if (model_type) {
#>     target += neg_binomial_2_lpmf(cases | reports, rep_phi[model_type]);
#>   }else{
#>     target += poisson_lpmf(cases | reports);
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
#>     // infection overdispersion
#>     if(model_type) {
#>       inf_phi ~ exponential(1);
#>     }
#> 
#>     // prior on R
#>     initial_R[estimate_r] ~ gamma(r_alpha, r_beta);
#>     
#>    // priors for R gp
#>    R_rho ~ lognormal(1.098612, 0.5); //log(3)
#>    R_alpha ~ std_normal();
#>    R_eta ~ std_normal();
#>     
#>     // penalised_prior on generation interval
#>     target += normal_lpdf(gt_mean | gt_mean_mean, gt_mean_sd) * t;
#>     target += normal_lpdf(gt_sd | gt_sd_mean, gt_sd_sd) * t;
#>     
#>     // Likelihood of Rt given infections
#>     if (model_type) {
#>       target += neg_binomial_2_lpmf(cases | branch_reports, inf_phi[model_type*estimate_r]);
#>     }else{
#>       target += poisson_lpmf(cases | branch_reports);
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
rts <- c(rep(2.5, 20), rep(2, 20), (2 - 1:15 * 0.1),
         rep(0.5, 10), (0.5 + 1:7 * 0.1), rep(1.2, 20),
         rep(1, 10), (1 + 1:10 * 0.05))
## Add noise
rts <- rts * rnorm(length(rts), 1, 0.05)

rts
#>   [1] 2.5623140 2.4980184 2.6991919 2.5367327 2.5573963 2.5170286 2.3292668
#>   [8] 2.4248246 2.3889409 2.5049281 2.3535154 2.5900088 2.4260866 2.4700671
#>  [15] 2.6694465 2.6434026 2.6261427 2.7828025 2.5531542 2.8023181 1.9909193
#>  [22] 2.0857080 1.8729938 2.0659323 1.8922478 2.0690301 1.9532624 1.9774066
#>  [29] 1.9665504 2.0914207 2.0584566 1.8094693 2.0640880 2.1685164 2.2914501
#>  [36] 2.0138209 2.0441051 1.8115966 2.0598780 2.0477991 1.8842648 1.8587768
#>  [43] 1.6208323 1.4799803 1.6198850 1.3957192 1.4105012 1.3260171 1.0910410
#>  [50] 1.1129816 0.9626194 0.8206260 0.7229856 0.5606589 0.4943175 0.4915417
#>  [57] 0.4713031 0.4812145 0.4829283 0.4862380 0.5058287 0.4970090 0.5148365
#>  [64] 0.4537586 0.4815822 0.6069229 0.6626104 0.8139226 0.9109149 0.9881556
#>  [71] 1.0760357 1.1699750 1.2105936 1.1676167 1.1505080 1.1426990 1.1710122
#>  [78] 1.2797216 1.1863771 1.2121202 1.2548112 1.1991243 1.1194449 1.2834417
#>  [85] 1.1181440 1.1757061 1.3343713 1.2809297 1.1067026 1.2093056 1.2906514
#>  [92] 1.2484025 0.9308492 0.9544376 1.0771419 0.9978504 1.0633707 0.9845634
#>  [99] 1.0444121 0.9344635 1.0370130 1.0232480 1.1055137 1.0575814 1.1685350
#> [106] 1.1866768 1.2802180 1.3925252 1.2614184 1.2557713 1.5352025 1.4243190
```

  - In order to simulate cases by date of report from a reproduction
    number trace an estimates the incubation period, reporting delay and
    generation time are required. Here defaults from `EpiNow` for
    covid-19 are
used.

<!-- end list -->

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

generation_time_defs <- EpiNow::gamma_dist_def(mean = generation_time$mean, 
                                               mean_sd = generation_time$mean_sd,
                                               sd = generation_time$sd,
                                               sd_sd = generation_time$sd_sd,
                                               max_value = generation_time$max, samples = 1)

generation_time_pdf <- c(0, EpiNow::dist_skel(n = 0:30, 
                        model = "gamma", 
                        params = generation_time_defs$params[[1]],
                        max_value = 30, 
                        dist = TRUE, cum = FALSE))

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
simulated_cases <- EpiNow::simulate_cases(rts, initial_cases = 10, initial_date = as.Date("2020-03-01"),
                                          generation_interval = generation_time_pdf,
                                          delay_def = delay_defs[10, ],
                                          incubation_def = incubation_defs[10, ],
                                          reporting_effect = c(1.6, 1.2, 1, 1, 0.8, 0.4, 1))
simulated_cases
#>            date cases reference
#>   1: 2020-03-02     5 infection
#>   2: 2020-03-03    11 infection
#>   3: 2020-03-04    14 infection
#>   4: 2020-03-05    22 infection
#>   5: 2020-03-06    26 infection
#>  ---                           
#> 327: 2020-06-16    55    report
#> 328: 2020-06-17    40    report
#> 329: 2020-06-18    37    report
#> 330: 2020-06-19    31    report
#> 331: 2020-06-20    12    report
```

### Compare approaches on simulated data

``` r
## Extract simulated infections
simulated_reports <- simulated_cases[reference == "report"][, confirm := cases][,
                                     cases := NULL][date >= as.Date("2020-03-11")]

## Rt prior
rt_prior <- list(mean = 2.6, sd = 2)

## Reconstruction via backwards sampling
sampling_cases <- nowcast_pipeline(reported_cases = simulated_reports[, import_status := "local"], 
                                   target_date = max(simulated_reports$date),
                                   delay_defs = delay_defs, 
                                   incubation_defs = incubation_defs,
                                   nowcast_lag = 0, approx_delay = TRUE)

## Non-parameteric reconstruction
non_parametric_cases <- nowcast(simulated_reports,
                                family = "negbin",
                                incubation_period = incubation_period,
                                reporting_delay = reporting_delay,
                                generation_time = generation_time, 
                                estimate_rt = TRUE,
                                rt_prior = rt_prior,
                                cores = 4, chains = 4,
                                samples = 1000, warmup = 500,
                                return_all = TRUE, model = model,
                                verbose = TRUE)
#> Running for 1500 samples and 112 time steps
#> Warning: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems
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
                                        estimate_rt = TRUE,
                                        samples = 1000, warmup = 500,
                                        cores = 4, chains = 4,
                                        return_all = TRUE, model = model, verbose = TRUE)
        
        return(list(sampling_cases, non_parametric_cases))
                                       })
#> Warning: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
#> http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
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
#> 1:    Monday 1.4576638 1.5779447 1.6073847 1.6079611 1.6379320 1.7746951
#> 2:   Tuesday 1.0467277 1.1585143 1.1823277 1.1836636 1.2088761 1.3153208
#> 3: Wednesday 0.9165018 1.0020279 1.0250456 1.0252520 1.0480361 1.1736922
#> 4:  Thursday 0.8962578 0.9883448 1.0097581 1.0102829 1.0312129 1.1450056
#> 5:    Friday 0.7131925 0.7847016 0.8024660 0.8030217 0.8209167 0.9084410
#> 6:  Saturday 0.3526790 0.3942070 0.4041530 0.4044520 0.4143264 0.4595099
#> 7:    Sunday 0.8646553 0.9438861 0.9648263 0.9653667 0.9864824 1.0800008
```

  - Recover reporting delays

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$rep_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$rep_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean          sd
#> 1:      mean 1.6087719 0.009307611
#> 2:        sd 0.6836172 0.016722758
```

  - Recover incubation period

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$inc_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$inc_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean          sd
#> 1:      mean 1.6212027 0.006118631
#> 2:        sd 0.4172562 0.006357252
```

  - Recover generation time

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$gt_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$gt_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter     mean         sd
#> 1:      mean 3.625585 0.06750276
#> 2:        sd 3.076138 0.07261434
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

  - Extract Rt
estimates

<!-- end list -->

``` r
rt_sim <- non_parametric_cases$R[, .(median = median(value), bottom = quantile(value, 0.025),
                 lower = quantile(value, 0.25), upper = quantile(value, 0.75),
                 top = quantile(value, 0.975), type = "Estimated"), by = c("date")]


rts <- data.table::data.table(
    date = seq(as.Date("2020-03-01"), 
               as.Date("2020-03-01") + lubridate::days(length(rts) - 1), by = "days"),
    median = rts, 
    type = "Truth"
  )

rt_sim <- data.table::rbindlist(list(rts, rt_sim), fill = TRUE)
rt_sim <- rt_sim[date >= as.Date("2020-03-01")][,
                 type := factor(type, levels = c("Truth", "Estimated"))]
```

  - Plot example simulation. *Note: Here we have cut-off cases prior to
    the start of March. This truncates the long tale observed in the
    sampling approach.*

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


plot_rt <- ggplot2::ggplot(rt_sim, ggplot2::aes(x = date, col = type, fill = type)) +
   ggplot2::geom_linerange(data = rt_sim[!type %in% "Truth"],
                          ggplot2::aes(ymin = bottom, ymax = top), 
                         alpha = 0.4, size = 1.5) +
 ggplot2::geom_linerange(data = rt_sim[!type %in% "Truth"],
                         ggplot2::aes(ymin = lower, ymax = upper), 
                         alpha = 1, size = 1.5) +
 ggplot2::geom_line(data =  rt_sim[type %in% "Truth"],
                    ggplot2::aes(y = median), size = 1.1, alpha = 0.7) +
 ggplot2::geom_point(data =  rt_sim[type %in% "Truth"],
                    ggplot2::aes(y = median), size = 1.1, alpha = 1, 
                    show.legend = FALSE) +
  ggplot2::geom_hline(yintercept = 1, linetype = 2) + 
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_color_manual(values = c("grey", "black")) +
  ggplot2::labs(y = "Rt", x = "Date", col = "Type")

plot_obs + 
  plot + 
  plot_rt +
  patchwork::plot_layout(ncol = 1)
```

<img src="figures/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

### Reported Covid-19 cases in the United Kingdom, United States of America and South Korea

  - Explore reporting effects by
country

<!-- end list -->

``` r
purrr::map(results, ~ .[[2]]$day_of_week[, as.list(summary(value)), by = "wday"])
#> $Austria
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.4905255 0.6709845 0.7092208 0.7115732 0.7484879 0.941267
#> 2:   Tuesday 0.7562244 0.9113330 0.9621470 0.9646137 1.0153546 1.330185
#> 3: Wednesday 0.7889382 0.9757100 1.0288871 1.0323755 1.0847925 1.342199
#> 4:  Thursday 0.7273523 0.9123487 0.9629140 0.9651904 1.0141713 1.281221
#> 5:    Friday 0.8350077 1.0340892 1.0901377 1.0935072 1.1484713 1.448119
#> 6:  Saturday 0.9507691 1.1765252 1.2380463 1.2423772 1.3056912 1.595338
#> 7:    Sunday 0.7397312 0.9366691 0.9859680 0.9903627 1.0404709 1.280013
#> 
#> $`United Kingdom`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7830954 1.0046778 1.0598064 1.0627915 1.1159369 1.406760
#> 2:   Tuesday 0.6198194 0.7735677 0.8171401 0.8213855 0.8646848 1.104779
#> 3: Wednesday 0.6499173 0.8649800 0.9131612 0.9152395 0.9615404 1.170625
#> 4:  Thursday 0.7188039 0.9102436 0.9594581 0.9642403 1.0148614 1.288638
#> 5:    Friday 0.7759915 0.9820179 1.0353002 1.0378992 1.0902040 1.421023
#> 6:  Saturday 0.7564302 1.0138360 1.0689119 1.0704922 1.1237354 1.429407
#> 7:    Sunday 0.8700791 1.0685839 1.1247688 1.1279517 1.1835092 1.400769
#> 
#> $`United States of America`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7710862 0.8779974 0.9031431 0.9043710 0.9296599 1.074304
#> 2:   Tuesday 0.8088956 0.9152615 0.9411576 0.9431975 0.9695850 1.172020
#> 3: Wednesday 0.8060068 0.8957220 0.9220903 0.9226625 0.9482166 1.071084
#> 4:  Thursday 0.8211490 0.9237426 0.9525485 0.9535643 0.9825141 1.143120
#> 5:    Friday 0.9143059 1.0379227 1.0684542 1.0684111 1.0983809 1.231315
#> 6:  Saturday 0.9705658 1.1105940 1.1434697 1.1458923 1.1788160 1.342993
#> 7:    Sunday 0.8833059 1.0320429 1.0606333 1.0619012 1.0922653 1.269133
#> 
#> $Russia
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7976414 0.9640627 1.0047855 1.0074416 1.0485824 1.225616
#> 2:   Tuesday 0.7584148 0.9277684 0.9663910 0.9694535 1.0074213 1.194287
#> 3: Wednesday 0.7421135 0.9064927 0.9495519 0.9499456 0.9902946 1.222153
#> 4:  Thursday 0.8234741 0.9788424 1.0213422 1.0234099 1.0633884 1.266413
#> 5:    Friday 0.7915128 0.9661205 1.0083554 1.0104975 1.0522869 1.294366
#> 6:  Saturday 0.8727061 1.0344841 1.0828686 1.0843704 1.1306687 1.418592
#> 7:    Sunday 0.7602219 0.9134921 0.9534680 0.9548815 0.9945370 1.207848
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

<img src="figures/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

  - Prepare R estimates for plotting

<!-- end list -->

``` r
country_r <- purrr::map2(results, names(results), ~ .x[[2]]$R[,
                  .(median = median(value), bottom = quantile(value, 0.025),
                    lower = quantile(value, 0.25), upper = quantile(value, 0.75),
                    top = quantile(value, 0.975), type = "Estimated"), by = c("date")][,region := .y])

country_r <- data.table::rbindlist(country_r)[date >= as.Date("2020-03-01")]
```

  - Plot R estimates (**note that the reporting delays are not based on
    data and so these results are for model checking only**).

<!-- end list -->

``` r
plot_rt_countries <- 
  ggplot2::ggplot(country_r , ggplot2::aes(x = date, col = region, fill = region)) +
   ggplot2::geom_linerange(ggplot2::aes(ymin = bottom, ymax = top), 
                         alpha = 0.4, size = 1.2) +
 ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), 
                         alpha = 0.8, size = 1.2) +
  ggplot2::geom_hline(yintercept = 1, linetype = 2) + 
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(y = "Rt", x = "Date", col = "Type") +
  facet_wrap(~ region) +
  ggplot2::theme(legend.position = "none")

plot_rt_countries
```

<img src="figures/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

``` r
ggsave("figures/rt_plot.png", plot_rt_countries, dpi = 400, width = 12, height = 12)
```

## Discussion

## References
