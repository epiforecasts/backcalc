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

<!-- end list -->

``` r
model
#> S4 class stanmodel 'nowcast' coded as follows:
#> functions {
#>   // Convolve a pdf and case vector using matrix multiplication
#>   vector convolve(vector cases, vector pdf) {
#>     int t = num_elements(cases);
#>     matrix[t, t] delay_mat = rep_matrix(0, t, t);
#>     int max_pdf = num_elements(pdf) + 1;
#>     row_vector[max_pdf] row_pdf = to_row_vector(append_row(pdf, 0.0));
#>     vector[t] convolved_cases;
#>     
#>     for (s in 1:t) {
#>       int max_length = min(s, max_pdf);
#>       delay_mat[s, (s - max_length + 1):s] = row_pdf[(max_pdf - max_length + 1):max_pdf];
#>     }
#>   
#>    convolved_cases = delay_mat * to_vector(cases);
#> 
#>    convolved_cases[1] = 0.000001;
#>    
#>    return convolved_cases;
#>   }
#> 
#>   real discretised_lognormal_pmf(int y, real mu, real sigma) {
#> 
#>     return((lognormal_cdf(y, mu, sigma) - lognormal_cdf(y - 1, mu, sigma)));
#>   }
#> }
#> 
#> 
#> data {
#>   int t; // number of time steps
#>   int max_rep; 
#>   int max_inc; 
#>   int day_of_week[t];
#>   int <lower = 0> cases[t];
#>   vector<lower = 0>[t] shifted_cases; 
#>   real inc_mean_sd;                  // prior sd of mean incubation period
#>   real inc_mean_mean;                // prior mean of mean incubation period
#>   real inc_sd_mean;                  // prior sd of sd of incubation period
#>   real inc_sd_sd;                    // prior sd of sd of incubation period
#>   real rep_mean_mean;                // prior mean of mean reporting delay
#>   real rep_mean_sd;                  // prior sd of mean reporting delay
#>   real rep_sd_mean;                  // prior mean of sd of reporting delay
#>   real rep_sd_sd;                    // prior sd of sd of reporting delay
#>   int model_type; //Type of model: 1 = Poisson otherwise negative binomial
#> }
#> 
#> transformed data{
#>   int<lower = 0> weekly_cases[t];
#>   int<lower = 0> cum_cases[t];
#>   
#>   //Calculate weekly cases
#>   cum_cases[1] = cases[1];
#>   weekly_cases[1] = cases[1];
#>   for (s in 2:t) { 
#>     cum_cases[s] = cum_cases[s - 1] + cases[s];
#>     weekly_cases[s] = cum_cases[s] - cum_cases[max(1, s - 7)];
#>   }
#>   
#> }
#> parameters{
#>   vector<lower = 0>[t] noise;
#>   real <lower = 0> inc_mean;         // mean of incubation period
#>   real <lower = 0> inc_sd;           // sd of incubation period
#>   real <lower = 0> rep_mean;         // mean of reporting delay
#>   real <lower = 0> rep_sd;           // sd of incubation period
#>   real<lower = 0> phi; 
#>   simplex[7] day_of_week_eff_raw;
#> }
#> 
#> transformed parameters {
#>   vector[max_rep] rev_delay;
#>   vector[max_inc] rev_incubation;
#>   vector<lower = 0>[t] infections;
#>   vector<lower = 0>[t] onsets;
#>   vector<lower = 0>[t] reports;
#>   vector<lower = 0>[t] cum_reports;
#>   vector<lower = 0>[t] weekly_reports;
#>   vector[7] day_of_week_eff;
#>   
#>   day_of_week_eff = day_of_week_eff_raw * 7;
#> 
#>   //Reverse the distributions to allow vectorised access
#>     for (j in 1:max_rep) {
#>       rev_delay[j] =
#>         discretised_lognormal_pmf(max_rep - j + 1, inc_mean, inc_sd);
#>         }
#>    
#>     for (j in 1:max_inc) {
#>       rev_incubation[j] =
#>         discretised_lognormal_pmf(max_inc - j + 1, rep_mean, rep_sd);
#>     }
#> 
#>   //Generation infections from median shifted cases and non-parameteric noise
#>   infections = shifted_cases .* noise;
#> 
#>   // Onsets from infections
#>   onsets = convolve(infections, rev_incubation);
#>      
#>   // Reports from onsets
#>   reports = convolve(onsets, rev_delay);
#>      
#>  //Calculate Cumulative reports
#>   cum_reports = cumulative_sum(reports);
#>   
#>   for (s in 1:t) {
#>     //Calculate weekly reports
#>     weekly_reports[s] = s == 1 ? cum_reports[1] : cum_reports[s] - cum_reports[max(1, s - 7)];
#>     // Add reporting effects
#>     reports[s] *= day_of_week_eff[day_of_week[s]];
#>     }
#> }
#> 
#> model {
#> 
#>   // Reporting overdispersion
#>   phi ~ exponential(1);
#> 
#>   // Noise on median shift
#>   for (i in 1:t) {
#>     noise[i] ~ normal(1, 0.4) T[0,];
#>   }
#>   
#>   // daily cases given reports
#>   if (model_type == 1) {
#>     target +=  poisson_lpmf(cases | reports);
#>   }else{
#>     target += neg_binomial_2_lpmf(cases | reports, phi);
#>   }
#>   
#>   // weekly cases given weekly reports
#>   target += poisson_lpmf(weekly_cases[7:t] | weekly_reports[7:t]);
#>   
#>   // penalised priors
#>   target += normal_lpdf(inc_mean | inc_mean_mean, inc_mean_sd) * t;
#>   target += normal_lpdf(inc_sd | inc_sd_mean, inc_sd_sd) * t;
#>   target += normal_lpdf(rep_mean | rep_mean_mean, rep_mean_sd) * t;
#>   target += normal_lpdf(rep_sd | rep_sd_mean, rep_sd_sd) * t;
#> }
#>   
#> generated quantities {
#>   int imputed_infections[t];
#>  
#>   imputed_infections = poisson_rng(infections);
#> 
#> }
```

## Analysis

### Simulate data

  - Define a realistic basic reproduction number estimate that starts at
    2, decreases linearly to 0.5, remains constant, increases linearly
    to 1.2 and then again remains constant.

<!-- end list -->

``` r
## Define an initial rt vector 
rts <- c(rep(2, 20), (2 - 1:15 * 0.1), rep(0.5, 10), (0.5 + 1:7 * 0.1), rep(1.2, 20))
rts
#>  [1] 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
#> [20] 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.5 0.5 0.5
#> [39] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.2 1.2 1.2 1.2 1.2
#> [58] 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2
```

  - In order to simulate cases by date of infection from a reproduction
    number trace an estimate of the generation interval is required. We
    use the current default estimate for Covid-19 from `EpiNow`.

<!-- end list -->

``` r
## Use the mean default generation interval for covid
generation_interval <- rowMeans(EpiNow::covid_generation_times)

generation_interval
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
simulated_cases <- EpiNow::simulate_cases(rts, initial_cases = 100 , initial_date = as.Date("2020-03-01"),
                                          generation_interval = generation_interval, delay_def = delay_defs[1, ],
                                          incubation_def = incubation_defs[1, ],
                                          reporting_effect = c(1.4, rep(1, 4), 0.8, 0.8))
simulated_cases
#>            date cases reference
#>   1: 2020-03-02    35 infection
#>   2: 2020-03-03    54 infection
#>   3: 2020-03-04    70 infection
#>   4: 2020-03-05    73 infection
#>   5: 2020-03-06    92 infection
#>  ---                           
#> 205: 2020-05-07   724    report
#> 206: 2020-05-08   744    report
#> 207: 2020-05-09   624    report
#> 208: 2020-05-10   648    report
#> 209: 2020-05-11  1163    report
```

### Compare approaches on simulated data

``` r
## Extract simulated infections
simulated_reports <- simulated_cases[reference == "report"][, confirm := cases][, cases := NULL]

# Median covid generation interval - used to smooth prior cases by date of infection
generation_interval <- rowMeans(EpiNow::covid_generation_times)
generation_interval <- sum(!(cumsum(generation_interval) > 0.5)) + 1   

## Reconstruction via backwards sampling
sampling_cases <- nowcast_pipeline(reported_cases = simulated_reports[, import_status := "local"], 
                                   target_date = max(simulated_reports$date),
                                   delay_defs = delay_defs,
                                   incubation_defs = incubation_defs,
                                   nowcast_lag = 0, approx_delay = TRUE)

## Non-parameteric reconstruction
non_parametric_cases <- nowcast(simulated_reports,
                                family = "negbin", incubation_period = incubation_period,
                                reporting_delay = reporting_delay,
                                generation_interval = 7, cores = 4, chains = 4,
                                samples = 2000, warmup = 1000,
                                return_all = TRUE, model = model,
                                verbose = TRUE)
#> Running for 3000 samples and 71 time steps
#> Warning: There were 5467 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
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
                                           approx_delay = TRUE)
        
        message("Non-parametric nowcasting for: ", country)
        ## Non-parametric reconstruction
        non_parametric_cases <- nowcast(cases,
                                        family = "negbin",
                                         incubation_period = incubation_period,
                                         reporting_delay = reporting_delay,
                                         generation_interval = 7, 
                                         samples = 2000, warmup = 1000,
                                         cores = 4, chains = 4,
                                         return_all = TRUE, model = model)
        
        return(list(sampling_cases, non_parametric_cases))
                                       })
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> http://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: There were 7504 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
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
#> Warning: There were 8000 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
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
#> Warning: There were 7988 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
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

names(results) <- countries
```

## Results

### Simulated data

  - Reporting effects.

<!-- end list -->

``` r
non_parametric_cases$day_of_week[, as.list(summary(value)), by = "wday"]
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 1.0612804 1.3379198 1.3979104 1.3999619 1.4601711 1.905373
#> 2:   Tuesday 0.6893334 0.9407431 0.9869856 0.9897406 1.0359467 1.323394
#> 3: Wednesday 0.7594623 0.9431678 0.9906221 0.9922409 1.0384694 1.342695
#> 4:  Thursday 0.7494983 0.9266544 0.9709438 0.9741616 1.0184370 1.268216
#> 5:    Friday 0.7619107 0.9822665 1.0303073 1.0344177 1.0831569 1.395764
#> 6:  Saturday 0.5860795 0.7632764 0.8020968 0.8048721 0.8443212 1.101315
#> 7:    Sunday 0.5934522 0.7638819 0.8024100 0.8046052 0.8418437 1.216177
```

  - Recover reporting delays

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$rep_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$rep_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter     mean         sd
#> 1:      mean 1.414904 0.01018756
#> 2:        sd 0.882964 0.01997711
```

  - Recover incubation period

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$inc_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$inc_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean          sd
#> 1:      mean 1.4960972 0.006826419
#> 2:        sd 0.4523496 0.007937645
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

out <- data.table::rbindlist(list(sampling_cases, non_parametric_infections), fill = TRUE)
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
plot_data <- simulated_cases[date >= as.Date("2020-03-01")][,
                             type := factor(type, levels = c("Simulated infections",
                                                             "Simulated reported cases",
                                                             "Non-parametric",
                                                             "Sampled"))]
plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = date, col = type, fill = type)) +
  ggplot2::geom_col(data = plot_data[type %in% "Simulated infections"],
                    ggplot2::aes(y = median), fill = "grey", col = "white", alpha = 0.8) +
  ggplot2::geom_line(data = plot_data[type %in% "Simulated reported cases"],
                      ggplot2::aes(y = median), size = 1.1) +
  ggplot2::geom_linerange(data = plot_data[!type %in% "Simulated infections"],
                          ggplot2::aes(ymin = bottom, ymax = top), 
                         alpha = 0.4, size = 1.5) +
 ggplot2::geom_linerange(data = plot_data[!type %in% "Simulated infections"],
                         ggplot2::aes(ymin = lower, ymax = upper), 
                         alpha = 0.6, size = 1.5) +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(y = "Cases", x = "Date", col = "Type")

plot
```

<img src="figures/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

### Reported Covid-19 cases in the United Kingdom, United States of America and South Korea

  - Explore reporting effects by
country

<!-- end list -->

``` r
purrr::map(results, ~ .[[2]]$day_of_week[, as.list(summary(value)), by = "wday"])
#> $Austria
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.4669745 0.6689007 0.7206945 0.7257355 0.7765343 1.114203
#> 2:   Tuesday 0.6382305 0.8742759 0.9398685 0.9443072 1.0083197 1.414285
#> 3: Wednesday 0.6668444 0.9571783 1.0234750 1.0294105 1.0964797 1.493047
#> 4:  Thursday 0.6334571 0.9009634 0.9662293 0.9720001 1.0377843 1.523515
#> 5:    Friday 0.7669587 1.0156918 1.0888544 1.0946703 1.1689886 1.556046
#> 6:  Saturday 0.8266859 1.1656979 1.2447583 1.2506084 1.3292541 1.737623
#> 7:    Sunday 0.6253702 0.9092676 0.9771725 0.9832678 1.0510455 1.445325
#> 
#> $`United Kingdom`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7388452 1.0083463 1.0885986 1.0948306 1.1731053 1.596544
#> 2:   Tuesday 0.5233288 0.7665553 0.8272902 0.8332671 0.8953187 1.276188
#> 3: Wednesday 0.6129372 0.8426309 0.9104190 0.9161048 0.9814633 1.344184
#> 4:  Thursday 0.5782797 0.8860045 0.9557673 0.9611856 1.0290832 1.493742
#> 5:    Friday 0.6916695 0.9559284 1.0269914 1.0331561 1.1053419 1.567411
#> 6:  Saturday 0.6863300 0.9763633 1.0483026 1.0552607 1.1243660 1.643847
#> 7:    Sunday 0.7350129 1.0248281 1.0994692 1.1061950 1.1807101 1.689201
#> 
#> $`United States of America`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7560342 0.8848475 0.9195014 0.9215024 0.9559357 1.133765
#> 2:   Tuesday 0.7744194 0.9194093 0.9554675 0.9574450 0.9929311 1.192861
#> 3: Wednesday 0.7793832 0.8968170 0.9325586 0.9341981 0.9702118 1.163909
#> 4:  Thursday 0.7428890 0.9058808 0.9430109 0.9442855 0.9806517 1.182489
#> 5:    Friday 0.8394010 1.0089842 1.0495585 1.0511028 1.0914278 1.291866
#> 6:  Saturday 0.9318711 1.0879335 1.1300293 1.1338715 1.1753839 1.389428
#> 7:    Sunday 0.8361640 1.0180410 1.0563070 1.0575948 1.0961447 1.296545
#> 
#> $Russia
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.5128884 0.8444914 0.9339511 0.9423040 1.0272943 1.784059
#> 2:   Tuesday 0.5466355 0.8736421 0.9610450 0.9716258 1.0591236 1.901598
#> 3: Wednesday 0.5129630 0.7561033 0.8372503 0.8471555 0.9266051 1.595401
#> 4:  Thursday 0.6175618 0.9888030 1.0913885 1.1027854 1.2051573 2.024116
#> 5:    Friday 0.5622044 0.8811856 0.9682046 0.9787580 1.0648672 1.580794
#> 6:  Saturday 0.6651156 1.0959534 1.2099567 1.2215629 1.3360609 2.335221
#> 7:    Sunday 0.5260889 0.8386355 0.9231650 0.9358084 1.0221042 1.604754
```

  - Prepare data for
plotting

<!-- end list -->

``` r
summarised_cases <- data.table::copy(reported_cases)[, median := confirm][, type := "Reported Cases"]


summarised_results <- purrr::map2(results, names(results),
                                  ~ summarise_nowcasting_approaches(.x[[1]], .x[[2]])[, region := .y])

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
                         alpha = 0.4, size = 1) +
 ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper), 
                         alpha = 0.8, size = 1) +
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
