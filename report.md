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
#>     int max_pdf = num_elements(pdf);
#>     row_vector[max_pdf] row_pdf = to_row_vector(pdf);
#>     vector[t] convolved_cases;
#>     
#>     for (s in 1:t) {
#>       int max_length = min(s, max_pdf);
#>       delay_mat[s, (s - max_length + 1):s] = row_pdf[(max_pdf - max_length + 1):max_pdf];
#>     }
#>   
#>    convolved_cases = delay_mat * to_vector(cases);
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
#>   vector[7] day_of_week_eff;
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
#> 
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
#>     weekly_reports[s] = cum_reports[s] - cum_reports[max(1, s - 7)];
#>     // Add reporting effects
#>     reports[s] *= day_of_week_eff[day_of_week[s]];
#>     }
#> }
#> 
#> model {
#>   // Week effect
#>   for (j in 1:7) {
#>     day_of_week_eff[j] ~ normal(1, 0.2) T[0,];
#>   }
#>   
#>   // Reporting overdispersion
#>   phi ~ exponential(1);
#> 
#>   // Noise on median shift
#>   for (i in 1:t) {
#>     noise[i] ~ normal(1, 0.2) T[0,];
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
#>   if (model_type == 1) {
#>     imputed_infections = poisson_rng(infections);
#>   }else{
#>     imputed_infections = neg_binomial_2_rng(infections, phi);
#>   }
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
#>   1: 2020-03-02    50 infection
#>   2: 2020-03-03    54 infection
#>   3: 2020-03-04    67 infection
#>   4: 2020-03-05   101 infection
#>   5: 2020-03-06   125 infection
#>  ---                           
#> 207: 2020-05-07   971    report
#> 208: 2020-05-08  1024    report
#> 209: 2020-05-09   819    report
#> 210: 2020-05-10   876    report
#> 211: 2020-05-11  1657    report
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
                                generation_interval = generation_interval, cores = 4, chains = 4,
                                samples = 1000, return_all = TRUE, model = model,
                                verbose = TRUE)
#> Running for 2000 samples and 72 time steps
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
                                         generation_interval = generation_interval, 
                                         samples = 1000, cores = 4, chains = 4,
                                         return_all = TRUE, model = model)
        
        return(list(sampling_cases, non_parametric_cases))
                                       })

names(results) <- countries
```

## Results

### Simulated data

  - Reporting effects.

<!-- end list -->

``` r
non_parametric_cases$day_of_week[, as.list(summary(value)), by = "wday"]
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 1.1190654 1.3301449 1.3856289 1.3879045 1.4431756 1.740700
#> 2:   Tuesday 0.8229757 0.9748764 1.0213324 1.0243398 1.0703676 1.324536
#> 3: Wednesday 0.8312469 0.9691534 1.0134297 1.0143509 1.0572296 1.304087
#> 4:  Thursday 0.8114591 0.9763905 1.0210896 1.0235535 1.0693038 1.344542
#> 5:    Friday 0.8617285 1.0056231 1.0486717 1.0524875 1.0962570 1.310068
#> 6:  Saturday 0.6489438 0.8215581 0.8597265 0.8650650 0.9040898 1.172037
#> 7:    Sunday 0.7030686 0.8436684 0.8843590 0.8876099 0.9285475 1.128053
```

  - Recover reporting delays

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$rep_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$rep_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean         sd
#> 1:      mean 1.3782441 0.02487498
#> 2:        sd 0.5689783 0.30851332
```

  - Recover incubation period

<!-- end list -->

``` r
data.table::rbindlist(list(
  non_parametric_cases$inc_mean[, .(parameter = "mean", mean = mean(value), sd = sd(value))],
  non_parametric_cases$inc_sd[, .(parameter = "sd", mean = mean(value), sd = sd(value))]
))
#>    parameter      mean          sd
#> 1:      mean 1.4990592 0.007714307
#> 2:        sd 0.4329428 0.012153610
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
#> 1:    Monday 0.5927318 0.7733996 0.8346302 0.8397165 0.8977349 1.263271
#> 2:   Tuesday 0.6910251 0.9539920 1.0191157 1.0240155 1.0882138 1.432329
#> 3: Wednesday 0.7147801 1.0365139 1.1092983 1.1132127 1.1843984 1.608283
#> 4:  Thursday 0.7603698 0.9684866 1.0358921 1.0391132 1.1060063 1.403142
#> 5:    Friday 0.7801284 1.0498981 1.1199424 1.1218265 1.1882261 1.556453
#> 6:  Saturday 0.7927248 1.1538379 1.2237964 1.2264513 1.2952796 1.623405
#> 7:    Sunday 0.7139299 0.9511790 1.0198030 1.0247374 1.0923124 1.426632
#> 
#> $`United Kingdom`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7512457 1.0458641 1.1199028 1.1262421 1.2026619 1.611525
#> 2:   Tuesday 0.5165016 0.8348616 0.9037527 0.9104707 0.9777852 1.517989
#> 3: Wednesday 0.6670868 0.9075144 0.9755122 0.9813382 1.0500160 1.450841
#> 4:  Thursday 0.7001166 0.9516667 1.0217846 1.0272663 1.0963654 1.454746
#> 5:    Friday 0.7782634 1.0108068 1.0854092 1.0904367 1.1633567 1.603845
#> 6:  Saturday 0.7511191 1.0299033 1.1025894 1.1054621 1.1747001 1.559523
#> 7:    Sunday 0.8159668 1.0591894 1.1377700 1.1432255 1.2226046 1.594999
#> 
#> $`United States of America`
#>         wday      Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> 1:    Monday 0.7664894 0.9108803 0.9503306 0.9518160 0.9900741 1.176280
#> 2:   Tuesday 0.7810508 0.9447378 0.9842269 0.9882454 1.0305802 1.227700
#> 3: Wednesday 0.7835703 0.9360037 0.9756947 0.9762278 1.0140231 1.193930
#> 4:  Thursday 0.8134383 0.9478350 0.9861615 0.9881642 1.0258865 1.236909
#> 5:    Friday 0.8788655 1.0425256 1.0871017 1.0883173 1.1308972 1.342752
#> 6:  Saturday 0.9301455 1.0918874 1.1389828 1.1434254 1.1935769 1.395370
#> 7:    Sunday 0.8476215 1.0357152 1.0750295 1.0782591 1.1195598 1.354348
#> 
#> $Russia
#>         wday      Min.   1st Qu.    Median      Mean  3rd Qu.     Max.
#> 1:    Monday 0.7173941 0.9959770 1.0556791 1.0603599 1.121556 1.418810
#> 2:   Tuesday 0.6806643 0.9224963 0.9815973 0.9862169 1.043192 1.374461
#> 3: Wednesday 0.7129118 0.9516524 1.0106674 1.0146315 1.070632 1.468565
#> 4:  Thursday 0.7686629 1.0418027 1.1107470 1.1149173 1.182414 1.563897
#> 5:    Friday 0.7733246 1.0115712 1.0755805 1.0785990 1.139767 1.511377
#> 6:  Saturday 0.8384414 1.0766263 1.1408221 1.1472256 1.210748 1.606291
#> 7:    Sunday 0.7136652 0.9641468 1.0227347 1.0266152 1.084384 1.379870
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
