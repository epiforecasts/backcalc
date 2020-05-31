Evaluating approaches to backcalculating cases counts by date of
infection from cases counts by date of report
================

**Authors:** EpiForecasts, CMMID Covid working group, Sebastian Funk

## Summary

## Introduction

## Methods

## Model

  - Orginal implementation here:
    <https://github.com/jhellewell14/rtconfirm/blob/master/R/backcalc.R>

  - This implementation uses median shifted reported cases (smoothed
    using a rolling average over the width of the median generation
    interval) as a prior and then fits independent gaussian noise on top
    of this. For future cases (i.e with no data to shift to into the
    last reported case count is used).

  - Weekend and monday effects are included as multiplicative terms

  - Reporting delays and incubation periods are passed in so uncertainty
    can only be generated by passing in multiple samples (causes a
    non-linear slow down) or fitting the model multiple times (each
    model run is ~ 3 - 4 seconds).

  - An alternative would be jointly fitting the delays and case counts
    but I can’t see how fitting the delays in the model in any way can
    be a good thing as I don’t think it will be identifiable.

## Analysis

  - Packages

<!-- end list -->

``` r
library(EpiNow)
library(data.table)
```

  - Simulate data

<!-- end list -->

``` r
## Define an initial rt vector 
rts <- c(rep(2, 20), (2 - 1:15 * 0.1), rep(0.5, 10))
rts
```

    ##  [1] 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
    ## [20] 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.5 0.5 0.5
    ## [39] 0.5 0.5 0.5 0.5 0.5 0.5 0.5

``` r
## Use the mean default generation interval for covid
generation_interval <- rowMeans(EpiNow::covid_generation_times)

## Sample a report delay as a lognormal
delay_def <- EpiNow::lognorm_dist_def(mean = 5, mean_sd = 1,
                                      sd = 3, sd_sd = 1, max_value = 30,
                                      samples = 1, to_log = TRUE)


## Sample a incubation period (again using the default for covid)
incubation_def <- EpiNow::lognorm_dist_def(mean = EpiNow::covid_incubation_period[1, ]$mean,
                                          mean_sd = EpiNow::covid_incubation_period[1, ]$mean_sd,
                                          sd = EpiNow::covid_incubation_period[1, ]$sd,
                                          sd_sd = EpiNow::covid_incubation_period[1, ]$sd_sd,
                                          max_value = 30, samples = 1)

## Simulate cases with a decrease in reporting at weekends and an incease on Monday                                     
simulated_cases <- simulate_cases(rts, initial_cases = 100 , initial_date = as.Date("2020-03-01"),
                    generation_interval = generation_interval, delay_def = delay_def,
                   incubation_def = incubation_def, reporting_effect = c(1.1, rep(1, 4), 0.95, 0.95))

simulated_cases
```

    ##            date cases reference
    ##   1: 2020-03-02    53 infection
    ##   2: 2020-03-03    59 infection
    ##   3: 2020-03-04    79 infection
    ##   4: 2020-03-05   107 infection
    ##   5: 2020-03-06   131 infection
    ##  ---                           
    ## 128: 2020-04-10  2258    report
    ## 129: 2020-04-11  2119    report
    ## 130: 2020-04-12  2252    report
    ## 131: 2020-04-13  2016    report
    ## 132: 2020-04-14  1910    report

  - Fit model and compare to simulated data (+ to current sampling) and
    ability to recover reporting effects.

  - Fit to observed data: Austria (see `nowcast.R`) example.

## Results

## Discussion

## References