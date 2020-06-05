source("nowcast.R")

rts <- c(rep(2, 20), (2 - 1:15 * 0.1), rep(0.5, 10), (0.5 + 1:7 * 0.1), rep(1.2, 10))
generation_interval <- rowMeans(EpiNow::covid_generation_times)

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

simulated_cases <- EpiNow::simulate_cases(rts, initial_cases = 100 , initial_date = as.Date("2020-03-01"),
                                          generation_interval = generation_interval, delay_def = delay_defs[1, ],
                                          incubation_def = incubation_defs[1, ],
                                          reporting_effect = c(1.4, rep(1, 4), 0.8, 0.8))

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

epitest <- sampling_cases[type == "infection_upscaled",.(date,confirm = cases)][,.(confirm = round(median(confirm))),by = "date"]

## Non-parameteric reconstruction
model <- rstan::stan_model("nowcast.stan")
non_parametric_cases <- nowcast(simulated_reports,
                                prior = epitest,
                                family = "poisson", incubation_period = incubation_period,
                                reporting_delay = reporting_delay,
                                generation_interval = generation_interval, cores = 1, chains = 1,
                                samples = 500, return_all = TRUE, model = model, verbose = TRUE)



res <- rstan::extract(non_parametric_cases$fit)

data.table(cases = apply(res$imputed_infections, 2, median),
           reports = simulated_reports$confirm,
           infections = subset(simulated_cases, reference == "infection")$cases,
           date = simulated_reports$date) %>%
  ggplot(aes(x = date)) +
  geom_bar(aes(y = infections), stat = "identity") +
  geom_line(aes(y = reports), col = "purple") +
  geom_line(aes(y = cases), col = "orange")
