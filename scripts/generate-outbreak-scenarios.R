
# Packages ----------------------------------------------------------------

library(data.table)
library(EpiNow)
library(purrr)
library(future)
library(future.apply)
library(here)


# Set up parallel processing ----------------------------------------------

future::plan("multisession")

# Functions ---------------------------------------------------------------

source(here::here("R", "generate_simulations.R"))

# Initial cases -----------------------------------------------------------

case_scenarios <- data.table::data.table(
  index = 1,
  case_scenario = c("10", "100", "1000"),
  case = c(10, 100, 100)
)

# Rt scenarios ------------------------------------------------------------

rt_scenarios <- data.table::data.table(
  index = 1,
  rt_scenario = c("Stable", "Stable subcriticial", "Step change", "Linear change", "Complex"),
  rt = list(
    rep(1.2, 50),
    rep(0.9, 50),
    c(rep(3, 10), rep(2, 10), rep(1.4, 10), rep(1.2, 10), rep(0.8, 10)),
    c(rep(2, 10), (2 - 1:15 * 0.1), rep(0.5, 10), (0.5 + 1:15 * 0.1)),
    c(rep(3, 10), rep(1, 10), (1 + 1:10 * 0.1), rep(0.5, 5), (0.5 + 1:5 * 0.1),
      rep(1, 10))
  )
)

# Distributions -----------------------------------------------------------

## Covid incubation period - see ?EpiNow::covid_incubation_period
incubation_period <- list(mean = EpiNow::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow::covid_incubation_period[1, ]$sd_sd,
                          max = 30)

## Covid generation time estimate - see ?EpiNow::covid_generation_times_summary
generation_time <- list(mean = EpiNow::covid_generation_times_summary[1, ]$mean,
                        mean_sd = EpiNow::covid_generation_times_summary[1, ]$mean_sd,
                        sd = EpiNow::covid_generation_times_summary[1, ]$sd,
                        sd_sd = EpiNow::covid_generation_times_summary[1, ]$sd_sd,
                        max = 30)


delay_scenarios <- data.table::data.table(
  index = 1,
  delay_scenario = c("2 days", "5 days", "10 days"),
  reporting_delay = list(
    list(mean = log(2),
         mean_sd = log(1),
         sd = log(1),
         sd_sd = log(0.5),
         max = 30),
    list(mean = log(5),
         mean_sd = log(1),
         sd = log(2),
         sd_sd = log(1),
         max = 30),
    list(mean = log(10),
         mean_sd = log(2),
         sd = log(4),
         sd_sd = log(2),
         max = 30)
  ),
  incubation_period = rep(list(incubation_period), 3),
  generation_time = rep(list(generation_time), 3))



# Reporting scenarios -----------------------------------------------------

reporting_scenarios <- data.table::data.table(
  index = 1,
  weekly_report_scenario = c("None", "Small", "Large"),
  reporting_effect = list(rep(1, 7),
                          c(1.2, 1.1, 1, 1, 0.9, 0.8, 1),
                          c(2, 0.5, 0.5, 1.5, 0.75, 0.75, 1)))




# Combine scenarios -------------------------------------------------------

scenarios <- data.table::merge.data.table(rt_scenarios, delay_scenarios,
                                          by = "index", all = TRUE, allow.cartesian = TRUE)

scenarios <- data.table::merge.data.table(scenarios, reporting_scenarios,
                                          by = "index", all = TRUE, allow.cartesian = TRUE)

scenarios <- data.table::merge.data.table(scenarios, case_scenarios, 
                                          by = "index", all = TRUE, allow.cartesian = TRUE)


# Reorder and add indexing to scenarios -----------------------------------

scenarios <- scenarios[, index := NULL][, scenario := 1:.N]

data.table::setcolorder(scenarios, c("scenario", "rt_scenario", "case_scenario", "delay_scenario", "weekly_report_scenario"))


saveRDS(scenarios, here::here("data", "scenarios.rds"))


# Run scenario simulations ------------------------------------------------


simulations <- future.apply::future_lapply(split(scenarios, by = "scenario"),
                                           function(scenario) {
                                             generate_simulations(rts = scenario$rt[[1]],
                                                                  initial_cases = scenario$case[[1]], 
                                                                  generation_time = scenario$generation_time[[1]],
                                                                  incubation_period = scenario$incubation_period[[1]],
                                                                  reporting_delay = scenario$reporting_delay[[1]],
                                                                  reporting_effect = scenario$reporting_effect[[1]],
                                                                  samples = 10)},
                                           future.scheduling = Inf)


simulations <- data.table::rbindlist(simulations, idcol = "scenario")

data.table::fwrite(simulations, here::here("data", "simulations.csv"))
