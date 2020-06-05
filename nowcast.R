#' @export
#' @importFrom rstan sampling extract
#' @importFrom data.table data.table copy merge.data.table setorder rbindlist setDTthreads melt .N
#' @importFrom purrr tranpose
#' @importFrom future.apply future_lapply
#' @importFrom lubridate wday
#' @examples
#' reported_cases <- NCoVUtils::get_ecdc_cases(countries = "Russia")
#' reported_cases <- NCoVUtils::format_ecdc_data(reported_cases)
#' reported_cases <- data.table::as.data.table(reported_cases)[, confirm := cases][, cases := NULL]
#'  
#' incubation_period <- list(mean = EpiNow::covid_incubation_period[1, ]$mean,
#'                           mean_sd = EpiNow::covid_incubation_period[1, ]$mean_sd,
#'                           sd = EpiNow::covid_incubation_period[1, ]$sd,
#'                           sd_sd = EpiNow::covid_incubation_period[1, ]$sd_sd,
#'                           max = 30)
#'                    
#' reporting_delay <- list(mean = log(5),
#'                         mean_sd = log(2),
#'                         sd = log(2),
#'                         sd_sd = log(1.5),
#'                         max = 30)
#'
#'  generation_interval <- rowMeans(EpiNow::covid_generation_times)
#'  generation_interval <- sum(!(cumsum(generation_interval) > 0.5)) + 1   
#'   
#'  ## Compile model
#'  model <- rstan::stan_model("nowcast.stan")
#' 
#' ## Run model
#' out <- nowcast(reported_cases, family = "poisson",
#'                generation_interval = generation_interval,
#'                incubation_period = incubation_period,
#'                reporting_delay = reporting_delay,
#'                model = model,
#'                cores = 4, chains = 4,
#'                verbose = TRUE, return_all = TRUE
#'                )
#'
#' out                                   
nowcast <- function(reported_cases, family = "poisson",
                    incubation_period, reporting_delay,
                    generation_interval, model,
                    cores = 1,
                    chains = 2,
                    samples = 1000,
                    warmup = 1000,
                    return_all = FALSE,
                    verbose = FALSE,
                    prior){
  
  suppressMessages(data.table::setDTthreads(threads = 1))
  
  # Make sure there are no missing dates and order cases --------------------
  reported_cases_grid <- data.table::copy(reported_cases)[, .(date = seq(min(date), max(date), by = "days"))]

  nowcast <-  data.table::merge.data.table(
    reported_cases , reported_cases_grid, 
    by = c("date"), all.y = TRUE)
  
  reported_cases <-  reported_cases[is.na(confirm), confirm := 0 ][,.(date = date, confirm)]
  reported_cases <- data.table::setorder(reported_cases, date)
  
  ## Filter out 0 reported cases
  reported_cases <- reported_cases[, cum_cases := cumsum(confirm)][cum_cases > 0][, cum_cases := NULL]

# Estimate the mean delay -----------------------------------------------
  
  mean_shift <- incubation_period$mean + reporting_delay$mean

# Add the mean delay and incubation period on as 0 case days ------------

  reported_cases <- data.table::rbindlist(list(
    data.table::data.table(date = seq(min(reported_cases$date) - mean_shift - generation_interval, 
                                      min(reported_cases$date) - 1, by = "days"),
                           confirm = 0),
    reported_cases
  ))  
  
  # Calculate smoothed prior cases ------------------------------------------

  shifted_reported_cases <- prior
  
  ##Drop median generation interval initial values
  # shifted_reported_cases <- shifted_reported_cases[-(1:generation_interval)]
  reported_cases <- reported_cases[-(1:generation_interval)]
  shifted_reported_cases <- shifted_reported_cases[date %in% reported_cases$date]
  reported_cases <- reported_cases[date %in% shifted_reported_cases$date]


# Add week day info -------------------------------------------------------

  reported_cases <- reported_cases[, day_of_week := lubridate::wday(date, week_start = 1)][,
                                     `:=`(wkd = ifelse(day_of_week >= 6, 1, 0),
                                          mon = ifelse(day_of_week == 1, 1, 0))]
# Define stan model parameters --------------------------------------------

  data <- list(
    day_of_week = reported_cases$day_of_week,
    cases = reported_cases$confirm,
    shifted_cases = shifted_reported_cases$confirm,
    t = length(reported_cases$date),
    max_rep = reporting_delay$max,
    max_inc = incubation_period$max,
    inc_mean_mean = incubation_period$mean,
    inc_mean_sd = incubation_period$mean_sd,
    inc_sd_mean = incubation_period$sd,
    inc_sd_sd = incubation_period$sd_sd,
    rep_mean_mean = reporting_delay$mean,
    rep_mean_sd = reporting_delay$mean_sd,
    rep_sd_mean = reporting_delay$sd,
    rep_sd_sd = reporting_delay$sd_sd
  )  
  
  ## Set model to poisson or negative binomial
  if (family %in% "poisson") {
    data$model_type <- 1
  }else if (family %in% "negbin"){
    data$model_type <- 2
  }
  

# Set up initial conditions fn --------------------------------------------

init_fun <- function(){list(noise = rnorm(data$t, 1, 0.1),
                            day_of_week_eff = as.vector(MCMCpack::rdirichlet(1, rep(1, 7))),
                            phi = rexp(1, 1),
                            scale_inf = rep(1, data$t))}
  
# Load and run the stan model ---------------------------------------------

  if (missing(model)) {
    model <- rstan::stan_model("nowcast.stan")
  }

  
  if (verbose) {
    message(paste0("Running for ",samples + warmup," samples and ", data$t," time steps"))
  }
  

  fit <- rstan::sampling(model,
                           data = data,
                           chains = chains,
                           init = init_fun,
                           iter = samples + warmup, 
                           warmup = warmup,
                           cores = cores,
                           # control = list(adapt_delta = 0.99,
                           #                max_treedepth = 12),
                           refresh = ifelse(verbose, 50, 0))
    
    
    
    # Extract parameters of interest from the fit -----------------------------
    
    ## Extract sample from stan object
    samples <- rstan::extract(fit)
    
  
  return(list(fit = fit, samples = samples))
}

