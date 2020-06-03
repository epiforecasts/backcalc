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
#'  
#'  
#' generation_time <- list(mean = EpiNow::covid_generation_times_summary[1, ]$mean,
#'                         mean_sd = EpiNow::covid_generation_times_summary[1, ]$mean_sd,
#'                         sd = EpiNow::covid_generation_times_summary[1, ]$sd,
#'                         sd_sd = EpiNow::covid_generation_times_summary[1, ]$sd_sd,
#'                         max = 30)
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
#' rt_prior <- list(mean = 2.6, sd = 2)
#'   
#' ## Compile model
#' model <- rstan::stan_model("nowcast.stan")
#' 
#' ## Run model
#' out <- nowcast(reported_cases, family = "poisson",
#'                generation_time = generation_time,
#'                incubation_period = incubation_period,
#'                reporting_delay = reporting_delay,
#'                rt_prior = rt_prior,
#'                model = model,
#'                cores = 4, chains = 4,
#'                verbose = TRUE, return_all = TRUE
#'                )
#'
#' out                                   
nowcast <- function(reported_cases, family = "poisson",
                    incubation_period, reporting_delay,
                    generation_time, rt_prior, model,
                    cores = 1,
                    chains = 2,
                    samples = 1000,
                    warmup = 1000,
                    return_all = FALSE,
                    verbose = FALSE){
  
  suppressMessages(data.table::setDTthreads(threads = 1))
  
  # Make sure there are no missing dates and order cases --------------------
  reported_cases_grid <- data.table::copy(reported_cases)[, .(date = seq(min(date), max(date), by = "days"))]

  nowcast <- data.table::merge.data.table(
    reported_cases , reported_cases_grid, 
    by = c("date"), all.y = TRUE)
  
  reported_cases <- reported_cases[is.na(confirm), confirm := 0 ][,.(date = date, confirm)]
  reported_cases <- data.table::setorder(reported_cases, date)
  
  ## Filter out 0 reported cases
  reported_cases <- reported_cases[, cum_cases := cumsum(confirm)][cum_cases > 0][, cum_cases := NULL]

# Estimate the mean delay -----------------------------------------------
  
  mean_shift <- incubation_period$mean + reporting_delay$mean

# Integer mean generation interval ----------------------------------------

  mean_int_gt <- round(generation_time$mean)

# Add the mean delay and incubation period on as 0 case days ------------

  reported_cases <- data.table::rbindlist(list(
    data.table::data.table(date = seq(min(reported_cases$date) - mean_shift - mean_int_gt, 
                                      min(reported_cases$date) - 1, by = "days"),
                           confirm = 0),
    reported_cases
  ))  
  
  # Calculate smoothed prior cases ------------------------------------------

  shifted_reported_cases <- data.table::copy(reported_cases)[,
                      confirm := data.table::shift(confirm, n = as.integer(mean_shift),
                                                   type = "lead", fill = data.table::last(confirm))][,
                      confirm := data.table::frollmean(confirm, n = mean_int_gt, 
                                                       align = "center", fill = data.table::last(confirm))][,
                      confirm := data.table::fifelse(confirm == 0, 1e-4, confirm)]
  
  ##Drop median generation interval initial values
  shifted_reported_cases <- shifted_reported_cases[-(1:mean_int_gt)]
  reported_cases <- reported_cases[-(1:mean_int_gt)]


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
    inc_mean_mean = incubation_period$mean,
    inc_mean_sd = incubation_period$mean_sd,
    inc_sd_mean = incubation_period$sd,
    inc_sd_sd = incubation_period$sd_sd,
    max_inc = incubation_period$max,
    rep_mean_mean = reporting_delay$mean,
    rep_mean_sd = reporting_delay$mean_sd,
    rep_sd_mean = reporting_delay$sd,
    rep_sd_sd = reporting_delay$sd_sd,
    max_rep = reporting_delay$max,
    gt_mean_mean = generation_time$mean,
    gt_mean_sd = generation_time$mean_sd,
    gt_sd_mean = generation_time$sd,
    gt_sd_sd = generation_time$sd_sd,
    max_gt = generation_time$max,
    r_mean = rt_prior$mean,
    r_sd = rt_prior$sd
  )  
  
  ## Set model to poisson or negative binomial
  if (family %in% "poisson") {
    data$model_type <- 1
  }else if (family %in% "negbin"){
    data$model_type <- 2
  }
  

# Set up initial conditions fn --------------------------------------------

init_fun <- function(){list(noise = rnorm(data$t, 1, 0.1),
                            day_of_week_eff= rnorm(7, 1, 0.1),
                            rep_phi = rexp(1, 1),
                            R = rgamma(n = data$t, 
                                       shape = (rt_prior$mean / rt_prior$sd)^2, 
                                       scale = (rt_prior$sd^2) / rt_prior$mean),
                            inf_phi = rexp(1, 1))}
  
# Load and run the stan model ---------------------------------------------

  if (missing(model)) {
    model <- rstan::stan_model("nowcast.stan")
  }

  
  if (verbose) {
    message(paste0("Running for ",samples + warmup," samples and ", data$t," time steps"))
  }
  

  fit <- suppressWarnings(
           rstan::sampling(model,
                           data = data,
                           chains = chains,
                           init = init_fun,
                           iter = samples + warmup, 
                           warmup = warmup,
                           cores = cores,
                           refresh = ifelse(verbose, 50, 0))
           )
    
    
    
    # Extract parameters of interest from the fit -----------------------------
    
    ## Extract sample from stan object
    samples <- rstan::extract(fit)
    
    ## Construct reporting list
    out <- list()
    
    ## Generic data.frame reporting function
    extract_parameter <- function(param, samples, dates) {
      param_df <- data.table::as.data.table(
        t(
          data.table::as.data.table(
            samples[[param]]
          )
        ))
      
      param_df <- param_df[, time := 1:.N]
      param_df <- 
        data.table::melt(param_df, id.vars = "time",
                         variable.name = "var")
      
      param_df <- param_df[, var := NULL][, sample := 1:.N, by = .(time)]
      param_df <- param_df[, date := dates, by = .(sample)]
      param_df <- param_df[, .(parameter = param, time, date, 
                               sample, value)]
      
      return(param_df)
    }
    
    ## Report infections, and R
    out$infections <- extract_parameter("imputed_infections", 
                                        samples,
                                        reported_cases$date)
    
    out$R <- extract_parameter("R", 
                                samples,
                                reported_cases$date)
    
    if (return_all) {
      out$noise <- extract_parameter("noise", 
                                     samples,
                                     reported_cases$date)
      
      out$day_of_week <- extract_parameter("day_of_week_eff", 
                                           samples,
                                           1:7)
      
      char_day_of_week <- data.table::data.table(wday = c("Monday", "Tuesday", "Wednesday",
                                                          "Thursday", "Friday", "Saturday",
                                                          "Sunday"),
                                                 time = 1:7)
      out$day_of_week <- out$day_of_week[char_day_of_week, on = "time"][, 
                                         wday_numeric := time][, 
                                         time := NULL]

      extract_static_parameter <- function(param) {
        data.table::data.table(
          parameter = param,
          sample = 1:length(samples[[param]]),
          value = samples[[param]])
      }
      
      out$inc_mean <- extract_static_parameter("inc_mean")
      
      out$inc_sd <- extract_static_parameter("inc_sd")
      
      out$rep_mean <- extract_static_parameter("rep_mean")
      
      out$rep_sd <- extract_static_parameter("rep_sd")
      
      out$gt_mean <- extract_static_parameter("gt_mean")
      
      out$gt_sd <- extract_static_parameter("gt_sd")
      
      out$infection_overdispersion <- extract_static_parameter("inf_phi")
      
      out$fit <- fit
      
      ## Add prior infections
      out$prior_infections <- shifted_reported_cases[, .(parameter = "prior_infections", time = 1:.N, 
                                                       date, value = confirm)]
  }
  
  
  return(out)
}

