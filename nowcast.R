#' @export
#' @importFrom rstan sampling extract
#' @importFrom data.table data.table copy merge.data.table setorder rbindlist setDTthreads melt .N
#' @importFrom purrr tranpose
#' @importFrom future.apply future_lapply
#' @importFrom lubridate wday
#' @importFrom truncnorm rtruncnorm
#' @examples
#' reported_cases <- NCoVUtils::get_ecdc_cases(countries = "Russia")
#' reported_cases <- NCoVUtils::format_ecdc_data(reported_cases)
#' reported_cases <- data.table::as.data.table(reported_cases)[, confirm := cases][, cases := NULL]
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
#'                estimate_rt = TRUE,
#'                verbose = TRUE, return_all = TRUE
#'                )
#'
#' out                                   
nowcast <- function(reported_cases, family = "poisson",
                    incubation_period, reporting_delay,
                    generation_time, rt_prior,
                    prior_smoothing_window = 7,
                    model, cores = 1, chains = 2,
                    samples = 1000, warmup = 1000,
                    estimate_rt = FALSE, adapt_delta = 0.95,
                    max_treedepth = 15, return_all = FALSE,
                    verbose = FALSE){
  
  suppressMessages(data.table::setDTthreads(threads = 1))

  # Add prior for R if missing ---------------------------------

  if (missing(rt_prior)) {
    rt_prior <- list(mean = 1, sd = 2)
  }
  
  # Make sure there are no missing dates and order cases --------------------
  reported_cases_grid <- data.table::copy(reported_cases)[, .(date = seq(min(date), max(date), by = "days"))]

  nowcast <- data.table::merge.data.table(
    reported_cases , reported_cases_grid, 
    by = c("date"), all.y = TRUE)
  
  reported_cases <- reported_cases[is.na(confirm), confirm := 0 ][,.(date = date, confirm)]
  reported_cases <- data.table::setorder(reported_cases, date)
  
  ## Filter out 0 reported cases from the end of the data
  reported_cases <- reported_cases[order(-date)][, 
                                cum_cases := cumsum(confirm)][cum_cases > 0][, cum_cases := NULL]
  ## Filter out 0 reported cases from the beginning of the data
  reported_cases <- reported_cases[order(date)][,
                                cum_cases := cumsum(confirm)][cum_cases > 0][, cum_cases := NULL]
  
# Estimate the mean delay -----------------------------------------------
  
  mean_shift <- exp(incubation_period$mean) + exp(reporting_delay$mean)

# Add the mean delay and incubation period on as 0 case days ------------

  reported_cases <- data.table::rbindlist(list(
    data.table::data.table(date = seq(min(reported_cases$date) - mean_shift - prior_smoothing_window, 
                                      min(reported_cases$date) - 1, by = "days"),
                           confirm = 0),
    reported_cases
  ))  

# Calculate smoothed prior cases ------------------------------------------
  
  shifted_reported_cases <- data.table::copy(reported_cases)[,
                          confirm := data.table::shift(confirm, n = as.integer(mean_shift),
                          type = "lead", fill = data.table::last(confirm))][,
                          confirm := data.table::frollmean(confirm, n = prior_smoothing_window, 
                                                           align = "right", fill = data.table::last(confirm))][,
                          confirm := data.table::fifelse(confirm == 0, 1e-4, confirm)]
  
  ##Drop median generation interval initial values
  shifted_reported_cases <- shifted_reported_cases[-(1:prior_smoothing_window)]
  reported_cases <- reported_cases[-(1:prior_smoothing_window)]
  
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
    r_sd = rt_prior$sd,
    estimate_r = ifelse(estimate_rt, 1, 0)
  )  
  
  ## Set model to poisson or negative binomial
  if (family %in% "poisson") {
    data$model_type <- 0
  }else if (family %in% "negbin"){
    data$model_type <- 1
  }
  
# Set up initial conditions fn --------------------------------------------

init_fun <- function(){out <- list(
                            initial_noise = truncnorm::rtruncnorm(1, a = 0, mean = 1, sd = 0.4),
                            noise_diff = truncnorm::rtruncnorm(data$t - 1, a = 0, mean = 1, sd = 0.1),
                            inc_mean = truncnorm::rtruncnorm(1, a = 0, mean = incubation_period$mean, sd = incubation_period$mean_sd),
                            inc_sd = truncnorm::rtruncnorm(1, a = 0, mean = incubation_period$sd, sd = incubation_period$sd_sd),
                            rep_mean = truncnorm::rtruncnorm(1, a = 0, mean = reporting_delay$mean, sd = reporting_delay$mean_sd),
                            rep_sd = truncnorm::rtruncnorm(1, a = 0, mean = reporting_delay$sd,  sd = reporting_delay$sd_sd))
                        
                        if (data$model_type == 1) {
                          out$rep_phi <- array(rexp(1, 1))
                        }

                        if (estimate_rt) {
                        out$R <- rgamma(n = data$t, shape = (rt_prior$mean / rt_prior$sd)^2, 
                                                    scale = (rt_prior$sd^2) / rt_prior$mean)
                        out$gt_mean <- array(truncnorm::rtruncnorm(1, a = 0, mean = generation_time$mean,  
                                                             sd = generation_time$mean_sd))
                        out$gt_sd <-  array(truncnorm::rtruncnorm(1, a = 0, mean = generation_time$sd,
                                                            sd = generation_time$sd_sd))
                        }

                return(out)
}
  
# Load and run the stan model ---------------------------------------------

  if (missing(model)) {
    model <- rstan::stan_model("nowcast.stan")
  }
  
  if (verbose) {
    message(paste0("Running for ",samples + warmup," samples and ", data$t," time steps"))
  }

  fit <-
    rstan::sampling(model,
                    data = data,
                    chains = chains,
                    init = init_fun,
                    iter = samples + warmup, 
                    warmup = warmup,
                    cores = cores,
                    control = list(adapt_delta = adapt_delta,
                                   max_treedepth = max_treedepth),
                    refresh = ifelse(verbose, 50, 0))

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
    
    if (estimate_rt) {
      out$R <- extract_parameter("R", 
                                 samples,
                                 reported_cases$date)
    }

    
    if (return_all) {
      ## Add prior infections
      out$prior_infections <- shifted_reported_cases[, .(parameter = "prior_infections", time = 1:.N, 
                                                         date, value = confirm, sample = 1)]
        
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
      
      if (estimate_rt) {
        out$gt_mean <- extract_static_parameter("gt_mean")
        out$gt_mean[, value := value.V1][, value.V1 := NULL]
        
        out$gt_sd <- extract_static_parameter("gt_sd")
        out$gt_sd[, value := value.V1][, value.V1 := NULL]
      }

      out$fit <- fit
  }
  
  
  return(out)
}