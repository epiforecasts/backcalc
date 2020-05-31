#' @export
#' @importFrom rstan sampling extract
#' @importFrom data.table data.table copy merge.data.table setorder rbindlist setDTthreads melt .N
#' @importFrom purrr tranpose
#' @importFrom future.apply future_lapply
#' @importFrom lubridate wday
#' @examples
#' reported_cases <- NCoVUtils::get_ecdc_cases(countries = "Austria")
#' reported_cases <- NCoVUtils::format_ecdc_data(reported_cases)
#' reported_cases <- data.table::as.data.table(reported_cases)[, confirm := cases][, cases := NULL]
#'   
#' ## Sample a report delay as a lognormal
#' delay_defs <- EpiNow::lognorm_dist_def(mean = 5, mean_sd = 1,
#'                                        sd = 2, sd_sd = 1, max_value = 30,
#'                                        samples = 10, to_log = TRUE)
#'                                       
#' 
#' ## Sample a incubation period (again using the default for covid)
#' incubation_defs <- EpiNow::lognorm_dist_def(mean = EpiNow::covid_incubation_period[1, ]$mean,
#'                                           mean_sd = EpiNow::covid_incubation_period[1, ]$mean_sd,
#'                                           sd = EpiNow::covid_incubation_period[1, ]$sd,
#'                                           sd_sd = EpiNow::covid_incubation_period[1, ]$sd_sd,
#'                                           max_value = 30, samples = 10)
#'
#'  generation_interval <- rowMeans(EpiNow::covid_generation_times)
#'  generation_interval <- sum(!(cumsum(generation_interval) > 0.5)) + 1   
#'                                      
#' out <- nowcast(reported_cases, family = "poisson",
#'                delay_defs = delay_defs, incubation_defs = incubation_defs,
#'                generation_interval = generation_interval,
#'                include_fit = TRUE, verbose = TRUE) 
#'                
#' out                                   
nowcast <- function(reported_cases, family = "poisson",
                    delay_defs,
                    incubation_defs,
                    generation_interval,
                    cores = 1,
                    chains = 2,
                    samples = 1000,
                    warmup = 1000,
                    batch = TRUE,
                    include_fit = FALSE,
                    verbose = FALSE) {
  
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

  # Balance input delay and incubation samples ------------------------------
  
  balance_dfs <- function(df1, df2) {
    if (nrow(df1) > nrow(df2)) {
      df2 <- data.table::rbindlist(list(
        df2,
        df2[sample(1:nrow(df2), (nrow(df1) - nrow(df2)), replace = TRUE), ]
      ))
    }
    return(df2)
  }
  
  incubation_defs <- balance_dfs(delay_defs, incubation_defs)
  delay_defs <- balance_dfs(incubation_defs, delay_defs)
  
  # Calculate CDFs of distributions -----------------------------------------

  generate_cdf <- function(dist, max_value) {
    ## Define sample delay fn
    sample_fn <- function(n, ...) {
      EpiNow::dist_skel(n = n, 
                        model = dist$model[[1]], 
                        params = dist$params[[1]],
                        max_value = dist$max_value[[1]], 
                        ...)
    }
    
    dist_cdf <- sample_fn(0:max_value, dist = TRUE, cum = FALSE)
    
    return(dist_cdf)
  }
  
  delay_cdfs <- purrr::map(split(delay_defs[, index := 1:.N], by = "index"),
                           generate_cdf, max_value = nrow(reported_cases))
  
  delay_cdfs <- do.call(rbind, delay_cdfs)
  
  incubation_cdfs <- purrr::map(split(incubation_defs[, index := 1:.N], by = "index"),
                           generate_cdf, max_value = nrow(reported_cases))
  
  incubation_cdfs <- do.call(rbind, incubation_cdfs)
    

# Estimate the median delay -----------------------------------------------
  median_delay <- sum(!(cumsum(colMeans(delay_cdfs)) > 0.5)) + 1
  median_incubation <- sum(!(cumsum(colMeans(incubation_cdfs)) > 0.5)) + 1
  
  median_shift <- median_delay + median_incubation 

# Add the median delay and incubation period on as 0 case days ------------

  reported_cases <- data.table::rbindlist(list(
    data.table::data.table(date = seq(min(reported_cases$date) - median_shift - generation_interval, 
                                      min(reported_cases$date) - 1, by = "days"),
                           confirm = 0),
    reported_cases
  ))  
  
  # Calculate smoothed prior cases ------------------------------------------

  shifted_reported_cases <- data.table::copy(reported_cases)[,
                      confirm := data.table::shift(confirm, n = as.integer(median_shift),
                                                   type = "lead", fill = data.table::last(confirm))][,
                      confirm := data.table::frollmean(confirm, n = generation_interval, 
                                                       align = "center", fill = data.table::last(confirm))]
  
  ##Drop median generation interval initial values
  shifted_reported_cases <- shifted_reported_cases[-(1:generation_interval)]
  reported_cases <- reported_cases[-(1:generation_interval)]


# Add week day info -------------------------------------------------------

  reported_cases <- reported_cases[, weekday := lubridate::wday(date)][,
                                     `:=`(wkd = ifelse(weekday >= 6, 1, 0),
                                          mon = ifelse(weekday == 1, 1, 0))]
# Define stan model parameters --------------------------------------------

  data <- list(
    wkd = reported_cases$wkd, 
    mon = reported_cases$mon,
    cases = reported_cases$confirm,
    shifted_cases = shifted_reported_cases$confirm,
    delay = delay_cdfs,
    incubation = incubation_cdfs,
    t = length(reported_cases$date),
    d = ncol(delay_cdfs),
    inc = ncol(incubation_cdfs),
    samples = nrow(delay_cdfs)
  )  
  
  ## Set model to poisson or negative binomial
  if (family %in% "poisson") {
    data$model_type <- 1
  }else if (family %in% "negbin"){
    data$model_type <- 2
  }
  

# Set up initial conditions fn --------------------------------------------

init_fun <- function(){list(noise = rnorm(data$t, 1, 0.2),
                            wkd_eff = rnorm(1, 0, 0.1),
                            mon_eff = rnorm(1, 0, 0.1),
                            phi = rexp(1, 1))}
  
# Load and run the stan model ---------------------------------------------

  model <- rstan::stan_model("nowcast.stan")
  
  if (verbose) {
    message(paste0("Running for ",data$samples," samples and ", data$t," time steps"))
  }
  
  run_model <- function(data) { 
    fit <- suppressWarnings(
           rstan::sampling(model,
                           data = data,
                           chains = chains,
                           init = init_fun,
                           iter = round(samples / (chains *  nrow(delay_cdfs))) + warmup, 
                           warmup = warmup,
                           cores = cores,
                           refresh = ifelse(verbose, 50, 0)))
    
    
    
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
      
      param_df <- param_df[, var := NULL][,
                                          sample := 1:.N, by = .(time)]
      param_df <- param_df[, date := dates, by = .(sample)]
      param_df <- param_df[, .(parameter = param, time, date, 
                               sample, value)]
      
      return(param_df)
    }
    
    ## Report infections, prior infections and noise
    out$infections <- extract_parameter("imputed_infections", 
                                        samples,
                                        reported_cases$date)
    
    out$noise <- extract_parameter("noise", 
                                   samples,
                                   reported_cases$date)
    
    out$wkd_eff <- data.table::data.table(
      parameter = "wkd_eff",
      sample = 1:length(samples$wkd_eff),
      value = samples$wkd_eff)
    
    
    out$mon_eff <- data.table::data.table(
      parameter = "mon_eff",
      sample = 1:length(samples$mon_eff),
      value = samples$mon_eff)
    
    if (include_fit) {
      out$fit <- fit
    }

    return(out)
  }
  
  
# Run either all at once or in batch --------------------------------------

  if(data$samples == 1) {
    batch <- FALSE
  }
  
  if (!batch) {
    out <- run_model(data)
  }else{
    out <- future.apply::future_lapply(1:data$samples,
                                       function(i, stan_data) {
                                         stan_data$samples <- 1
                                         stan_data$delay <- t(as.matrix(stan_data$delay[i, ]))
                                         stan_data$incubation <- t(as.matrix(stan_data$incubation[i,]))
                                         
                                         out <- run_model(stan_data)
                                         
                                         return(out)
                                       }, stan_data = data)
    
    out <- purrr::transpose(out)
    
    bind_out <- function(out_list, by_var) {
      out_list <- data.table::rbindlist(out_list)
      out_list <- out_list[, sample := 1:.N, by = by_var]
      
      return(out_list)
    }
    
    
    out$infections <- bind_out(out$infections, by_var = "time")
    out$prior_infections <- out$prior_infections[[1]]
    out$noise <- bind_out(out$noise, by_var = "time")
    out$wkd_eff <- bind_out(out$wkd_eff, by_var = "parameter")
    out$mon_eff <- bind_out(out$mon_eff, by_var = "parameter")
  }
  
  ## Add prior infections
  out$prior_infections <- shifted_reported_cases[, .(parameter = "prior_infections", time = 1:.N, 
                                                     date, value = confirm)]
  
  
  return(out)
}