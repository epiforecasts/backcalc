
generate_simulations <- function(rts, initial_cases, 
                     generation_time, incubation_period,
                     reporting_delay, reporting_effect,
                     samples = 10) {
  

  ## Generation time
  generation_defs <- EpiNow::gamma_dist_def(mean = generation_time$mean,
                                            mean_sd = generation_time$mean_sd,
                                            sd = generation_time$sd,
                                            sd_sd = generation_time$sd_sd,
                                            max_value = generation_time$max, samples = samples)
  
  generate_pdf <- function(dist, max_value) {
    ## Define with 0 day padding
    sample_fn <- function(n, ...) {
      c(0, EpiNow::dist_skel(n = n, 
                        model = dist$model[[1]], 
                        params = dist$params[[1]],
                        max_value = dist$max_value[[1]] - 1, 
                        ...))
    }
    
    dist_pdf <- sample_fn(0:max_value, dist = TRUE, cum = FALSE)
    
    return(dist_pdf)
  }
  
  generation_pdfs <- purrr::map(split(generation_defs[, index := 1:.N], by = "index"),
                           generate_pdf, max_value = generation_time$max)
  
  ## Sample a report delay as a lognormal
  delay_defs <- EpiNow::lognorm_dist_def(mean = reporting_delay$mean, mean_sd = reporting_delay$mean_sd,
                                         sd = reporting_delay$sd, sd_sd = reporting_delay$sd_sd,
                                         max_value = reporting_delay$max, samples = samples)
  
  
  ## Sample a incubation period (again using the default for covid) 
  incubation_defs <- EpiNow::lognorm_dist_def(mean = incubation_period$mean,
                                              mean_sd = incubation_period$mean_sd,
                                              sd = incubation_period$sd,
                                              sd_sd = incubation_period$sd_sd,
                                              max_value = incubation_period$max,
                                              samples = samples)
  
  
  
  ## Run for each sample im turn
  dist_defs <- list(delay = split(delay_defs[, index := 1:.N], by = "index"),
                    incubation = split(incubation_defs[, index := 1:.N], by = "index"),
                    generation_time = generation_pdfs)
  
  dist_defs <- purrr::transpose(dist_defs)
  
  safe_sim <- purrr::safely(EpiNow::simulate_cases)
  
  out <- purrr::map(dist_defs,
                    function(dist) {
                      safe_sim(rts, initial_cases, initial_date = as.Date("2020-03-01"),
                               generation_interval = dist$generation_time,
                               delay_def = dist$delay,
                               incubation_def = dist$incubation,
                               reporting_effect = reporting_effect)[[1]]
                    })
  
  out <- data.table::rbindlist(out, idcol = "sample")                   

  return(out)
}


