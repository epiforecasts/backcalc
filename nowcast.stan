functions {
  // convolve a pdf and case vector using matrix multiplication
  vector convolve(vector cases, vector pdf, int direction) {
    int t = num_elements(cases);
    matrix[t, t] delay_mat = rep_matrix(0, t, t);
    int max_pdf = num_elements(pdf);
    row_vector[max_pdf] row_pdf = to_row_vector(pdf);
    vector[t] convolved_cases;

    for (s in 1:t) {
      if (direction) {
        int max_length = min(s, max_pdf);
        delay_mat[s, (s - max_length + 1):s] = row_pdf[(max_pdf - max_length + 1):max_pdf];
      }else{
        int max_length = min(t - s, max_pdf - 1);
        delay_mat[s, s:(s + max_length)] = row_pdf[1:(max_length + 1)];
      }
    }
    
   convolved_cases = delay_mat * to_vector(cases);

   return(convolved_cases);
  }

  // apply backsampling and upscaling based on a pdf
  vector backsample(vector cases, vector pdf) {
    int t = num_elements(cases);
    vector[t] backsampled_cases;
    int max_upscale = min(t, num_elements(pdf));
    int pdf_length = num_elements(pdf);
    vector[pdf_length] cdf;
    
    backsampled_cases = convolve(cases, pdf, 0);
    
    // apply upscaling
    cdf = cumulative_sum(pdf);
    
    for (i in  1:max_upscale) {
      backsampled_cases[(t - i + 1)] = (backsampled_cases[(t - i + 1)] + 1) / cdf[i];
    }
    
    //bound last day to equal day before
    backsampled_cases[t] = backsampled_cases[t -1];
    
    return(backsampled_cases);
  }
  
  // discretised lognormal pmf
  real discretised_lognormal_pmf(int y, real mu, real sigma) {
    return((lognormal_cdf(y, mu, sigma) - lognormal_cdf(y - 1, mu, sigma)));
  }
  
  // discretised gamma pmf
  real discretised_gamma_pmf(int y, real mu, real sigma) {
    // calculate alpha and beta for gamma distribution
    real alpha = (mu / sigma)^2;
    real beta = (sigma^2) / mu;
    return((gamma_cdf(y, alpha, beta) - gamma_cdf(y - 1, alpha, beta)));
  }
}


data {
  int t;                             // number of time steps
  int day_of_week[t];                // day of the week indicator (1 - 7)
  int <lower = 0> cases[t];          // observed cases
  vector<lower = 0>[t] shifted_cases;// median shifted smoothed cases
  real inc_mean_sd;                  // prior sd of mean incubation period
  real inc_mean_mean;                // prior mean of mean incubation period
  real inc_sd_mean;                  // prior sd of sd of incubation period
  real inc_sd_sd;                    // prior sd of sd of incubation period
  int max_inc;                       // maximum incubation period
  real rep_mean_mean;                // prior mean of mean reporting delay
  real rep_mean_sd;                  // prior sd of mean reporting delay
  real rep_sd_mean;                  // prior mean of sd of reporting delay
  real rep_sd_sd;                    // prior sd of sd of reporting delay
  int max_rep;                       // maximum report delay
  real <lower = 0> r_mean;           // prior mean of reproduction number
  real <lower = 0> r_sd;             // prior standard deviation of reproduction number
  real gt_mean_sd;                   // prior sd of mean generation time
  real gt_mean_mean;                 // prior mean of mean generation time
  real gt_sd_mean;                   // prior sd of sd of generation time
  real gt_sd_sd;                     // prior sd of sd of generation time
  int max_gt;                        // maximum generation time
  int model_type;                    // type of model: 1 = poisson otherwise negative binomial
  int estimate_r;                    // should the reproduction no be estimated (1 = yes)
}

transformed data{
  real r_alpha;                      // alpha parameter of the R gamma prior
  real r_beta;                       // beta parameter of the R gamma prior
  real initial_noise_mean;           // initial noise estimate
  real initial_R_mean;               // initial R estimate

  // calculate alpha and beta for gamma distribution
  r_alpha = (r_mean / r_sd)^2;
  r_beta = (r_sd^2) / r_mean;
  
  // precalculate parameters for noise
  initial_noise_mean = log(1);
  
  // precalculate parameters for R
  initial_R_mean = log(1);
}
parameters{
  real<lower = 0> initial_noise;                      // initial noise on the mean shifted observed cases
  vector<lower = 0>[t -1] noise_diff;                 // temporal scaling on noise
  simplex[7] day_of_week_eff_raw;                     // day of week reporting effect + control parameters
  real <lower = 0> inc_mean;                          // mean of incubation period
  real <lower = 0> inc_sd;                            // sd of incubation period
  real <lower = 0> rep_mean;                          // mean of reporting delay
  real <lower = 0> rep_sd;                            // sd of incubation period
  real<lower = 0> rep_phi;                            // overdispersion of the reporting process
  vector<lower = 0>[estimate_r] initial_R;            // initial reproduction number
  vector<lower = 0>[estimate_r > 0 ? t -1 : 0] R_diff;// scaling of reproduction number over time
  real<lower = 0> gt_mean[estimate_r];                // mean of generation time
  real <lower = 0> gt_sd[estimate_r];                 // sd of generation time
  real<lower = 0> inf_phi[estimate_r];                // overdispersion of the infection process
}

transformed parameters {
  // stored transformed parameters
  vector<lower = 0>[t] noise;                             // noise on the mean shifted observed cases
  vector<lower = 0>[t] infections;                        // infections over time
  vector<lower = 0>[t] reports;                           // reports over time
  vector[7] day_of_week_eff;                              // day of the week effect
  vector<lower = 0>[estimate_r > 0 ? t : 0] R;            // reproduction number over time
  vector[estimate_r > 0 ? t : 0] branch_reports;          // reports generated by the branching process
 {
  // temporary transformed parameters
  vector[max_rep] rev_delay;                              // reversed report delay pdf
  vector[max_inc] rev_incubation;                         // reversed incubation period pdf
  vector[t] onsets;                                       // onsets over time
  vector[estimate_r > 0 ? max_gt : 0] rev_generation_time;// reversed generation time pdf
  vector[estimate_r > 0 ? t : 0] infectiousness;          // infections over time
  vector[estimate_r > 0 ? t : 0] branch_onsets;           // onsets generated by the branching process
  vector[estimate_r > 0 ? t : 0] branch_infections;       // infections generated by the branching process
  
  // reverse the distributions to allow vectorised access
  for (j in 1:max_inc) {
    rev_incubation[j] =
        discretised_lognormal_pmf(max_inc - j + 1, inc_mean, inc_sd);
  }
  
  for (j in 1:max_rep) {
    rev_delay[j] =
        discretised_lognormal_pmf(max_rep - j + 1, rep_mean, rep_sd);
  }
    
  // define day of the week effect
  day_of_week_eff = 7 * day_of_week_eff_raw;

  // walk over noise terms
  noise[1] = initial_noise;
  for (s in 2:t) {
    noise[s] = noise[s - 1] * noise_diff[s - 1];
  }
  
  // generate infections from backcalculated and non-parameteric noise
  for (s in 1:t) {
    infections[s] = shifted_cases[s] * noise[s];
  }
  
  // onsets from infections
  onsets = convolve(infections, rev_incubation, 1);

  // reports from onsets
  reports = convolve(onsets, rev_delay, 1);

  for (s in 1:t) {
    // add reporting effects (adjust for simplex scale)
    reports[s] *= day_of_week_eff[day_of_week[s]];
  }
    
  ////////////////////////////////////////////////////// 
  // estimate reproduction no from a branching process
  //////////////////////////////////////////////////////
  if (estimate_r) {
    // calculate pdf of generation time from distribution
    for (j in 1:(max_gt)) {
       rev_generation_time[j] =
           discretised_gamma_pmf(max_gt - j + 1, gt_mean[estimate_r], gt_sd[estimate_r]);
     }
     
     // infectiousness from infections
     infectiousness = convolve(infections, rev_generation_time, 1);
  
     // reproduction number is a scaling of previous timestep
     R[1] = initial_R[estimate_r];
     for (s in 2:t) {
       R[s] = R[s - 1] * R_diff[s - 1];
     }
     
     // Estimate infections using branching process
     branch_infections = R .* infectiousness;
     
     // onsets from infections
     branch_onsets = convolve(branch_infections, rev_incubation, 1);

     // reports from onsets
     branch_reports = convolve(branch_onsets, rev_delay, 1);

     for (s in 1:t) {
      // add reporting effects (adjust for simplex scale)
      branch_reports[s] *= day_of_week_eff[day_of_week[s]];
    }
  }
 }
}

model {
  // reporting overdispersion
  rep_phi ~ exponential(1);

  // noise on median shift
  initial_noise ~ lognormal(initial_noise_mean, 0.1);
  noise_diff ~ lognormal(initial_noise_mean, 0.1);

  // daily cases given reports
  if (model_type == 1) {
    target += poisson_lpmf(cases | reports);
  }else{
    target += neg_binomial_2_lpmf(cases | reports, rep_phi);
  }

  // penalised priors for incubation period, and report delay
  target += normal_lpdf(inc_mean | inc_mean_mean, inc_mean_sd) * t;
  target += normal_lpdf(inc_sd | inc_sd_mean, inc_sd_sd) * t;
  target += normal_lpdf(rep_mean | rep_mean_mean, rep_mean_sd) * t;
  target += normal_lpdf(rep_sd | rep_sd_mean, rep_sd_sd) * t;
  
  ////////////////////////////////////////////////////// 
  // estimate reproduction no from a branching process
  //////////////////////////////////////////////////////
  if (estimate_r) {
    
    // infection overdispersion
    inf_phi ~ exponential(1);

    // initial prior on R
    initial_R ~ gamma(r_alpha, r_beta);
    R_diff ~ lognormal(initial_R_mean, 1);
     
    // penalised_prior on generation interval
    target += normal_lpdf(gt_mean | gt_mean_mean, gt_mean_sd) * t;
    target += normal_lpdf(gt_sd | gt_sd_mean, gt_sd_sd) * t;
  
    // Likelihood of Rt given infections
    if (model_type == 1) {
      target += poisson_lpmf(cases | branch_reports);
    }else{
      target += neg_binomial_2_lpmf(cases | branch_reports, inf_phi[estimate_r]);
    }
  }
}
  
generated quantities {
  int imputed_infections[t];
  
  // simulated infections - assume poisson (with negative binomial reporting)
  imputed_infections = poisson_rng(infections);
}