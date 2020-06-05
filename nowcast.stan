functions {
  // convolve a pdf and case vector using matrix multiplication
  vector convolve(vector cases, vector pdf, int direction) {
    int t = num_elements(cases);
    matrix[t, t] delay_mat = rep_matrix(0, t, t);
    int max_pdf = num_elements(pdf) + 1;
    row_vector[max_pdf] row_pdf = direction ? to_row_vector(append_row(pdf, 0.0)) :
                                              to_row_vector(append_row(0.0, pdf)) ;
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
  
   // initialise all entries as non-zero 
   if (direction) {
    convolved_cases[1] = 0.00001;
   }else{
    convolved_cases[t] =  0.00001; 
   }
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
    
    print("Pre-upscale")
    print(backsampled_cases)
    
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
  vector<lower = 0>[t] shifted_cases;// backsampled cases
  vector[max_rep] delay;             // report delay pdf
  vector[max_inc] incubation;        // incubation period pdf
  int<lower = 0> weekly_cases[t];    // weekly observed cases
  int<lower = 0> cum_cases[t];       // cumulative cases
  real r_alpha;                      // alpha parameter of the R gamma prior
  real r_beta;                       // beta parameter of the R gamma prior

  // calculate alpha and beta for gamma distribution
  r_alpha = (r_mean / r_sd)^2;
  r_beta = (r_sd^2) / r_mean;
  
  // calculate weekly cases
  cum_cases[1] = cases[1];
  weekly_cases[1] = cases[1];
  for (s in 2:t) { 
    cum_cases[s] = cum_cases[s - 1] + cases[s];
    weekly_cases[s] = cum_cases[s] - cum_cases[max(1, s - 7)];
  }
  
  // initialise cases using a backwards convolution from report
  // upscale to adjust for right truncation
  {
    for (j in 1:max_rep) {
      delay[j] =
        discretised_lognormal_pmf(j, rep_mean_mean, rep_sd_mean);
    }
    print("Cases")
    print(cases)
    shifted_cases = backsample(to_vector(cases) + 0.00001, delay);

    print("Onsets")
    print(shifted_cases)
  }
  // then from onset (also with upscaling)
  {
    for (j in 1:max_inc) {
      incubation[j] =
           discretised_lognormal_pmf(j, inc_mean_mean, inc_sd_mean);
     }
    print("Incubation")
    print(incubation)
    shifted_cases = backsample(shifted_cases, incubation);

    print("Infections")
    print(shifted_cases)
  }
}
parameters{
  vector<lower = 0>[t] noise;                      // noise on the mean shifted observed cases
  vector[7] day_of_week_eff_raw;                   // day of week reporting effect + control parameters
  vector[7] mu;
  cholesky_factor_corr[7] L_Sigma;
  real <lower = 0> inc_mean;                       // mean of incubation period
  real <lower = 0> inc_sd;                         // sd of incubation period
  real <lower = 0> rep_mean;                       // mean of reporting delay
  real <lower = 0> rep_sd;                         // sd of incubation period
  real<lower = 0> rep_phi;                         // overdispersion of the reporting process
  vector<lower = 0>[estimate_r > 0 ? t : 0] R;     // effective reproduction number over time
  real<lower = 0> gt_mean;                         // mean of generation time
  real <lower = 0> gt_sd;                          // sd of generation time
}

transformed parameters {
  vector[max_rep] rev_delay;                              // reversed report delay pdf
  vector[max_inc] rev_incubation;                         // reversed incubation period pdf
  vector<lower = 0>[t] infections;                        // infections over time
  vector<lower = 0>[t] onsets;                            // onsets over time
  vector<lower = 0>[t] reports;                           // reports over time
  vector<lower = 0>[t] cum_reports;                       // cumulative reported cases (unadjusted for reporting)
  vector<lower = 0>[t] weekly_reports;                    // weekly reported cases
  simplex[7] day_of_week_eff_unscaled;                    // unscaled day of week effect
  vector[7] day_of_week_eff;                              // day of the week effect
  vector[estimate_r > 0 ? max_gt : 0] rev_generation_time;// reversed generation time pdf
  vector[estimate_r > 0 ? t : 0] infectiousness;          // infections over time
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
  day_of_week_eff_unscaled = softmax(day_of_week_eff_raw);
  day_of_week_eff = 7 * day_of_week_eff_unscaled;

  // generate infections from backcalculated and non-parameteric noise (squared)
  for (s in 1:t) {
      infections[s] = shifted_cases[s] * noise[s];
  }
  
  // onsets from infections
  onsets = convolve(infections, rev_incubation, 1);

  // reports from onsets
  reports = convolve(onsets, rev_delay, 1);

  // calculate cumulative reports
  cum_reports = cumulative_sum(reports);
  
  for (s in 1:t) {
    // calculate weekly reports
    weekly_reports[s] = s == 1 ? cum_reports[1] : cum_reports[s] - cum_reports[max(1, s - 7)];
    
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
           discretised_gamma_pmf(max_gt - j + 1, gt_mean, gt_sd);
     }
     // infectiousness from infections
     infectiousness = convolve(infections, rev_generation_time, 1);
  
     // Estimate infections using branching process
     branch_infections = R .* infectiousness;
  }
}

model {
  //prior day of the week effect
  day_of_week_eff_raw ~ multi_normal_cholesky(mu, L_Sigma);
  mu ~ normal(0, 1);
  // this is uniform over all correlation matrices
  L_Sigma ~ lkj_corr_cholesky(1.0);
  
  // reporting overdispersion
  rep_phi ~ exponential(1);

  // noise on median shift
  for (i in 1:t) {
    noise[i] ~ normal(1, 0.4) T[0,];
  }

  // daily cases given reports
  if (model_type == 1) {
    target += poisson_lpmf(cases | reports);
    target += poisson_lpmf(weekly_cases[7:t] | weekly_reports[7:t]);
  }else{
    target += neg_binomial_2_lpmf(cases | reports, rep_phi);
    target += neg_binomial_2_lpmf(weekly_cases[7:t] | weekly_reports[7:t], rep_phi);
  }

  // penalised priors for incubation period, and report delay
  target += normal_lpdf(inc_mean | inc_mean_mean, inc_mean_sd) * t;
  target += normal_lpdf(inc_sd | inc_sd_mean, inc_sd_sd) * t;
  target += normal_lpdf(rep_mean | rep_mean_mean, rep_mean_sd) * t;
  target += normal_lpdf(rep_sd | rep_sd_mean, rep_sd_sd) * t;
  
  ////////////////////////////////////////////////////// 
  // estimtate reproduction no from a branching process
  //////////////////////////////////////////////////////
  if (estimate_r) {
    
    // initial prior on R
    R[1] ~ gamma(r_alpha, r_beta);
    
   for (s in 2:t) {
     {
       // rescale previous R and overall sd prior for gamma
       real r_mean_ = R[s - 1];
       // assumption here is that todays R is like yesterdays
       real r_sd_ = 1;
       real r_alpha_ = (r_mean_ / r_sd_)^2;
       real r_beta_ = (r_sd_^2) / r_mean_;
       // reproduction number prior dependent on previous time step estimate
       R[s] ~ gamma(r_alpha_, r_beta_);
       }
     }
     
    // penalised_prior on generation interval
    target += normal_lpdf(gt_mean | gt_mean_mean, gt_mean_sd) * t;
    target += normal_lpdf(gt_sd | gt_sd_mean, gt_sd_sd) * t;
  
    // Likelihood of Rt given infections
    target += normal_lpdf(infections | branch_infections, 0.1);
  }
}
  
generated quantities {
  int imputed_infections[t];
  vector[t] prior_infections;
  
  // simulated infections - assume poisson (with negative binomial reporting)
  imputed_infections = poisson_rng(infections);
  
  // prior infections from backsampling
  prior_infections = shifted_cases;
}