functions {
  // Convolve a pdf and case vector using matrix multiplication
  vector convolve(vector cases, vector pdf) {
    int t = num_elements(cases);
    matrix[t, t] delay_mat = rep_matrix(0, t, t);
    int max_pdf = num_elements(pdf);
    row_vector[max_pdf] row_pdf = to_row_vector(pdf);
    vector[t] convolved_cases;
    
    for (s in 1:t) {
      int max_length = min(s, max_pdf);
      delay_mat[s, (s - max_length + 1):s] = row_pdf[(max_pdf - max_length + 1):max_pdf];
    }
  
   convolved_cases = delay_mat * to_vector(cases);

   return convolved_cases;
  }

  real discretised_lognormal_pmf(int y, real mu, real sigma) {

    return((lognormal_cdf(y, mu, sigma) - lognormal_cdf(y - 1, mu, sigma)));
  }
}


data {
  int t;                             // number of time steps
  int day_of_week[t];                // day of the week indicator (1 - 7)
  int <lower = 0> cases[t];          // observed cases
  vector<lower = 0>[t] shifted_cases;// mean delay + incubation shifted cases (smoothed)
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
  int model_type;                    //type of model: 1 = poisson otherwise negative binomial
}

transformed data{
  int<lower = 0> weekly_cases[t];    // weekly observed cases
  int<lower = 0> cum_cases[t];       // cumulative cases
  real r_alpha;                      //alpha parameter of the R gamma prior
  real r_beta;                       //beta parameter of the R gamma prior

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
  
}
parameters{
  vector<lower = 0>[t] noise;        // noise on the mean shifted observed cases
  real <lower = 0> inc_mean;         // mean of incubation period
  real <lower = 0> inc_sd;           // sd of incubation period
  real <lower = 0> rep_mean;         // mean of reporting delay
  real <lower = 0> rep_sd;           // sd of incubation period
  real<lower = 0> rep_phi;               // overdispersion of the reporting process
  vector[7] day_of_week_eff;         // day of week reporting effect
}

transformed parameters {
  vector[max_rep] rev_delay;         // reversed report delay pdf
  vector[max_inc] rev_incubation;    // reversed incubation period pdf
  vector<lower = 0>[t] infections;   // infections over time
  vector<lower = 0>[t] onsets;       // onsets over time
  vector<lower = 0>[t] reports;      // reports over time
  vector<lower = 0>[t] cum_reports;  // cumulative reported cases (unadjusted for reporting)
  vector<lower = 0>[t] weekly_reports; //weekly reported cases

   
  // reverse the distributions to allow vectorised access
    for (j in 1:max_rep) {
      rev_delay[j] =
        discretised_lognormal_pmf(max_rep - j + 1, inc_mean, inc_sd);
        }
   
    for (j in 1:max_inc) {
      rev_incubation[j] =
        discretised_lognormal_pmf(max_inc - j + 1, rep_mean, rep_sd);
    }

  // generate infections from median shifted cases and non-parameteric noise
  infections = shifted_cases .* noise;

  // onsets from infections
  onsets = convolve(infections, rev_incubation);
     
  // reports from onsets
  reports = convolve(onsets, rev_delay);
     
 // calculate Cumulative reports
  cum_reports = cumulative_sum(reports);
  
  for (s in 1:t) {
    // calculate weekly reports
    weekly_reports[s] = cum_reports[s] - cum_reports[max(1, s - 7)];
    // add reporting effects
    reports[s] *= day_of_week_eff[day_of_week[s]];
    }
}

model {
  // day of week reporting effect
  for (j in 1:7) {
    day_of_week_eff[j] ~ normal(1, 0.2) T[0,];
  }
  
  // reporting overdispersion
  rep_phi ~ exponential(1);

  // noise on median shift
  for (i in 1:t) {
    noise[i] ~ normal(1, 0.5) T[0,];
  }
  
  // daily cases given reports
  if (model_type == 1) {
    target += poisson_lpmf(cases | reports);
  }else{
    target += neg_binomial_2_lpmf(cases | reports, rep_phi);
  }
  
  // weekly cases given weekly reports
  target += poisson_lpmf(weekly_cases[7:t] | weekly_reports[7:t]);


  // penalised priors
  target += normal_lpdf(inc_mean | inc_mean_mean, inc_mean_sd) * t;
  target += normal_lpdf(inc_sd | inc_sd_mean, inc_sd_sd) * t;
  target += normal_lpdf(rep_mean | rep_mean_mean, rep_mean_sd) * t;
  target += normal_lpdf(rep_sd | rep_sd_mean, rep_sd_sd) * t;
}
  
generated quantities {
  int imputed_infections[t];
  
    // Simulated infections
  if (model_type == 1) {
    imputed_infections = poisson_rng(infections);
  }else{
    imputed_infections = neg_binomial_2_rng(infections, rep_phi);
  }

}
