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
}


data {
  int t; // number of time steps
  int d; 
  int inc; 
  int samples;
  int wkd[t];
  int mon[t];
  int <lower = 0> cases[t];
  vector<lower = 0>[t] shifted_cases; 
  int model_type; //Type of model: 1 = Poisson otherwise negative binomial
  vector[d] delay[samples]; 
  vector[inc] incubation[samples];
}

transformed data{
  vector[d] rev_delay[samples];
  vector[inc] rev_incubation[samples];
  //Reverse the 
  for (h in 1:samples) {
    for (j in 1:d) {
      rev_delay[h][j] = delay[h][d - j + 1];
    }
   
    for (j in 1:inc) {
        rev_incubation[h][j] = incubation[h][inc - j + 1];
    }
  }
  
}

parameters{
  vector<lower = 0>[t] noise;
  real <lower = 0> phi; 
  real wkd_eff;
  real mon_eff;
}

transformed parameters {
  vector<lower = 0>[t] infections;
  vector<lower = 0>[t] onsets[samples];
  vector<lower = 0>[t] reports[samples];
  
  //Generation infections from median shifted cases and non-parameteric noise
  infections = shifted_cases .* noise;

  
  for(h in 1:samples) {
     // Onsets from infections
     onsets[h] = convolve(infections, rev_incubation[h]);
     
     // Reports from onsets
     reports[h] = convolve(onsets[h], rev_delay[h]);
     
     // Add reporting effects
     for (s in 1:t) {
      reports[h, s] = reports[h, s] + (1 + (wkd_eff * wkd[s]) + (mon_eff * mon[s]));
     }
  }
}

model {
  // Week effects
  wkd_eff ~ normal(0, 0.1);
  mon_eff ~ normal(0, 0.1);
  
  // Reporting overdispersion
  phi ~ exponential(1);

  // Noise on median shift
  for (i in 1:t) {
    noise[i] ~ normal(1, 0.2) T[0,];
  }
  
  for (h in 1:samples) {
    // Log likelihood of reports
     if (model_type == 1) {
       target +=  poisson_lpmf(cases | reports[h]);
     }else{
       target += neg_binomial_2_lpmf(cases | reports[h], phi);
     }
  }

}
  
generated quantities {
  int imputed_infections[t];
  imputed_infections = poisson_rng(infections);
}
