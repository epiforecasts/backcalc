data {
  int t; // number of time steps
  int d; 
  int inc; 
  int samples;
  int wkd[t];
  int mon[t];
  int <lower = 0> cases[t];
  real <lower = 0> shifted_cases[t]; 
  int model_type; //Type of model: 1 = Poisson otherwise negative binomial
  real delay[samples, d]; 
  real incubation[samples, inc];
}

parameters{
  real <lower = 0> noise[t];
  real wkd_eff;
  real mon_eff;
}

transformed parameters {
  real<lower = 0> infections[t];
  real<lower = 0> onsets[samples, t];
  real<lower = 0> reports[samples, t];
  
  for (s in 1:t) {
     infections[s] = shifted_cases[s] * noise[s];
  }
  
  for(h in 1:samples) {
     for (s in 1:t){
       onsets[h, s] = 0;
        for(i in 0:(min(s - 1, inc - 1))){
            onsets[h, s] += infections[s - i] * incubation[h, i + 1];
        }
     }
    
    for (s in 1:t){
      reports[h, s] = 0;
        for(i in 0:(min(s - 1, d - 1))){
             reports[h, s] += onsets[h, s - i] * delay[h, i + 1];
        }
     }
   }
}

model {
  real adj_reports[t];
  wkd_eff ~ normal(0, 0.1);
  mon_eff ~ normal(0, 0.1);
  
  for (i in 1:t) {
    noise[i] ~ normal(1, 0.1) T[0,];
  }
  
  for (h in 1:samples) {
    for (s in 1:t) {
      adj_reports[s] = reports[h, s] * (1 + (wkd_eff * wkd[s]) + (mon_eff * mon[s]));
    }
    
    target += poisson_lpmf(cases | adj_reports);
  }

}
  
generated quantities {
  int imputed_infections[t];
  imputed_infections = poisson_rng(infections);
}
