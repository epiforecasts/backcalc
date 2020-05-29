data {
  int t; // number of time steps
  int d; 
  int inc; 
  int samples;
  int <lower = 0> cases[t];
  real <lower = 0> shifted_cases[t]; 
  int model_type; //Type of model: 1 = Poisson otherwise negative binomial
  real delay[samples, d]; 
  real incubation[samples, inc];
}

parameters{
  real<lower = 0> infections[t];
  real <lower = 0> noise[t];
}

transformed parameters {
  real<lower = 0> onsets[samples, t];
  real<lower = 0> reports[samples, t];
  
  for(h in 1:samples) {
     for (s in 1:t){
       onsets[h, s] = 0;
        for(i in 0:(min((s - 1), inc - 1))){
            onsets[h, s] += infections[s - i] * incubation[h, i + 1];
        }
     }
    
    for (s in 1:t){
      reports[h, s] = 0;
        for(i in 0:(min((s - 1), d - 1))){
             reports[h, s] += onsets[h, s - i] * delay[h, i + 1];
        }
     }
   }
}

model {

  for (i in 1:t) {
    noise[t] ~ normal(0, 0.1) T[0,];
    infections[t] ~ normal(shifted_cases[t], shifted_cases[t] * noise[t]) T[0,];
  }
  
  for (h in 1:samples) {
    target += poisson_lpmf(cases | reports[h, 1:t]);
  }

}
  
generated quantities {
  int imputed_infections[t];
  imputed_infections = poisson_rng(infections);
}
