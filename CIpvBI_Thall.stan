data {
  int<lower=0> x_se;
  int<lower=0> x_sp;
  int<lower=0> n_se;
  int<lower=0> n_sp;
  real<lower=0, upper=1> Prev;
  real mean_mu; 
  real<lower=0> sigma_mu; 
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
}
parameters{
  real<lower=0, upper=1> Se;
  real<lower=0, upper=1> Sp;
  real mu_se;
  real<lower=0> sigma_se;
  real mu_sp;
  real<lower=0> sigma_sp; 
}
transformed parameters{
  real rho_se;
  real rho_sp;
  real<lower=0, upper=1> PPV;
  real<lower=0, upper=1> NPV;
  
  rho_se = log(Se/(1-Se));
  rho_sp = log(Sp/(1-Sp));
  PPV = Se*Prev / ( Se*Prev + (1-Sp)*(1-Prev) );
  NPV = Sp*(1-Prev) / ( (1-Se)*Prev + Sp*(1-Prev) );
}
model{
  //layer1 prior
  x_se ~ binomial(n_se, Se);
  x_sp ~ binomial(n_sp, Sp);
  //layer2 prior
  rho_se ~ normal(mu_se,sigma_se);
  rho_sp ~ normal(mu_sp,sigma_sp);
  //layer3 prior
  mu_se ~ normal(mean_mu, sigma_mu);
  sigma_se ~ gamma(sigma_alpha, sigma_beta);
  mu_sp ~ normal(mean_mu, sigma_mu);
  sigma_sp ~ gamma(sigma_alpha, sigma_beta);
}