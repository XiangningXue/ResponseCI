library(rstan)
# library(ggplot2)
# mean.paras = c(0, 0.0001)
# tau.paras = c(0.001, 0.001)
# mean.paras2 = c(0, 10^4)
# sigma.paras = c(10^3, 10^3)

CIpvBI_Thall_stan = function(x1,x0, pr, conf.level = 0.95, B = 5000, 
                        mean.paras = c(0, 10^4),  sigma.paras = c(10^3, 10^3)){
  
  #rstan code
  model.stan = "
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
"
observed.data = list("x_se" = x1[1], 
                     "x_sp" = x0[2], 
                     "n_se" = sum(x1), 
                     "n_sp" = sum(x0), 
                     "Prev" = pr,
                     "mean_mu" = mean.paras[1], 
                     "sigma_mu" = mean.paras[2], 
                     "sigma_alpha" = sigma.paras[1], 
                     "sigma_beta" = sigma.paras[2])

stanmodel = stan_model(model_code = model.stan, model_name = "CIpvBI_Thall_stan")
stanfit <- sampling(stanmodel, data = observed.data,
                    chains = 4, cores = 4, iter = B, warmup = B/5)
# stanfit = stan(file = "CIpvBI_Thall.stan",model_name = "CIpvBI_Thall_stan",data = observed.data,
#                chains = 4, cores = 4, iter = B, warmup = B/5) 
# print(stanfit)
# results<- extract(stanfit)
# ggplot(data.frame(results$PPV), aes(x = results$PPV)) + 
#   xlab("PPV") + 
#   geom_density(alpha = 0.5) +
#   theme_bw()
# ggplot(data.frame(posterior$mcmc[[1]][, "PPV"]), aes(x = posterior$mcmc[[1]][, "PPV"])) + 
#   xlab("PPV") + 
#   geom_density(alpha = 0.5) +
#   theme_bw()
CI.alpha = 1-conf.level
probs.CI = c(CI.alpha/2, 1-CI.alpha/2)
stanfit.summary = summary(stanfit, probs = probs.CI)$summary
PPV.est = stanfit.summary["PPV", "mean"]
NPV.est = stanfit.summary["NPV", "mean"]

PPV.CI = stanfit.summary["PPV", c(4, 5)]
NPV.CI = stanfit.summary["NPV", c(4, 5)]

CIs = rbind(NPV.CI, PPV.CI); rownames(CIs) = c("NPV", "PPV")
ests = c(NPV.est, PPV.est); names(ests) = c("NPV", "PPV")
return(list(conf.int = CIs, 
            estimate = ests))

  
}

#debug 
stan.model = "
transformed data {
  int N = 5;
}
parameters {
  vector[N] theta;
}
model {
  1 ~ categorical(softmax(theta));
}
"
stanmodel = stan_model(model_code = stan.model, model_name = "CIpvBI_Thall_stan")