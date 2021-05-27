#' Title
#'
#' @param x1 A vector of two (integer) values, specifying the observed number of positive (x1[1]) and negative (x1[2]) test results in the group of true positives.
#' @param x0 A vector of two (integer) values, specifying the observed number of positive (x0[1]) and negative (x0[2]) test results in the group of true negatives.
#' @param pr A single numeric value between 0 and 1, defining an assumed fixed (known) prevalence (for CIpvBI), where prevalence is the proportion of positives in the population.
#' @param conf.level The confidence level, a single numeric value between 0 and 1, defaults to 0.95
#' @param B A single integer, the number of samples from the posterior to be drawn.
#' @param mean.paras A vector of two values for mean and precision of the prior of mean of the logit of sensitivity and specificity
#' @param tau.paras A vector of two values for two parameters pass to dgamma of the prior of variance of the logit of sensitivity and specificity
#'
#' @return
#' @export
#'
#' @examples
CIpvBI_Thall = function(x1,x0, pr, conf.level = 0.95, B = 5000, 
                        mean.paras = c(0, 0.0001),  tau.paras = c(0.001, 0.001)){
  
  modelstring="
model{

  ## Sensitivity
  x_se ~ dbin(Se,n_se) 
  logit(Se) <- rho_se
  rho_se ~ dnorm(mu_se,tau_se) 

  # Priors
  mu_se ~ dnorm(mean.Mu, perc.Mu) 
  tau_se ~ dgamma(tau.alpha, tau.beta) 
  
  
  ## Specificity
  x_sp ~ dbin(Sp,n_sp) 
  logit(Sp) <- rho_sp
  rho_sp ~ dnorm(mu_sp,tau_sp) 

  # Priors
  mu_sp ~ dnorm(mean.Mu, perc.Mu) 
  tau_sp ~ dgamma(tau.alpha, tau.beta) 
  
  ## PPV and NPV
  PPV <- Se*Prev / ( Se*Prev + (1-Sp)*(1-Prev) )
  NPV <- Sp*(1-Prev) / ( (1-Se)*Prev + Sp*(1-Prev) )

}
"


# data set
data <- list(x_se=x1[1], n_se=sum(x1), 
             x_sp=x0[2], n_sp=sum(x0),
             mean.Mu=mean.paras[1], perc.Mu=mean.paras[2], 
             tau.alpha=tau.paras[1], tau.beta=tau.paras[2],
             Prev=pr)

#initial values for Bayesian fitting
inits1 <- dump.format((list(mu_se=0.01, tau_se=.01,mu_sp=0.01, tau_sp=.01,.RNG.name="base::Wichmann-Hill",.RNG.seed=1)))
inits2 <- dump.format((list(mu_se=0.02, tau_se=.02,mu_sp=0.01, tau_sp=.01,.RNG.name="base::Wichmann-Hill",.RNG.seed=2)))


set.seed(42)
params <- c('Se','rho_se','mu_se','tau_se','Sp','rho_sp','mu_sp','tau_sp',"PPV","NPV")

posterior <- run.jags(modelstring,
                      data=data,
                      monitor=params,
                      n.chains = 2,
                      inits=c(inits1,inits2),
                      adapt=1000,
                      burnin=as.integer(B/5),
                      sample=B,
                      thin=5)
posterior.summary = summary(posterior, confidence = conf.level)
PPV.est = posterior.summary["PPV", "Mean"]
NPV.est = posterior.summary["NPV", "Mean"]

PPV.CI = posterior.summary["PPV", c(1, 3)]
NPV.CI = posterior.summary["NPV", c(1, 3)]

CIs = rbind(NPV.CI, PPV.CI); rownames(CIs) = c("NPV", "PPV")
ests = c(NPV.est, PPV.est); names(ests) = c("NPV", "PPV")
return(list(conf.int = CIs, 
            estimate = ests))

}