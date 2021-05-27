rm(list=ls())
library(runjags)
library(coda)
library(bdpv)

## Test outcome
# cnt <- as.data.frame(matrix (c(240,87,178,288),2,2,byrow = T))
cnt <- as.data.frame(matrix (c(200,20,10,300),2,2,byrow = F))

rownames(cnt) <- c("Test positive","Test negative")
colnames(cnt) <- c("True positive","True negative")

cnt

## Prevalence
# prevalence <- 0.03
prevalence <- 0.003


## point estimates
Se=cnt[1,1]/(cnt[1,1]+cnt[2,1])
Sp=cnt[2,2]/(cnt[2,2]+cnt[1,2])
Prev=prevalence
PPV <- Se*Prev / ( Se*Prev + (1-Sp)*(1-Prev) )
NPV <- Sp*(1-Prev) / ( (1-Se)*Prev + Sp*(1-Prev) )

Se
Sp
PPV
NPV

binom.test(cnt[1,1], cnt[1,1]+cnt[2,1])
binom.test(cnt[2,2], cnt[2,2]+cnt[1,2])

## CI for PPV/NPV by Stamey 2010 beta-binomial model 1
CI_BB1 <- CIpvBI( x0=cnt$`True negative`, x1=cnt$`True positive`, pr=prevalence)

## CI for PPV/NPV by Mercaldo 2010 Wald method 
CI_NPV <- CIlnpv(x0=cnt$`True negative`, x1=cnt$`True positive`, p=prevalence,
       conf.level = 0.95, alternative = "two.sided")

CI_PPV <- CIlppv(x0=cnt$`True negative`, x1=cnt$`True positive`, p=prevalence,
       conf.level = 0.95, alternative = "two.sided")



#######################################
## Revised model based on Thall 2003 ##
######## Kui Shen, 5/19/2021 ##########
#######################################


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
data <- list(x_se=cnt[1,1], n_se=cnt[1,1]+cnt[2,1], 
             x_sp=cnt[2,2], n_sp=cnt[2,2]+cnt[1,2],
             mean.Mu=0, perc.Mu=0.0001, 
             tau.alpha=0.001, tau.beta=0.001,
             Prev=prevalence)


inits1 <- dump.format((list(mu_se=0.01, tau_se=.01,mu_sp=0.01, tau_sp=.01,.RNG.name="base::Wichmann-Hill",.RNG.seed=1)))
inits2 <- dump.format((list(mu_se=0.02, tau_se=.02,mu_sp=0.01, tau_sp=.01,.RNG.name="base::Wichmann-Hill",.RNG.seed=2)))


set.seed(42)
params <- c('Se','rho_se','mu_se','tau_se','Sp','rho_sp','mu_sp','tau_sp',"PPV","NPV")

posterior <- run.jags(modelstring,
                      data=data,
                      monitor=params,
                      n.chains = 2,
                      inits=c(inits1,inits2),
                      adapt=10000,
                      burnin=10000,
                      sample=50000,
                      thin=5)

res <- round(summary(posterior),4)
PPV_est <- res["PPV",1:4]
NPV_est <- res["NPV",1:4]
res

binom.test(cnt[1,1], cnt[1,1]+cnt[2,1])
binom.test(cnt[2,2], cnt[2,2]+cnt[1,2])



##################
## Comparisions ##

PPV_est
NPV_est

CI_PPV
CI_NPV

CI_BB1
###########################


## diagnostics

plot(posterior,vars="Se")
plot(posterior,vars="Sp")

traceplot(posterior$mcmc)
autocorr(posterior$mcmc)
effectiveSize(posterior$mcmc)
gelman.diag(posterior$mcmc)



