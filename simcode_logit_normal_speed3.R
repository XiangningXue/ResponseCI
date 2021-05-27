setwd("/home/yuf31/xxn/ResponseCI")
source("CIpvBI_Thall_stan.R")
r.2by2 = function(n1, n0, S, C){
  x11 = rbinom(1, n1, S)
  x00 = rbinom(1, n0, C)
  tab.2by2 = cbind.data.frame(c(x11, n1-x11), c(n0-x00, x00))
  colnames(tab.2by2) = c("Disease1", "Disease0")
  rownames(tab.2by2) = c("Test1", "Test0")
  return(tab.2by2)
}
library(bdpv)
library(parallel)
library(ggplot2)

#simulation
n1 = seq(20, 100, by = 20); #n0 = seq(20, 100, by = 20)
prev = c(0.005, 0.05, seq(0.1, 0.9, by = 0.2)) 
Se = seq(0.1, 0.9, by = 0.2)
Sp = seq(0.1, 0.9, by = 0.2)
params = expand.grid(n1 = n1, prev = prev, Se = Se, Sp = Sp, rep = seq_len(1000))

time.start = Sys.time()
res = mclapply(1:nrow(params), function(a){
  n1 = params$n1[a]
  n0 = params$n1[a]
  prev = params$prev[a]
  Se = params$Se[a]
  Sp = params$Sp[a]
  
  set.seed(params$rep)
  tab.2by2 = r.2by2(n1, n0, Se, Sp)
  
  logit.ppv = CombCIppv(x0 = tab.2by2$Disease0, x1 = tab.2by2$Disease1, prev, conf.level = 0.95, alternative = "two.sided")
  logit.npv = CombCInpv(x0 = tab.2by2$Disease0, x1 = tab.2by2$Disease1, prev, conf.level = 0.95, alternative = "two.sided")
  
  Bayes1 = CIpvBI(x0 = tab.2by2$Disease0, x1 = tab.2by2$Disease1, prev)
  Bayes2 = CIpvBII(x0 = tab.2by2$Disease0, x1 = tab.2by2$Disease1, prev, xpr = c(ceiling(500*prev), ceiling(500*(1-prev))))
  
  Bayes3 = CIpvBI_Thall_stan(x0 = tab.2by2$Disease0, x1 = tab.2by2$Disease1, prev)
  
  one.row = params[a, ]
  one.row$PPV.truth = Se*prev/(Se*prev+(1-Sp)*(1-prev))
  one.row$NPV.truth = 1-Sp*(1-prev)/((1-Se)*prev+1-Sp*(1-prev))
  res.out = rbind.data.frame(cbind.data.frame(one.row, method = "logit", 
                                              PPV.est = logit.ppv$estimate, PPV.lower = logit.ppv$conf.int[1], PPV.upper = logit.ppv$conf.int[2], 
                                              NPV.est = logit.npv$estimate, NPV.lower = logit.npv$conf.int[1], NPV.upper = logit.npv$conf.int[2]), 
                             cbind.data.frame(one.row, method = "Bayes_beta_known_prev", 
                                              PPV.est = Bayes1$estimate["PPV"], PPV.lower = Bayes1$conf.int["PPV", 1], PPV.upper = Bayes1$conf.int["PPV", 2],
                                              NPV.est = Bayes1$estimate["NPV"], NPV.lower = Bayes1$conf.int["NPV", 1], NPV.upper = Bayes1$conf.int["NPV", 2]), 
                             cbind.data.frame(one.row, method = "Bayes_beta_unknown_prev", 
                                              PPV.est = Bayes2$estimate["PPV"], PPV.lower = Bayes2$conf.int["PPV", 1], PPV.upper = Bayes2$conf.int["PPV", 2],
                                              NPV.est = Bayes2$estimate["NPV"], NPV.lower = Bayes2$conf.int["NPV", 1], NPV.upper = Bayes2$conf.int["NPV", 2]), 
                             cbind.data.frame(one.row, method = "Bayes_logit_norm", 
                                              PPV.est = Bayes3$estimate["PPV"], PPV.lower = Bayes3$conf.int["PPV", 1], PPV.upper = Bayes3$conf.int["PPV", 2],
                                              NPV.est = Bayes3$estimate["NPV"], NPV.lower = Bayes3$conf.int["NPV", 1], NPV.upper = Bayes3$conf.int["NPV", 2]))
  res.out$PPV.cover = (res.out$PPV.truth>res.out$PPV.lower)&(res.out$PPV.truth<res.out$PPV.upper)
  res.out$PPV.CI.width = res.out$PPV.upper-res.out$PPV.lower
  res.out$NPV.cover = (res.out$NPV.truth>res.out$NPV.lower)&(res.out$NPV.truth<res.out$NPV.upper)
  res.out$NPV.CI.width = res.out$NPV.upper-res.out$NPV.lower
  print(paste0(a/nrow(params)))
  tdiff = Sys.time()-time.start
  units(tdiff) = "mins"
  print(paste0("Time spent = ", tdiff, " mins."))
  print(paste0("Time remain = ", (tdiff)/a*(nrow(params)-a), " mins."))
  
  return(res.out)
}, mc.cores = 10)
saveRDS(res, "simres_logit_normal.rds")