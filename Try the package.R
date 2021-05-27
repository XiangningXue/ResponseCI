install.packages("bdpv")
library(bdpv)


# Package example ---------------------------------------------------------
# 1) Example data: Mercaldo et al.(2007), Table VIII:
Tab8<-matrix(c(240, 178, 87, 288), ncol=2)
colnames(Tab8)<-c("Case","Control")
rownames(Tab8)<-c("ApoEe4plus","ApoEe4minus")
Tab8
# Assuming prevalence=0.03
BDtest(xmat=Tab8, pr=0.03, conf.level = 0.95)  
# Assuming prevalence=0.5
BDtest(xmat=Tab8, pr=0.5, conf.level = 0.95)

# 2) Experimental design acc. to Steinberg et al.(2009)
TEST<-nPV(se=c(0.76, 0.78, 0.80, 0.82, 0.84),
          sp=c(0.93, 0.94, 0.95, 0.96, 0.97),
          pr=0.0625, NPV0=0.98, PPV0=0.25, NPVpower = 0.8, PPVpower = 0.8,
          rangeP = c(0.10, 0.9), nsteps = 20, alpha = 0.05)
TEST
plotnPV(TEST, log="y", legpar=list(x=0.6))

# 3) Simulation of power and coverage probability
simPVmat(se=0.8, sp=0.95, pr=0.0625, n1=c(177, 181),
         n0=c(554, 87), NPV0=0.98, PPV0=c(0.4, 0.25))


# The Stamey paper example ------------------------------------------------
# example data: Stamey and Holt, Table 8 (page 108)
# Diseased
# Test D=1 D=0
# T=1 240 87
# T=0 178 288
#n1,n0: 418 375
# reproduce the results for the Bayes I method
# in Stamey and Holt (2010), Table 9, page 108
# assuming known prevalence 0.03
# ppv 0.0591, 0.0860
# npv 0.9810, 0.9850
CIpvBI( x1=c(240,178), x0=c(87,288), pr=0.03)
# assuming known prevalence 0.04
# ppv 0.0779, 0.1111
# npv 0.9745, 0.9800
CIpvBI( x1=c(240,178), x0=c(87,288), pr=0.04)
# compare with standard logit intervals
tab <- cbind( x1=c(240,178), x0=c(87,288))
tab
BDtest(tab, pr=0.03)
BDtest(tab, pr=0.04)
# reproduce the results for the Bayes II method
# in Stamey and Holt (2010), Table 9, page 108
CIpvBII( x1=c(240,178), x0=c(87,288),  shapespr=c(16,486))#16/(16+486) = 0.03187251 #sum 502
CIpvBII( x1=c(240,178), x0=c(87,288), shapespr=c(21,481))


#The simulation in the paper
prevalence = 0.05
Se = 0.55
n0 = n1 = 25
C = 0.95
set.seed(15213)
tab.2by2 = r.2by2(n0, n1, Se, C)
tab.2by2 = cbind(c(10, 15), c(0, 25))
knitr::kable(tab.2by2, caption = "A random table")
BDtest(tab.2by2, prevalence)

get.PPV.NPC(tab.2by2, prevalence)
