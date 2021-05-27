#CI test
#a lucky n 
z = qnorm(0.975)
n = 17
p = 0.5
set.seed(15213)
x.vec = rbinom(10000, n, p)
p.hat = x.vec /n
val1 = p.hat-z*sqrt(p.hat*(1-p.hat)/n)
val2 = p.hat+z*sqrt(p.hat*(1-p.hat)/n)
mean(p<val1)+mean(p>val2)
#an unlucky n
n = 18
p = 0.5
set.seed(15213)
x.vec = rbinom(10000, n, p)
p.hat = x.vec /n
val1 = p.hat-z*sqrt(p.hat*(1-p.hat)/n)
val2 = p.hat+z*sqrt(p.hat*(1-p.hat)/n)
mean(p<val1)+mean(p>val2)


n = 591
p = 0.005
set.seed(15213)
x.vec = rbinom(10000, n, p)
p.hat = x.vec /n
val1 = p.hat-z*sqrt(p.hat*(1-p.hat)/n)
val2 = p.hat+z*sqrt(p.hat*(1-p.hat)/n)
mean(p<val1)+mean(p>val2)

n = 592
p = 0.005
set.seed(15213)
x.vec = rbinom(10000, n, p)
p.hat = x.vec /n
val1 = p.hat-z*sqrt(p.hat*(1-p.hat)/n)
val2 = p.hat+z*sqrt(p.hat*(1-p.hat)/n)
mean(p<val1)+mean(p>val2)


z = qnorm(0.995)
n = 16
p = 0.5
set.seed(15213)
x.vec = rbinom(10000, n, p)
p.hat = x.vec /n
val1 = p.hat-z*sqrt(p.hat*(1-p.hat)/n)
val2 = p.hat+z*sqrt(p.hat*(1-p.hat)/n)
mean(p<val1)+mean(p>val2)

##exact CI
install.packages("LaplacesDemon")
library(LaplacesDemon)
