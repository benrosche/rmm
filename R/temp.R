# govsurv <-
#   data.frame(pid=c(10,1,2,3,4,5,6,7,8,9,10,1,9,4,6,8,3,9,9,2,4),
#              gid=rep(1:10, times=c(2, 3, 1, 1, 2, 2, 2, 3, 2, 3))) %>%
#   mutate(survival=rweibull(21, 1, 1),
#          event=rbinom(21, 1, 0.5),
#          Y=rnorm(21, 0, 1),
#          var1=rnorm(21, 0, 1),
#          var2=rnorm(21, 0, 1),
#          var3=rnorm(21, 0, 1),
#          var4=rnorm(21, 0, 1),
#          var5=rnorm(21, 0, 1)) %>%
#   group_by(gid) %>% add_count(name="N")
# save(govsurv, file = "data/govsurv.Rdata")
# 
# rmm(formula=Surv(survival, event) ~ 1 + var1 + var2 + var3 + mm(id(pid, gid), mmc(var4 + var5), mmw(w~1/var6, constraint=1)),
#    family="Weibull",
#    priors=list(var1="dnorm(0,0.0005)", var2="dunif(-100,100)", tau.l1="dscaled.gamma(25, 1)", sigma.l2="dexp(0.0001)"),
#    data=data)
