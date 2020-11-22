
rm(list=ls())

library(rmm)
library(survival)

setwd("rmm/")

source("create/createData.R")

# ================================================================================================ #
# Simulation function
# ================================================================================================ #

statSim <- function(iter, model, message=F) {
  
  # Iteration 1
  sink("null") # hides output
  fit <- eval(parse(text=model))
  sink()
  
  ret <- data.frame(cbind(i=1,rbind(fit$coefficients))) # first iteration
  
  if(message) message("1")
  
  if(any(names(fit) == "variable"))  colnames(ret) <- c("i", fit$variable)
    
  if(iter > 1) {
    # Iteration 2 to iter
    for(i in 2:iter) {
      sink("null") # hides output
      fit <- eval(parse(text=model))
      sink()
      ret[i,] <- cbind(i=i,rbind(fit$coefficients))
      if(message) message(i)
    }
    
  }

  ret <- ret %>% summarise_all(mean) %>% mutate(i=iter)
  ret<- list(ret, model)
  
  return(ret)
  
}

# HERE


# ================================================================================================ #
# Linear outcome (lm)
# ================================================================================================ #

# Estimation at level 1 -------------------------------------------------------------------------- #

# Level 2
statSim(10, "lm(y ~ 1 + majority + mwc, dat=crY(dat=crX(), party=c(0,0,0), gov=c(1,3,3,0), country=c(0,0), weight=c(0)))")

statSim(30, "lm(y ~ 1 + majority + mwc, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(0,0), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,10,0, 0,0,1),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(0,0), weight=c(0,0,0)))")

# Level 3
statSim(30, "lm(y ~ 1 + majority + mwc, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(0,1), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,10),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(0,1), weight=c(0,0,0)))")

statSim(30, "lm(y ~ 1 + majority + mwc + investiture, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(3,1), weight=c(0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,10),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(3,1), weight=c(0)))")

# Level 1
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3)), party=c(3,3,0), gov=c(1,0,0,0), country=c(0,0), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(10,0,0, 0,1,0, 0,0,1),3,3)), party=c(3,3,1), gov=c(1,0,0,0), country=c(0,0), weight=c(0,0,0)))")

statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(), party=c(3,3,1), gov=c(1,0,0,1), country=c(0,1), weight=c(0)))")

statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3)), party=c(0,0,0), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(10,0,0, 0,1,0, 0,0,1),3,3)), party=c(0,0,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")

# Weight function
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,2, 0,2,0, 2,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")

statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,2, 0,2,0, 2,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")

statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0))))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,2, 0,2,0, 2,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0))))")

# Comments:
# Level 2
# - party=c(0,0,0), gov=c(1,3,3,0), country=c(0,0), weightf=c(0): perfect estimation
# - party=c(0,0,0), gov=c(1,3,3,1), country=c(0,0), weightf=c(0): no attenuation bias going from re.gov=1 to 3
# Level 3
# - party=c(0,0,0), gov=c(1,3,3,1), country=c(0,1), weight=c(0))): small overestimation of l2 predictors going from re.country=1 to 5
# - party=c(0,0,0), gov=c(1,3,3,1), country=c(1,1), weight=c(0))): sizable underestimation of l3 predictors going from re.country=1 to 5 (perfect estimation when re.country=0)
# Level 1
# - party=c(3,3,0), gov=c(1,0,0,0), country=c(0,0), weight=c(0))): sizable underestimation of l1 predictors that does not depend on re.party, re.gov, or re.country
# - re.party does change the effects of l2 and l3 predictors in ambiguous ways
# - weight>0: without correlation between party predictors and higher-level predictors, we have ambiguous effects on outcome
#   -> positive effects seem to decrease higher-level effects and increase party-level effects
#   -> negative effects seem to increase higher-level effects and decrease party-level effects

# 2do: Make graphs out of this ...

# Estimation at level 2 -------------------------------------------------------------------------- #

statSim(10, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(), party=c(3,3,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0)), level=2))")
statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)), level=2))")

statSim(30, "lm(y ~ 1 + majority + mwc + investiture + ipd + fdep, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)), level=2))")

# ================================================================================================ #
# Survival outcome (coxph)
# ================================================================================================ #

# Estimation at level 1 -------------------------------------------------------------------------- #

# Level 2+3
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(0,0,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0))))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(0,0,0), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(0,0,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")

# Level 1
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(3,3,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")

statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")

# Weights
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")

statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,2, 0,2,0, 2,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,2, 0,2,0, 2,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,2, 0,2,0, 2,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")

# Comments:
# - only l2 + l3 predictors: perfect prediction
# - adding l2+l3 RE seriously attenuates l2+l3 predictor effects
# - adding l1 RE also attenuatesl2+l3 predictor effects
# - l1 predictors seriously underestimated, especially when including RE
# - adding weight regressors has most significant attenuating effect when negative 

# Estimation at level 2 -------------------------------------------------------------------------- #

statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0), level=2))")

statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0), level=2))")

statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0), level=2))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0), level=2))")

# Comments:
# - Party effects close to 0 when w=-10

# ================================================================================================ #
# Linear outcome (RMM vs lm)
# ================================================================================================ #

statSim(1, "rmm(y ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), data=crY(dat=crX(), party=c(0,0,0), gov=c(1,3,3,1), country=c(3,0), weight=c(0,0,0)))")

statSim(1, "rmm(y ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0), level=1))")
statSim(30, "lm(y ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0), level=2))")

statSim(1, "rmm(y ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0), level=2))")

statSim(10, "rmm(y ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")
statSim(30, "lm(y ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,2,0, 2,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0), level=2))")

# Comments:

# ================================================================================================ #
# Survival outcome (RMM vs coxph)
# ================================================================================================ #

# No party level, no random effects
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(), party=c(0,0,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3)), party=c(0,0,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0), level=2))")

# Standard model
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0)))")
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0), level=2))")

# Weights on ------------------------------------------------------------------------------------- #

# Level 1 predictor (positive)
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(10,0,0), level=2))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(10,0,0)))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n)^exp(-(pseatrel)),c=2)) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(10,0,0)))")

# Level 2 predictor (positive)
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0), level=2))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n)^exp(-hetero))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,10,0)))")

# Level 3 predictor (positive)
statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,10), level=2))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,10)))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n)^exp(-pmpower))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3), seed=1), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,10)))")

statSim(30, "coxph(Surv(survtime, event) ~ 1 + majority + mwc + ipd + fdep + investiture, dat=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0), level=2))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")
statSim(1, "rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n)^exp(-hetero))) + (1|cid), family=\"Weibull\", iter=5000, burnin=500, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))")


rmm(Surv(survtime, event) ~ 1 + majority + mwc + investiture + mm(id(pid, gid), mmc(ipd + fdep), mmw(w ~ 1/offset(n)^exp(-hetero))) + (1|cid), family="Weibull", iter=1000, burnin=100, monitor=F, data=crY(dat=crX(Sigma=matrix(c(2,0,0, 0,2,0, 0,0,2),3,3)), party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,-10,0)))

# Check level 3 again. It seems like the weight coefficient is not correctly estimated. But l2 and l1 are?

# Comments:
# - If there is no variation at l2 (ie perfect prediction), I get "Error in update.jags(model, n.iter, ...) : LOGIC ERROR: Non-finite boundary in truncated normal"
# - DIC can be negative: https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/DIC-slides.pdf
