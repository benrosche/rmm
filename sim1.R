
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

# ================================================================================================ #
# Estimating the variance at each level
# ================================================================================================ #

rmm(sim_y ~ 1 + majority + investiture + mwc + mm(id(pid, gid), mmc(ipd+ fdep), mmw(w ~ 1/offset(n))) + hm(id=cid, type=RE), data=crDat(party=c(0,0,0), gov=c(1,3,3,1), country=c(3,0), weight=c(0,0,0), Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), level=1, seed=NULL))
