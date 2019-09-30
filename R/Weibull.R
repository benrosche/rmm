Model_Weibull <- function() {
  
  # Likelihood
  for(j in 1:n.l2) {
    
    # Weibull regression
    censored[j] ~ dinterval(t[j], t.cen[j])  
    t[j] ~ dweib(sigma.l2, mu[j]) # shape, scale
    
    # log(mu) = aggregated(l1) + l2
    log(mu[j]) <- inprod(l1[l1i1[j]:l1i2[j]], w[l1i1[j]:l1i2[j]]) + l2[j] 
    
    # Level 2: design matrix and weights
    l2[j] <- inprod(X.l2[j,], b.l2)
  }
  
  # Level 1: design matrix and random effect
  for(i in 1:n.l1) {
    l1[i] <- inprod(X.l1[i,], b.l1) + e.l1[l1id[i]] 
    
    uw[i]
    
    w[i] 
  }
  
  for (i in 1:n.ul1) {
    e.l1[i] ~ dnorm(0, tau.l1)
  }
  
  # Priors design matrix
  for(b in 1:n.Xl1) {
    b.l1[b] ~ dnorm(0,0.0001)
  }
  
  for(b in 1:n.Xl2) {
    b.l2[b] ~ dnorm(0,0.0001)
  }
  
  for(b in 1:n.Xlw) {
    b.w[b] ~ dnorm(0,0.0001)
  }
  
  # Priors variance terms
  tau.l1 ~ dscaled.gamma(25, 1)
  sigma.l1 <- 1/sqrt(tau.l1)
  
  sigma.l2 ~ dexp(0.0001) 
  
  # Posterior predictive p-values
  for (x in 1:n.Xl1) {             
    ppp.bl1[x] <- step(b.l1[x])
  }
  for (x in 1:n.Xl2) {             
    ppp.bl2[x] <- step(b.l2[x])
  }
}

