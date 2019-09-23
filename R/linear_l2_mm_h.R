modelstring <- function() {
  
  # Likelihood
  for(j in 1:n.l2) {
    
    # Outcome: level 2, linear function, level12=MM, l23=H
    Y[j] ~ dnorm(mu[j], tau.l2)
    
    mu[j] <- inprod(l1[i1[j]:i2[j]], w[i1[j]:i2[j]]) + l2 + l3
    
    # Level 1: design matrix and random effect
    for(i in 1:n.l1) {
      l1[i] <- inprod(X.l1[i,], b.l1) + e.l1[i] 
    }
    
    for (i in 1:n.ul1) {
      e.l1[i] ~ dnorm(0, tau.l1)
    }
    
    # Level 2: design matrix and weights
    for(j in 1:n.l2) {
      
      l2[j] <- inprod(X.l2[j,], b.l2)
      
      # weights
      for (h in (i1[j]):(i2[j])) {
        w[h] <- 1 / nl1[j]
      }
    }
    
    # Level 3: design matrix and random effect
    for(k in 1:n.l3) {
      l3[k] <- inprod(X.l3[k,], b.l3) + e.l3[k] 
    }
    
    for (k in 1:n.ul3) {
      e.l3[k] ~ dnorm(0, tau.l3)
    }
  
    # Priors design matrix
    for(b in 1:n.Xl1) {
      b.l1[b] ~ dnorm(0,0.0001)
    }
      
    for(b in 1:n.Xl2) {
      b.l2[b] ~ dnorm(0,0.0001)
    }
      
    for(b in 1:n.Xl3) {
      b.l3[b] ~ dnorm(0,0.0001)
    }
    
    # Priors variance terms
    tau.l1 ~ dscaled.gamma(25, 1)
    sigma.l1 <- 1/sqrt(tau.l1)
    
    tau.l2 ~ dscaled.gamma(25, 1)
    sigma.l2 <- 1/sqrt(tau.l2)
    
    tau.l3 ~ dscaled.gamma(25, 1)
    sigma.l3 <- 1/sqrt(tau.l3)
    
    # Posterior predictive p-values
    for (x in 1:n.Xl1) {             
      ppp.bl1[x] <- step(b.l1[x])
    }
    for (x in 1:n.Xl2) {             
      ppp.bl2[x] <- step(b.l2[x])
    }
    for (x in 1:n.Xl3) {             
      ppp.bl3[x] <- step(b.l3[x])
    }
}