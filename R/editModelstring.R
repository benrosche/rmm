# ================================================================================================ #
# Function editModelstring 
# ================================================================================================ #

editModelstring <- function(family, priors, mm, level1, level2, level3, DIR, modelfile) {
  
  # Unpack lists --------------------------------------------------------------------------------- #
  
  hm <- !is.null(level3$dat)
  
  l1vars <- level1$vars
  l2vars <- level2$vars
  l3vars <- level3$vars
  
  mmwfunction <- mm$mmwfunction
  mmwcoefstring <- mm$mmwcoefstring
  mmwconstraint <- mm$mmwconstraint
  mmwar <- mm$mmwar
  
  # Model family --------------------------------------------------------------------------------- #
  
  if(family=="Gaussian") {
    modelstring <- if(hm) readr::read_file(paste0(DIR, "/JAGS/Gaussian_l123.txt")) else readr::read_file(paste0(DIR, "/JAGS/Gaussian_l12.txt")) # levels
  } else if(family=="Weibull") {
    modelstring <- if(hm) readr::read_file(paste0(DIR, "/JAGS/Weibull_l123.txt")) else readr::read_file(paste0(DIR, "/JAGS/Weibull_l12.txt")) # levels
  } else if(family=="Cox") {
    modelstring <- if(hm) readr::read_file(paste0(DIR, "/JAGS/Cox_l123.txt")) else readr::read_file(paste0(DIR, "/JAGS/Cox_l12.txt")) # levels
  }
  
  win2unix <- function(str) {
    gsub("\r\n", "\n", str, fixed = TRUE)
  }
  
  modelstring <- win2unix(modelstring)
  
  # No covariates at level 1?
  if(length(l1vars)==0) { 
    modelstring <- stringr::str_remove(modelstring,  fixed("  for(x in 1:n.Xl1) {\n    b.l1[x] ~ dnorm(0,0.0001)\n    ppp.b.l1[x] <- step(b.l1[x])\n  }\n  \n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l1[i] <- inprod(X.l1[i,], b.l1) + re.l1[l1id[i]]"), "l1[i] <- re.l1[l1id[i]]")
  }
  
  # Auto-regressive l1 effect?
  if(mmwar==T) { 
    modelstring <- stringr::str_replace(modelstring, fixed("re.l1[l1id[i]]"), "re.l1[l1id[i], n.GPn[i]]")
    modelstring <- stringr::str_replace(modelstring, fixed("re.l1[i] ~ dnorm(0, tau.l1)\n  "), "re.l1[i,1] ~ dnorm(0, tau.l1)\n    for(j in 2:n.GPNi[ i ]) {\n      re.l1[i,j] ~ dnorm(re.l1[i,j-1], tau.l1)\n    }\n  ")
  }
  
  # No covariates at level 3?
  if(hm & length(l3vars)==0) { 
    modelstring <- stringr::str_remove(modelstring,  fixed("for(x in 1:n.Xl3) {\n    b.l3[x] ~ dnorm(0,0.0001)\n    ppp.b.l3[x] <- step(b.l3[x])\n  }\n  \n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l3[k] <- inprod(X.l3[k,], b.l3) + re.l3[k]"), "l3[k] <- re.l3[k]")
  }
  
  # Weight function
  if(mmwfunction != "uw[i] <- 1/X.w[i,1]") modelstring <- stringr::str_replace(modelstring, fixed("uw[i] <- 1/X.w[i,1]"), mmwfunction)
  if(mmwconstraint == 2) modelstring <- stringr::str_replace(modelstring, fixed("w[i] <- uw[i] / sum(uw[l1i1.l1[i]:l1i2.l1[i]])"), "w[i] <- uw[i] * n.l2/sum(uw[]) # rescale to sum up to n.l1 overall")
  if(stringr::str_detect(mmwfunction, "b.w\\[.\\]")) { 
    modelstring <- stringr::str_replace(modelstring, fixed("b.w\n"), paste0(mmwcoefstring, collapse = "")) 
  } else {
    modelstring <- stringr::str_replace(modelstring, fixed("b.w\n"), "") 
  }
  
  # Priors
  if(!is.null(priors)) {
    
    vars <- names(priors) # variables in list of priors
    prior <- priors       # prior for each variable
    
    # change priors in modelstring
    for(i in 1:length(vars)) {
      if(vars[i] %in% l1vars) {
        modelstring <- stringr::str_replace(modelstring, fixed("b.l1[b] ~ dnorm(0,0.0001)"), paste0("b.l1[b] ~ ", prior[i]))
      } else if(vars[i] %in% l2vars) {
        modelstring <- stringr::str_replace(modelstring, fixed("b.l2[b] ~ dnorm(0,0.0001)"), paste0("b.l2[b] ~ ", prior[i]))
      } else if(vars[i] %in% l3vars) {
        modelstring <- stringr::str_replace(modelstring, fixed("b.l3[b] ~ dnorm(0,0.0001)"), paste0("b.l3[b] ~ ", prior[i]))
      } else if(vars[i] %in% lwvars) {
        modelstring <- stringr::str_replace(modelstring, fixed("b.w[b] ~ dnorm(0,0.0001)"), paste0("b.w[b] ~ ", prior[i]))
      } else if(vars[i] == "tau.l1") {
        modelstring <- stringr::str_replace(modelstring, fixed("tau.l1 ~ dscaled.gamma(25, 1)"), paste0("tau.l1 ~ ", prior[i]))
      } else if(vars[i] == "tau.l2") { 
        modelstring <- stringr::str_replace(modelstring, fixed("tau.l2 ~ dscaled.gamma(25, 1)"), paste0("tau.l2 ~ ", prior[i]))
      } else if(vars[i] == "tau.l3") {
        modelstring <- stringr::str_replace(modelstring, fixed("tau.l3 ~ dscaled.gamma(25, 1)"), paste0("tau.l3 ~ ", prior[i]))
      } 
    }
  }
  
  if(modelfile==T) readr::write_file(modelstring, paste0(DIR, "/temp/modelstring.txt"))
  
  return(modelstring)
  
}

