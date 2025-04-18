# ================================================================================================ #
# Function editModelstring 
# ================================================================================================ #

editModelstring <- function(family, priors, l1, l3, level1, level2, level3, weight, DIR, monitor, modelfile) {
  
  # Unpack variables ----------------------------------------------------------------------------- #
  
  lhs <- level2[["lhs"]]
  
  l1vars <- level1[["vars"]]
  l2vars <- level2[["vars"]]
  l3vars <- level3[["vars"]]
  wvars <- weight[["vars"]]
  wvars_p <- weight[["vars_p"]]
  wparams <- weight[["params"]]
  
  mm <- l1[["mm"]]
  mmwfunction <- l1[["mmwfunction"]]
  mmwconstraint <- l1[["mmwconstraint"]]
  mmwar <- l1[["mmwar"]]
  
  hm <- l3[["hm"]]
  l3type <- l3[["l3type"]]
  
  # Model family --------------------------------------------------------------------------------- #
  
  if(family=="Gaussian") {
    
    if(isTRUE(mm)&isTRUE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Gaussian_l123.txt"))
    } else if(isTRUE(mm)&isFALSE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Gaussian_l12.txt")) 
    } else if(isFALSE(mm)&isTRUE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Gaussian_l23.txt")) 
    } else if(isFALSE(mm)&isFALSE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Gaussian_l2.txt")) 
    }
    
  } else if(family=="Binomial") {
    
    if(isTRUE(mm)&isTRUE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Binomial_l123.txt"))
    }  else if(isTRUE(mm)&isFALSE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Binomial_l12.txt")) 
    } else if(isFALSE(mm)&isTRUE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Binomial_l23.txt")) 
    } else if(isFALSE(mm)&isFALSE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Binomial_l2.txt")) 
    }
    
  } else if(family=="Weibull") {
    
    if(isTRUE(mm)&isTRUE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Weibull_l123.txt")) 
    } else if(isTRUE(mm)&isFALSE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Weibull_l12.txt")) 
    } else if(isFALSE(mm)&isTRUE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Weibull_l23.txt")) 
    } else if(isFALSE(mm)&isFALSE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Weibull_l2.txt")) 
    }
    
  } else if(family=="Cox") {
    
    if(isTRUE(mm)&isTRUE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Cox_l123.txt")) 
    } else if(isTRUE(mm)&isFALSE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Cox_l12.txt")) 
    } else if(isFALSE(mm)&isTRUE(hm)) {
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Cox_l23.txt")) 
    } else if(isFALSE(mm)&isFALSE(hm)) { 
      modelstring <- readr::read_file(paste0(DIR, "/JAGS/Cox_l2.txt")) 
    }
    
  }
  
  win2unix <- function(str) {
    gsub("\r\n", "\n", str, fixed = TRUE)
  }
  
  modelstring <- win2unix(modelstring)
  
  # Add predicted values of the dependent variable?
  if(monitor) {
    if(family=="Gaussian") modelstring <- stringr::str_replace(modelstring, "(Y\\[j\\] \\~) (dnorm\\(mu\\[j\\], tau.l2\\))", "\\1 \\2\n    pred[j] ~ \\2")
    if(family=="Weibull" & length(lhs)==2) modelstring <- stringr::str_replace(modelstring, "(t\\[j\\] \\~) (dweib\\(shape, lambda\\[j\\]\\))", "\\1 \\2\n\n    pred[j] ~ \\2")
    if(family=="Weibull" & length(lhs)==3) modelstring <- stringr::str_replace(modelstring, "(t\\[j\\] \\~) (dweib\\(shape, lambda\\[j\\]\\))", "\\1 \\2\n\n    ones[j] ~ dinterval(pred[j], ct.lbub[j,])\n    pred[j] ~ \\2") # 2do: Remove?
    # 2do: predictions for Binomial and Cox model
  }
  
  # No covariates at level 1?
  if(mm & length(l1vars)==0) { 
    modelstring <- stringr::str_remove(modelstring,  fixed("  for(x in 1:n.Xl1) {\n    b.l1[x] ~ dnorm(0,0.0001)\n  }\n  \n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l1[i] <- inprod(X.l1[i,], b.l1) + re.l1[l1id[i]]"), "l1[i] <- re.l1[l1id[i]]")
  }
  
  # Auto-regressive l1 effect?
  if(mm & mmwar) { 
    modelstring <- stringr::str_replace(modelstring, fixed("re.l1[l1id[i]]"), "re.l1[l1id[i], n.GPn[i]]")
    modelstring <- stringr::str_replace(modelstring, fixed("re.l1[i] ~ dnorm(0, tau.l1)\n  "), "re.l1[i,1] ~ dnorm(0, tau.l1)\n    for(j in 2:n.GPNi[ i ]) {\n      re.l1[i,j] ~ dnorm(re.l1[i,j-1], tau.l1)\n    }\n  ")
  }
  
  # No covariates at level 3?
  if(hm & length(l3vars)==0) { 
    modelstring <- stringr::str_remove(modelstring,  fixed("for(x in 1:n.Xl3) {\n    b.l3[x] ~ dnorm(0,0.0001)\n  }\n  \n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l3[k] <- inprod(X.l3[k,], b.l3) + re.l3[k]"), "l3[k] <- re.l3[k]")
  }
  
  # Fixed effect at level 3?
  if(hm & l3type=="FE") {
    modelstring <- stringr::str_remove(modelstring,  fixed("for (k in 1:n.l3) {\n    re.l3[k] ~ dnorm(0, tau.l3)\n  }\n  \n"))
    modelstring <- stringr::str_remove(modelstring,  fixed("tau.l3 ~ dscaled.gamma(25, 1)\n  sigma.l3 <- 1/sqrt(tau.l3)\n  \n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l3[k] <- inprod(X.l3[k,], b.l3) + re.l3[k]"), "l3[k] <- inprod(X.l3[k,], b.l3)")
  }
  
  # Weight function
  if(mm) {
    
    # Translate weight function into JAGS format

    # Replace parameters in mmwfunction with JAGS code
    if(length(wparams)>0) {
      mmwpriors <- c()
      for(i in 1:length(wparams)) {
        mmwfunction <- str_replace_all(mmwfunction, fixed(wparams[i]), paste0("b.w[",i,"]")) 
        mmwpriors  <- append(mmwpriors, paste0("b.w[", i, "] ~ dnorm(0,0.0001)\n  ")) # collect priors
      }
    }

    # Replace variables
    for(i in 1:length(wvars)) {
      mmwfunction <- str_replace_all(mmwfunction, paste0("\\b", wvars[i], "\\b"), paste0("X.w[i,", i, "]")) # replace variables 
    }
    
    # Replace w
    mmwfunction <- stringr::str_replace(mmwfunction, "w ~", "uw[i] <- ")

    # Insert mmwfunction into modelstring
    modelstring <- stringr::str_replace(modelstring, fixed("uw[i] <- 1/X.w[i,1]"), mmwfunction)
    
    # Insert constraint into modelstring
    if(mmwconstraint == 2) {
      modelstring <- stringr::str_replace(modelstring, fixed("w[i] <- uw[i] / sum(uw[l1i1.l1[i]:l1i2.l1[i]])"), "w[i] <- uw[i] * n.l2/sum(uw[]) # rescale to sum up to n.l1")
    }
    
    # Insert priors into modelstring
    if(stringr::str_detect(mmwfunction, "b.w\\[.\\]")) { 
      modelstring <- stringr::str_replace(modelstring, fixed("b.w # placeholder\n"), paste0(mmwpriors, collapse = "")) 
    } else {
      modelstring <- stringr::str_replace(modelstring, fixed("b.w # placeholder\n"), "") 
    }
    
  }
  
  # Priors
  
  if(!is.null(priors)) {
    
    for(i in 1:length(priors)) {
    
      # Extract full param (e.g., b.w or b.w[1])
      full_param <- str_extract(priors[i], "^[^~<]+") %>% str_squish()
      
      # Check if full_param is indexed (e.g., b.w[1])
      if (str_detect(full_param, "\\[")) {
        
        # Replace prior for one parameter
        full_param_escaped <- str_replace_all(full_param, "(\\[|\\])", "\\\\\\1")
        pattern <- paste0("(?m)^([ \\t]*)", full_param_escaped, "\\s*(~|<-)\\s*[^\\n]*")
        modelstring <- str_replace(modelstring, pattern, paste0("\\1", priors[i]))
        
      } else {
        
        # Replace prior for all indices of a parameter
        base_param <- str_extract(priors[i], "^[^~<]+") %>% str_squish()
        operator <- str_extract(priors[i], "(~|<-)")
        rhs <- str_extract(priors[i], "(?<=~|<-).*") %>% str_squish()
        escaped <- str_replace_all(base_param, "(\\.|\\[|\\])", "\\\\\\1") # for regex
        # pattern <- paste0("(?m)^\\s*(", escaped, "\\[[^\\]]+\\])\\s*(~|<-)\\s*[^\\n]*")
        pattern <- paste0("(?m)^\\s*(", escaped, "(?:\\[[^\\]]+\\])?)\\s*(~|<-)\\s*[^\\n]*")
        modelstring <- str_replace_all(modelstring, pattern, paste0("\\1 ", operator, " ", rhs))

      }
      
    }
  
  }
  
  # Save or read modelstring
  if(isTRUE(modelfile)) {
    readr::write_file(modelstring, paste0(DIR, "/temp/modelstring.txt")) # save model to file
  } else if(!isFALSE(modelfile) & length(modelfile)>0) {
    modelstring <- readr::read_file(modelfile) # read model from file
  }
  
  return(modelstring)
  
}
