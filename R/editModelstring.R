# ================================================================================================ #
# Function editModelstring 
# ================================================================================================ #

editModelstring <- function(family, priors, l1, l3, level1, level2, level3, DIR, monitor, modelfile) {
  
  # Unpack lists --------------------------------------------------------------------------------- #
  
  lhs <- level2$lhs
  
  l1vars <- level1$vars
  l2vars <- level2$vars
  l3vars <- level3$vars
  
  mm <- l1$mm
  mmwfunction <- l1$mmwfunction
  mmwcoefstring <- l1$mmwcoefstring
  mmwconstraint <- l1$mmwconstraint
  mmwar <- l1$mmwar
  
  hm <- l3$hm
  l3type <- l3$l3type
  
  # Model family --------------------------------------------------------------------------------- #
  
  if(family=="Gaussian") {
    modelstring <- if(isTRUE(mm)&isTRUE(hm)) readr::read_file(paste0(DIR, "/JAGS/Gaussian_l123.txt")) else if(isTRUE(mm)&isFALSE(hm)) readr::read_file(paste0(DIR, "/JAGS/Gaussian_l12.txt")) else if(isFALSE(mm)&isTRUE(hm)) readr::read_file(paste0(DIR, "/JAGS/Gaussian_l23.txt")) else if(isFALSE(mm)&isFALSE(hm)) readr::read_file(paste0(DIR, "/JAGS/Gaussian_l2.txt"))
  } else if(family=="Weibull") {
    modelstring <- if(isTRUE(mm)&isTRUE(hm)) readr::read_file(paste0(DIR, "/JAGS/Weibull_l123.txt")) else if(isTRUE(mm)&isFALSE(hm)) readr::read_file(paste0(DIR, "/JAGS/Weibull_l12.txt")) else if(isFALSE(mm)&isTRUE(hm)) readr::read_file(paste0(DIR, "/JAGS/Weibull_l23.txt")) else if(isFALSE(mm)&isFALSE(hm)) readr::read_file(paste0(DIR, "/JAGS/Weibull_l2.txt"))
  } else if(family=="Cox") {
    modelstring <- if(isTRUE(mm)&isTRUE(hm)) readr::read_file(paste0(DIR, "/JAGS/Cox_l123.txt")) else if(isTRUE(mm)&isFALSE(hm)) readr::read_file(paste0(DIR, "/JAGS/Cox_l12.txt")) else if(isFALSE(mm)&isTRUE(hm)) readr::read_file(paste0(DIR, "/JAGS/Cox_l23.txt")) else if(isFALSE(mm)&isFALSE(hm)) readr::read_file(paste0(DIR, "/JAGS/Cox_l2.txt"))
  }
  
  win2unix <- function(str) {
    gsub("\r\n", "\n", str, fixed = TRUE)
  }
  
  modelstring <- win2unix(modelstring)
  
  # Add predicted values of the dependent variable?
  if(monitor) {
    if(family=="Gaussian") modelstring <- stringr::str_replace(modelstring, "(Y\\[j\\] \\~) (dnorm\\(mu\\[j\\], tau.l2\\))", "\\1 \\2\n    pred[j] ~ \\2")
    if(family=="Weibull" & length(lhs)==2) modelstring <- stringr::str_replace(modelstring, "(t\\[j\\] \\~) (dweib\\(shape, lambda\\[j\\]\\))", "\\1 \\2\n\n    pred[j] ~ \\2")
    if(family=="Weibull" & length(lhs)==3) modelstring <- stringr::str_replace(modelstring, "(t\\[j\\] \\~) (dweib\\(shape, lambda\\[j\\]\\))", "\\1 \\2\n\n    ones[j] ~ dinterval(pred[j], ct.lbub[j,])\n    pred[j] ~ \\2")
    # 2do: predictions for the Cox model
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
    
    if(mmwfunction != "uw[i] <- 1/X.w[i,1]") modelstring <- stringr::str_replace(modelstring, fixed("uw[i] <- 1/X.w[i,1]"), mmwfunction)
    if(mmwconstraint == 2) modelstring <- stringr::str_replace(modelstring, fixed("w[i] <- uw[i] / sum(uw[l1i1.l1[i]:l1i2.l1[i]])"), "w[i] <- uw[i] * n.l2/sum(uw[]) # rescale to sum up to n.l1 overall")
    if(stringr::str_detect(mmwfunction, "b.w\\[.\\]")) { 
      modelstring <- stringr::str_replace(modelstring, fixed("b.w # placeholder\n"), paste0(mmwcoefstring, collapse = "")) 
    } else {
      modelstring <- stringr::str_replace(modelstring, fixed("b.w # placeholder\n"), "") 
    }
    
  }
  
  # Priors
  if(!is.null(priors)) {
    
    # Change parameters
    params <- unlist(sapply(names(priors), function(x) { 
      if(x %in% c("b.l1", "b.l2", "b.l3")) {
        paste0(x, "[x] ~ ") 
      } else if(x %in% c("b.w")) {
        if(isFALSE(mm)) stop("Priors cannot be set for b.w if no mm() construct is specified.")
        paste0(sapply(mmwcoefstring, function(y) stringr::str_extract(y, "b.w\\[.\\]"), USE.NAMES = F), " ~ ")
      } else { 
        paste0(x, " ~ ")
      }
    }, simplify = F), use.names = T)
    
    # Change priors 
    newpriors <- unlist(mapply(function(x, n) if(n == "b.w") rep(list(x), length(mmwcoefstring)) else list(x), priors, names(priors), SIMPLIFY=F), recursive=F)
  
    # Change priors in modelstring
    for(i in 1:length(params)) {
      modelstring <- stringr::str_replace(modelstring, fixed(as.vector(params[i])), paste0(params[i],  newpriors[names(params[i])][[1]], " # "))
    }
  }
  
  if(isTRUE(modelfile)) readr::write_file(modelstring, paste0(DIR, "/temp/modelstring.txt")) # save model to file
  else if (!isFALSE(modelfile) & length(modelfile)>0) modelstring <- readr::read_file(modelfile) # read model from file
  
  return(modelstring)
  
}
