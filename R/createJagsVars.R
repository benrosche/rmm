# ================================================================================================ #
# Function createJagsVars
# ================================================================================================ #

createJagsVars <- function(data, family, level1, level2, level3, weight, ids, l1, l3, monitor, modelfile, chains, inits) {
   
  # Unpack lists --------------------------------------------------------------------------------- #
  
  lhs <- level2[["lhs"]]
  
  mm <- l1[["mm"]]
  mmwconstraint <- l1[["mmwconstraint"]]
  mmwar <- l1[["mmwar"]]
  
  hm <- l3[["hm"]]
  l3name <- l3[["l3name"]]
  l3type <- l3[["l3type"]]
  showFE <- l3[["showFE"]]
  
  l1vars <- level1[["vars"]]
  l1dat  <- level1[["dat"]]
  
  l2vars <- level2[["vars"]]
  l2dat  <- level2[["dat"]]
  
  l3vars <- level3[["vars"]]
  l3dat  <- level3[["dat"]]
  
  wvars <- weight[["vars"]]
  wvars_p <- weight[["vars_p"]]
  wparams <- weight[["params"]]
  wdat <- weight[["dat"]]
  
  # ---------------------------------------------------------------------------------------------- #
  # Create IDs, Xs, and Ns
  # ---------------------------------------------------------------------------------------------- #

  # IDs ------------------------------------------------------------------------------------------ #
  
  l1id <- if(mm) l1dat %>% pull(l1id) else c() # length = rows @ level1, must be sorted by l2id
  l2id <- l2dat %>% pull(l2id) # length = rows @ level2
  l3id <- if(hm) l2dat %>% pull(l3id) else c() # length = rows @ level2
  
  l1n <- if(mm) l2dat %>% pull(l1n) else c() # number of l1-members per l2-unit
  
  l1i1 <- if(mm) l2dat %>% pull(l1i1) else c() # first index of l1-members per l2-unit (@ level2)
  l1i1.l1 <- if(mm) rep(l1i1, l1n) else c() # first index of l1-members per l2-unit (@ level1)
  
  l1i2 <- if(mm) l2dat %>% pull(l1i2) else c() # last index of l1-members per l2-unit (@ level2)
  l1i2.l1 <- if(mm) rep(l1i2, l1n) else c() # last index of l1-members per l2-unit (@ level1)
  
  # Xs ------------------------------------------------------------------------------------------- #
  
  X.l1 <- if(length(l1vars)>0) l1dat %>% dplyr::select(all_of(l1vars)) %>% as.matrix() else c()
  X.l2 <- if(length(l2vars)>0) l2dat %>% dplyr::select(all_of(l2vars)) %>% as.matrix() else c()
  X.l3 <- if(length(l3vars)>0) l3dat %>% dplyr::select(all_of(l3vars)) %>% as.matrix() else c()
  X.w  <- if(length(wvars)>0)   wdat %>% dplyr::select(all_of(wvars))  %>% as.matrix() else c()
  
  # Ns ------------------------------------------------------------------------------------------- #
  
  n.l1  <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2  <- length(l2id)
  n.l3  <- length(unique(l3id))
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  n.Xl3 <- dim(X.l3)[2]
  n.Xw  <- dim(X.w)[2]
  
  n.GPN  <- if(mm) l1dat %>% dplyr::count(l1id) %>% pull(n) %>% max() else c() # max number of gov participations
  n.GPNi <- if(mm) l1dat %>% dplyr::count(l1id) %>% pull(n) else c() # number of gov partipations per party, sorted l1id=1,2,3,...
  n.GPn  <- if(mm) l1dat %>% group_by(l1id) %>% dplyr::mutate(n=row_number()) %>% pull(n) else c() # participation index, sorted l1id=2,6,2,...
  
  # ---------------------------------------------------------------------------------------------- #
  # Create jags.params and jags.data
  # ---------------------------------------------------------------------------------------------- #
  
  # Level 1 -------------------------------------------------------------------------------------- #
  
  l1.params <- c()
  l1.data <- c()
  
  if(mm) l1.params <- c(l1.params, "sigma.l1")
  if(mm) l1.data  <- c(l1.data, c("l1id", "l1i1", "l1i2", "n.l1", "n.ul1"))
  if(mm & !is.null(n.Xl1)) { l1.data <- c(l1.data, c("X.l1", "n.Xl1")); l1.params <- c(l1.params, c("b.l1")) }
  
  if(mm & monitor) l1.params <- c(l1.params, "re.l1")
  if(mm & mmwar)   l1.data  <- c(l1.data, c("n.GPn", "n.GPNi"))

  # Level 2 -------------------------------------------------------------------------------------- #
  
  l2.params <- c()
  l2.data <- c()
  
  if(family!="Cox") l2.params <- c(l2.params, "sigma.l2") # Cox model does not have a variance term at level 2
  l2.data <- c(l2.data, "n.l2")
  if(!is.null(n.Xl2)) { l2.params <- c(l2.params, c("b.l2")); l2.data <- c(l2.data, c("X.l2", "n.Xl2")) }
  
  if(monitor) l2.params <- c(l2.params, "pred") # predicted values of the dependent variable

  # Level 3 -------------------------------------------------------------------------------------- #
  
  l3.params <- c()
  l3.data <- c()
  
  if(hm & l3type=="RE" & isTRUE(monitor)) l3.params <- c(l3.params, "sigma.l3", "re.l3") 
  if(hm & l3type=="RE" & isFALSE(monitor)) l3.params <- c(l3.params, "sigma.l3") 
  if(hm) l3.data <- c(l3.data, "l3id", "n.l3")
 
  if(!is.null(n.Xl3)) l3.data <- c(l3.data, c("X.l3", "n.Xl3"))
  if(!is.null(n.Xl3) & (l3type=="RE" | (l3type=="FE" & showFE==T))) l3.params <- c(l3.params, c("b.l3")) 
  
  # Weight function ------------------------------------------------------------------------------ #
  
  w.params <- c()
  w.data <- c()
  
  if(length(wparams)>0 & monitor == T) w.params <- c(w.params, "b.w", "w")
  if(length(wparams)==0 & length(wvars)>0 & monitor==T) w.params <- c(w.params, "w") 
  if(length(wparams)>0 & monitor == F) w.params <- c(w.params, "b.w") 

  if(mm) w.data <- c(w.data, "X.w")
  if(mmwconstraint==T) w.data <- c(w.data, c("l1i1.l1", "l1i2.l1"))
  
  # Collect terms -------------------------------------------------------------------------------- #
  
  jags.params <- c(l1.params, l2.params, l3.params, w.params)
  jags.data   <- c(l1.data, l2.data, l3.data, w.data)
  
  # ---------------------------------------------------------------------------------------------- #
  # Add model-specific params and data to jags.params and jags.data
  # ---------------------------------------------------------------------------------------------- #
  
  if(family %in% c("Gaussian", "Binomial")) {
    
    # Dependent variable
    Y <- l2dat %>% dplyr::rename(Y = all_of(lhs)) %>% pull(Y)
    
    # jags.data
    jags.data <- c(jags.data, "Y")
    
    # jags.inits
    jags.inits <- list() 
    
    # for return
    Ys <- list("Y"=Y)
    
  } else if(family=="Weibull") {
    
    # Survival time and censoring bounds
    t <- l2dat %>% dplyr::rename(t=all_of(lhs[1]), ev=all_of(lhs[2])) %>% dplyr::mutate(t=case_when(ev==0 ~ NA_real_, TRUE ~ t)) %>% pull(t) # NA if no event / censored 
    ct.lb <- l2dat %>% dplyr::rename(t=all_of(lhs[1]), ev=all_of(lhs[2])) %>% dplyr::mutate(ct.lb = t + ev) %>% pull(ct.lb) # ct.lb > t if ev
    
    # Event and censoring status
    event <- l2dat %>% dplyr::rename(ev=all_of(lhs[2])) %>% pull(ev)
    censored <- 1-event
    
    # jags.data
    jags.data  <- c(jags.data, c("t", "ct.lb", "censored")) 
    
    # jags.params
    jags.params <- c(jags.params, "shape")
    
    # jags.inits
    t.init <- t
    t.init[] <- NA
    t.init[censored==1] <- ct.lb[censored==1] + 1 
    
    jags.inits <- list(t=t.init, shape=1) 
    
    # for return
    Ys <- list("t"=t, "ct.lb"=ct.lb, "event"=event, "censored"=censored)
    
  } else if(family=="Cox") {
    
    t <- l2dat %>% dplyr::rename(t=all_of(lhs[1]), ev=all_of(lhs[2])) %>% pull(t)
    t.unique <- c(sort(unique(t)), max(t)+1) 
    n.tu <- length(t.unique)-1
    event <- l2dat %>% dplyr::rename(ev=all_of(lhs[2])) %>% pull(ev)
    
    # Counting process data
    Y <- matrix(data = NA, nrow = n.l2, ncol = n.tu)
    dN <- matrix(data = NA, nrow = n.l2, ncol = n.tu)
    for(i in 1:n.l2) {
      for(j in 1:n.tu) {
        Y[i,j] <- as.numeric(t[i] - t.unique[j] + 1e-05 >= 0) # Risk set: Y[i,j] = 1 if t[j] >= t.unique[r] 
        dN[i,j] <- Y[i,j] * event[i] * as.numeric(t.unique[j+1] - t[i] >= 1e-05)  # Number of failures in each time interval: dN[i,j] = 1 if t in [ t.unique[j], t.unique[j+1] )
      }
    }
    
    # jags.data
    jags.data   <- c(jags.data, c("Y", "dN", "t.unique", "n.tu", "c", "d"))
    
    # jags.inits
    jags.inits <- list(dL0 = rep(1.0, n.tu)) 
    
    # for return
    Ys <- list("Y"=Y, "dN"=dN, "t"=t, "t.unique"=t.unique, "event"=event, "c"=0.001, "d"=0.1, "n.tu"=n.tu)
    
  }
  
  # ---------------------------------------------------------------------------------------------- #
  # Modify jags.inits 
  # ---------------------------------------------------------------------------------------------- #
  
  # Add user-defined inits
  jags.inits <- c(jags.inits, inits)
  
  # Repeat inits n.chains times
  jags.inits <- lapply(1:chains, function(x) { jags.inits } )
  
  # ---------------------------------------------------------------------------------------------- #
  # Read model in if provided
  # ---------------------------------------------------------------------------------------------- #
  
  if(is.character(modelfile)) modelstring <- readr::read_file(modelfile) 
  
  # ---------------------------------------------------------------------------------------------- #
  # Collect return
  # ---------------------------------------------------------------------------------------------- #

  return(
    list(
      "ids"= list("l1id"=l1id, "l2id"=l2id, "l3id"=l3id, "l1i1"=l1i1, "l1i1.l1"=l1i1.l1, "l1i2"=l1i2, "l1i2.l1"=l1i2.l1),
      "Ns" = list("n.l1"=n.l1, "l1n"=l1n, "n.ul1"=n.ul1, "n.l2"=n.l2, "n.l3"=n.l3, "n.Xl1"=n.Xl1, "n.Xl2"=n.Xl2, "n.Xl3"=n.Xl3, "n.Xw"=n.Xw, "n.GPN"=n.GPN, "n.GPNi"=n.GPNi, "n.GPn"=n.GPn),
      "Xs" = list("X.l1"=X.l1, "X.l2"=X.l2, "X.l3"=X.l3, "X.w"=X.w),
      "Ys" = Ys,
      "jags.params" = jags.params,
      "jags.inits"  = jags.inits,
      "jags.data"   = jags.data
    )
  ) 
  
}
