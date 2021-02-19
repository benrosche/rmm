# ================================================================================================ #
# Function createJagsVars
# ================================================================================================ #

createJagsVars <- function(data, family, level1, level2, level3, weightf, ids, l1, l3, monitor, modelfile, chains, inits) {
   
  # Unpack lists --------------------------------------------------------------------------------- #
  
  lhs <- level2$lhs
  
  mm <- l1$mm
  mmwconstraint <- l1$mmwconstraint
  mmwar <- l1$mmwar
  
  hm <- l3$hm
  l3name <- l3$l3name
  l3type <- l3$l3type
  showFE <- l3$showFE
  
  l1vars <- level1$vars
  l1dat  <- level1$dat
  
  l2vars <- level2$vars
  l2dat  <- level2$dat
  
  l3vars <- level3$vars
  l3dat  <- level3$dat
  
  lwvars <- weightf$vars
  offsetvars <- weightf$offsetvars
  lwdat <- weightf$dat
  
  # ---------------------------------------------------------------------------------------------- #
  # Create IDs, Xs, and Ns
  # ---------------------------------------------------------------------------------------------- #

  # IDs ------------------------------------------------------------------------------------------ #
  
  l1id <- if(mm) l1dat %>% .$l1id else c() # length = rows @ level1, must be sorted by l2id
  l2id <- l2dat %>% .$l2id # length = rows @ level2
  l3id <- if(hm) l2dat %>% .$l3id else c() # length = rows @ level2
  
  l1n <- if(mm) l2dat %>% .$l1n else c() # number of l1-members per l2-unit
  
  l1i1 <- if(mm) l2dat %>% .$l1i1 else c() # first index of l1-members per l2-unit (@ level2)
  l1i1.l1 <- if(mm) rep(l1i1, l1n) else c() # first index of l1-members per l2-unit (@ level1)
  
  l1i2 <- if(mm) l2dat %>% .$l1i2 else c() # last index of l1-members per l2-unit (@ level2)
  l1i2.l1 <- if(mm) rep(l1i2, l1n) else c() # last index of l1-members per l2-unit (@ level1)
  
  # Xs ------------------------------------------------------------------------------------------- #
  
  X.l1 <- if(length(l1vars)>0) as.matrix(l1dat %>% dplyr::select(!!l1vars)) else c()
  X.l2 <- if(length(l2vars)>0) as.matrix(l2dat %>% dplyr::select(!!l2vars)) else c()
  X.l3 <- if(length(l3vars)>0) as.matrix(l3dat %>% dplyr::select(!!l3vars)) else c()
  X.w  <- if(length(lwvars)>0) as.matrix(lwdat %>% dplyr::select(!!lwvars)) else c()
  
  # Ns ------------------------------------------------------------------------------------------- #
  
  n.l1  <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2  <- length(l2id)
  n.l3  <- length(unique(l3id))
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  n.Xl3 <- dim(X.l3)[2]
  n.Xw  <- dim(X.w)[2]
  
  n.GPN  <- if(mm) l1dat %>% group_by(l1id) %>% count() %>% .$n %>% as.numeric() %>% max() else c() # max number of gov participations
  n.GPNi <- if(mm) l1dat %>% arrange(l1id) %>% group_by(l1id) %>% count() %>% .$n %>% as.numeric() else c() # number of gov partipations per party, sorted l1id=1,2,3,...
  n.GPn  <- if(mm) l1dat %>% group_by(l1id) %>% dplyr::mutate(n=row_number()) %>% .$n %>% as.numeric() else c() # participation index, sorted l1id=2,6,2,...
  
  # ---------------------------------------------------------------------------------------------- #
  # Create jags.params and jags.data
  # ---------------------------------------------------------------------------------------------- #
  
  # Level 1 -------------------------------------------------------------------------------------- #
  
  l1.param <- if(mm) c("sigma.l1") else c()
  l1.data  <- if(mm) c("l1id", "l1i1", "l1i2", "n.l1", "n.ul1") else c()
  
  if(mm & monitor) l1.param <- append(l1.param, c("re.l1"))
  if(mm & mmwar)   l1.data  <- append(l1.data, c("n.GPn", "n.GPNi"))
  if(mm & !is.null(n.Xl1)) { l1.data <- append(l1.data, c("X.l1", "n.Xl1")); l1.param <- append(l1.param, c("b.l1", "ppp.b.l1")) }
  
  # Level 2 -------------------------------------------------------------------------------------- #
  
  pred <- if(monitor) "pred" else NULL # predicted values of the dependent variable
  sigma.l2 <- if(!family=="Cox") "sigma.l2" else c() # Cox model does not have a variance term at level 2
  
  l2.param <- c(pred, sigma.l2)
  l2.data  <- c("n.l2")
  if(!is.null(n.Xl2)) { l2.data <- append(l2.data, c("X.l2", "n.Xl2")); l2.param <- append(l2.param, c("b.l2", "ppp.b.l2")) }
  
  # Level 3 -------------------------------------------------------------------------------------- #
  
  l3.param <- if(hm & l3type=="RE" & isTRUE(monitor)) c("sigma.l3", "re.l3") else if(hm & l3type=="RE" & isFALSE(monitor)) c("sigma.l3") else c()
  l3.data  <- if(hm) c("l3id", "n.l3") else c()
  if(!is.null(n.Xl3)) l3.data <- append(l3.data, c("X.l3", "n.Xl3"))
  if(!is.null(n.Xl3) & (l3type=="RE" | (l3type=="FE" & showFE==T))) l3.param <- append(l3.param, c("b.l3", "ppp.b.l3")) 
  
  # Weight function ------------------------------------------------------------------------------ #
  
  lw.param <- if(length(lwvars)>length(offsetvars) & monitor == T) c("b.w", "ppp.b.w", "w") else if(length(lwvars)>length(offsetvars) & monitor == F) c("b.w", "ppp.b.w") else if(length(lwvars)>0 & monitor==T) c("w") else c()
  lw.data  <- if(mm) c("X.w") else c()
  if(mmwconstraint==1) lw.data <- append(lw.data, c("l1i1.l1", "l1i2.l1"))
  
  # Collect terms -------------------------------------------------------------------------------- #
  
  jags.params <- c(l1.param, l2.param, l3.param, lw.param)
  jags.data   <- c(l1.data, l2.data, l3.data, lw.data)
  
  # ---------------------------------------------------------------------------------------------- #
  # Add model-specific params and data to jags.params and jags.data
  # ---------------------------------------------------------------------------------------------- #
  
  if(family=="Gaussian") {
    
    # Dependent variable
    Y <- l2dat %>% dplyr::rename(Y = !!lhs) %>% .$Y
    
    # jags.data 
    jags.data <- append(jags.data, "Y")
    
    # jags.inits
    jags.inits <- list() 
    
    # for return
    Ys <- list("Y"=Y)
    
  } else if(family=="Weibull") {
    
    # Survival time and censoring bounds
    t <- l2dat %>% dplyr::rename(t=!!lhs[1], ev=!!lhs[2]) %>% dplyr::mutate(t=case_when(ev==0 ~ NA_real_, TRUE ~ t)) %>% .$t # NA if no event / censored 
    ct.lb <- l2dat %>% dplyr::rename(t=!!lhs[1], ev=!!lhs[2]) %>% dplyr::mutate(ct.lb = t + ev) %>% .$ct.lb # ct.lb > t if ev
    ct.lbub <- if(length(lhs)==3) l2dat %>% dplyr::mutate(ct.lb=0) %>% dplyr::rename(ct.ub=!!lhs[3]) %>% dplyr::select(ct.lb, ct.ub) %>% as.matrix() else c() # lower and upper limit for predictions (if upper limit is provided)
    
    # Event and censoring status
    event <- l2dat %>% dplyr::rename(ev=!!lhs[2]) %>% .$ev
    censored <- 1-event
    if(monitor) ones <- rep(1, n.l2) 
    
    # jags.data
    jags.data  <- append(jags.data, c("t", "ct.lb", "censored")) 
    if(monitor & length(lhs)==3) jags.data <- append(jags.data, c("ct.lbub", "ones")) # predictions with upper bounds
    
    # jags.params
    jags.params <- append(jags.params, "shape")
    
    # jags.inits
    t.init <- t
    t.init[] <- NA
    t.init[censored==1] <- ct.lb[censored==1] + 1 
    
    jags.inits <- list(t=t.init, shape=1) 
    
    # for return
    Ys <- list("t"=t, "ct.lb"=ct.lb, "ct.lbub"=ct.lbub, "event"=event, "censored"=censored, "ones"=ones)
    
  } else if(family=="Cox") {
    
    t <- l2dat %>% dplyr::rename(t=!!lhs[1], ev=!!lhs[2]) %>% .$t
    t.unique <- append(sort(unique(t)), max(t)+1) 
    n.tu <- length(t.unique)-1
    event <- l2dat %>% dplyr::rename(ev=!!lhs[2]) %>% .$ev
    
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
    jags.data   <- append(jags.data, c("Y", "dN", "t.unique", "n.tu", "c", "d"))
    
    # jags.inits
    jags.inits <- list(dL0 = rep(1.0, n.tu)) 
    
    # for return
    Ys <- list("Y"=Y, "dN"=dN, "t"=t, "t.unique"=t.unique, "event"=event, "c"=0.001, "d"=0.1, "n.tu"=n.tu)
    
  }
  
  # ---------------------------------------------------------------------------------------------- #
  # Modify jags.inits 
  # ---------------------------------------------------------------------------------------------- #
  
  # Add user-defined inits
  jags.inits <- append(jags.inits, inits)
  
  # Repeat inits n.chains times
  jags.inits <- lapply(1:chains, function(x) { jags.inits } )
  
  # ---------------------------------------------------------------------------------------------- #
  # Read model in if provided
  # ---------------------------------------------------------------------------------------------- #
  
  if(is.character(modelfile)) modelstring <- readr::read_file(modelfile) 
  
  # ---------------------------------------------------------------------------------------------- #
  # Collect return
  # ---------------------------------------------------------------------------------------------- #

  return(list("ids"= list("l1id"=l1id, "l2id"=l2id, "l3id"=l3id, "l1i1"=l1i1, "l1i1.l1"=l1i1.l1, "l1i2"=l1i2, "l1i2.l1"=l1i2.l1),
              "Ns" = list("n.l1"=n.l1, "l1n"=l1n, "n.ul1"=n.ul1, "n.l2"=n.l2, "n.l3"=n.l3, "n.Xl1"=n.Xl1, "n.Xl2"=n.Xl2, "n.Xl3"=n.Xl3, "n.Xw"=n.Xw, "n.GPN"=n.GPN, "n.GPNi"=n.GPNi, "n.GPn"=n.GPn),
              "Xs" = list("X.l1"=X.l1, "X.l2"=X.l2, "X.l3"=X.l3, "X.w"=X.w),
              "Ys" = Ys,
              "jags.params" = jags.params, 
              "jags.inits"  = jags.inits,
              "jags.data"   = jags.data)) 
  
}
