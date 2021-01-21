# ================================================================================================ #
# Function createJagsVars
# ================================================================================================ #

createJagsVars <- function(family, data, level1, level2, level3, ids, vars, mm, l3, monitor, modelfile, seed, chains, inits) {
   
  # Unpack lists --------------------------------------------------------------------------------- #
  
  lhs <- vars$lhs
  
  mmwconstraint <- mm$mmwconstraint
  mmwar <- mm$mmwar
  lwvars <- vars$lwvars
  offsetvars <- vars$offsetvars
  
  l1vars <- level1$vars
  l1dat  <- level1$dat
  
  l2vars <- level2$vars
  l2dat  <- level2$dat
  
  l3vars <- level3$vars
  l3dat  <- level3$dat
  
  hm <- l3$hm
  l3name <- l3$l3name
  l3type <- l3$l3type
  showFE <- l3$showFE
  
  # IDs ------------------------------------------------------------------------------------------ #
  
  l1id <- l1dat %>% .$l1id # length = rows @ level1, must be sorted by l2id
  l2id <- l2dat %>% .$l2id # length = rows @ level2
  l3id <- if(hm) data %>% dplyr::group_by(l2id) %>% dplyr::filter(row_number()==1) %>% .$l3id else c() # length = rows @ level2
  
  l1n <- l2dat %>% .$l1n # number of l1-members per l2-unit
  
  l1i1 <- l2dat %>% .$l1i1 # first index of l1-members per l2-unit (@ level2)
  l1i1.l1 <- rep(l1i1, l1n) # first index of l1-members per l2-unit (@ level1)
  
  l1i2 <- l2dat %>% .$l1i2 # last index of l1-members per l2-unit (@ level2)
  l1i2.l1 <- rep(l1i2, l1n) # last index of l1-members per l2-unit (@ level1)
  
  # Xs ------------------------------------------------------------------------------------------- #
  
  X.l1 <- if(length(l1vars)>0) as.matrix(l1dat %>% dplyr::select(!!l1vars)) else c()
  X.l2 <- if(length(l2vars)>0) as.matrix(l2dat %>% dplyr::select(!!l2vars)) else c()
  X.l3 <- if(length(l3vars)>0) as.matrix(l3dat %>% dplyr::select(!!l3vars)) else c()
  X.w  <- if(length(lwvars)>0) as.matrix(data %>% dplyr::select(!!lwvars)) else c()
  
  # Ns ------------------------------------------------------------------------------------------- #
  
  n.l1  <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2  <- length(l2id)
  n.l3  <- length(unique(l3id))
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  n.Xl3 <- dim(X.l3)[2]
  n.Xw  <- dim(X.w)[2]
  
  n.GPN  <- l1dat %>% group_by(l1id) %>% count() %>% .$n %>% as.numeric() %>% max() # max number of gov participations
  n.GPNi <- l1dat %>% arrange(l1id) %>% group_by(l1id) %>% count() %>% .$n %>% as.numeric() # number of gov partipations per party, sorted l1id=1,2,3,...
  n.GPn  <- l1dat %>% group_by(l1id) %>% dplyr::mutate(n=row_number()) %>% .$n %>% as.numeric() # participation index, sorted l1id=2,6,2,...
  
  # ---------------------------------------------------------------------------------------------- #
  # Model-specific outcomes
  # ---------------------------------------------------------------------------------------------- #
  
  if(family=="Gaussian") {
    
    Y <- l2dat %>% dplyr::rename(Y = !!lhs) %>% .$Y

    # for return container
    Ys <- list("Y"=Y)
    sigma.l2 <- "sigma.l2"
    
  } else if(family=="Weibull") {
    
    t <- l2dat %>% dplyr::rename(t=lhs[1], ev=lhs[2]) %>% dplyr::mutate(t=case_when(ev==0 ~ NA_real_, TRUE ~ t)) %>% .$t
    t.cen <- l2dat %>% dplyr::rename(t=lhs[1], ev=lhs[2]) %>% dplyr::mutate(t.cens = t + ev) %>% .$t.cens
    event <- l2dat %>% dplyr::rename(ev=lhs[2]) %>% .$ev
    censored <- 1-event
    
    # for return container
    Ys <- list("t"=t, "t.cen"=t.cen, "event"=event, "censored"=censored)
    sigma.l2 <- "sigma.l2"
    
  } else if(family=="Cox") {
    
    t <- l2dat %>% dplyr::rename(t=lhs[1], ev=lhs[2]) %>% .$t
    t.unique <- append(sort(unique(t)), max(t)+1) 
    n.tu <- length(t.unique)-1
    event <- l2dat %>% dplyr::rename(ev=lhs[2]) %>% .$ev
    
    # Counting process data
    Y <- matrix(data = NA, nrow = n.l2, ncol = n.tu)
    dN <- matrix(data = NA, nrow = n.l2, ncol = n.tu)
    for(i in 1:n.l2) {
      for(j in 1:n.tu) {
        Y[i,j] <- as.numeric(t[i] - t.unique[j] + 1e-05 >= 0) # Risk set: Y[i,j] = 1 if t[j] >= t.unique[r] 
        dN[i,j] <- Y[i,j] * event[i] * as.numeric(t.unique[j+1] - t[i] >= 1e-05)  # Number of failures in each time interval: dN[i,j] = 1 if t in [ t.unique[j], t.unique[j+1] )
      }
    }
    
    # for return container
    Ys <- list("Y"=Y, "dN"=dN, "t"=t, "t.unique"=t.unique, "event"=event, "c"=0.001, "d"=0.1, "n.tu"=n.tu)
    sigma.l2 <- c()
    
  }
  
  # ---------------------------------------------------------------------------------------------- #
  # Create jags.params and jags.data
  # ---------------------------------------------------------------------------------------------- #
  
  # Level 1 -------------------------------------------------------------------------------------- #
  
  l1.param <- c("sigma.l1"); if(monitor==T) l1.param <- append(l1.param, c("re.l1")) 
  l1.data  <- c("l1id", "l1i1", "l1i2", "n.l1", "n.ul1"); if(mmwar==T) l1.data <- append(l1.data, c("n.GPn", "n.GPNi"))

  if(!is.null(n.Xl1)) { l1.data <- append(l1.data, c("X.l1", "n.Xl1")); l1.param <- append(l1.param, c("b.l1", "ppp.b.l1")) }
  
  # Level 2 -------------------------------------------------------------------------------------- #
  
  l2.param <- sigma.l2
  l2.data  <- c("n.l2")
  if(!is.null(n.Xl2)) { l2.data <- append(l2.data, c("X.l2", "n.Xl2")); l2.param <- append(l2.param, c("b.l2", "ppp.b.l2")) }
  
  # Level 3 -------------------------------------------------------------------------------------- #
  
  l3.param <- if(hm & l3type=="RE" & monitor==T) c("sigma.l3", "re.l3") else if(hm & l3type=="RE" & monitor==F) c("sigma.l3") else c()
  l3.data  <- if(hm) c("l3id", "n.l3") else c()
  if(!is.null(n.Xl3)) l3.data <- append(l3.data, c("X.l3", "n.Xl3"))
  if(!is.null(n.Xl3) & (l3type=="RE" | (l3type=="FE" & showFE==T))) l3.param <- append(l3.param, c("b.l3", "ppp.b.l3")) 
  
  # Weight function ------------------------------------------------------------------------------ #
  
  lw.param <- if(length(lwvars)>length(offsetvars) & monitor == T) c("b.w", "ppp.b.w", "w") else if(length(lwvars)>length(offsetvars) & monitor == F) c("b.w", "ppp.b.w") else if(length(lwvars)>0 & monitor==T) c("w") else c()
  lw.data  <- c("X.w")
  if(mmwconstraint==1) lw.data <- append(lw.data, c("l1i1.l1", "l1i2.l1"))
  
  jags.params <- c(l1.param, l2.param, l3.param, lw.param)
  jags.data   <- c(l1.data, l2.data, l3.data, lw.data)
  
  # Seed ----------------------------------------------------------------------------------------- #
  
  if(is.null(seed)) seed <- round(runif(1, 0, 1000))
  
  # ---------------------------------------------------------------------------------------------- #
  # Other model-specifics 
  # ---------------------------------------------------------------------------------------------- #
  
  if(family=="Gaussian") {
    
    # Dependent variable
    jags.data <- append(jags.data, "Y")
    
    # Initial values
    if(is.null(inits)) jags.inits <- list(".RNG.seed" = seed, tau.l1=1, tau.l2=1)
    
    if(hm & l3type=="RE") jags.inits <- lapply(jags.inits, FUN=function(x) { append(x, list(tau.l3=1.0)) }) 

  } else if(family=="Weibull") {
    
    # Dependent variable
    jags.data   <- append(jags.data, c("t", "t.cen", "censored"))
    
    # Initial values
    t.init <- t
    t.init[] <- NA
    t.init[censored==1] <- t.cen[censored==1] + 1 
    
    if(is.null(inits)) jags.inits <- list(".RNG.seed" = seed, t=t.init, shape=1.3, tau.l1=0.5)
    
  } else if(family=="Cox") {
    
    # Dependent variable
    jags.data   <- append(jags.data, c("Y", "dN", "t.unique", "n.tu", "c", "d"))
    
    # Initial values
    if(is.null(inits)) jags.inits <- list(".RNG.seed" = seed, dL0 = rep(1.0, n.tu), tau.l1=0.5)

  }
  
  # Repeat inits n.chains times
  jags.inits <- lapply(1:chains, function(x) { jags.inits[".RNG.seed"] <- seed+x; jags.inits } )
  
  # Read model in if provided
  if(is.character(modelfile)) modelstring <- readr::read_file(modelfile) 
  
  # Collect return ------------------------------------------------------------------------------- #

  return(list("ids"= list("l1id"=l1id, "l2id"=l2id, "l3id"=l3id, "l1i1"=l1i1, "l1i1.l1"=l1i1.l1, "l1i2"=l1i2, "l1i2.l1"=l1i2.l1),
              "Ns" = list("n.l1"=n.l1, "l1n"=l1n, "n.ul1"=n.ul1, "n.l2"=n.l2, "n.l3"=n.l3, "n.Xl1"=n.Xl1, "n.Xl2"=n.Xl2, "n.Xl3"=n.Xl3, "n.Xw"=n.Xw, "n.GPN"=n.GPN, "n.GPNi"=n.GPNi, "n.GPn"=n.GPn),
              "Xs" = list("X.l1"=X.l1, "X.l2"=X.l2, "X.l3"=X.l3, "X.w"=X.w),
              "Ys" = Ys,
              "jags.params" = jags.params, 
              "jags.inits"  = jags.inits,
              "jags.data"   = jags.data)) 
  
}
