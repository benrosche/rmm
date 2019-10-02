#' @title Fit Bayesian micro-macro regression models using JAGS
#'
#' @description Fit Bayesian micro-macro regression models using JAGS
#'
#' @param formula A symbolic description of the model in form of an R formula. More details below.
#' @param family Character vector, either "Gaussian" or "Weibull"
#' @param priors A list with parameter or variable names as tags and their prior specification as values. More details below.
#' @param iter Total number of iterations
#' @param burnin Number of iterations that will be discarded 
#' @param chain Number of chains
#' @param seed A random number
#' @param data The datasets must have level 1 as unit of analysis. More details below.
#'
#' @return JAGS output
#'
#' @examples rmm(Y ~ 1 + X1 + mm(id(l1id, l2id), mmc(X2 + X3), mmw(w ~ 1/N, constraint=1)),
#'     family="Gaussian",
#'     priors=list(var1="dnorm(0,0.0001)", tau.l1="dscaled.gamma(25, 1)"),
#'     data=data)
#'
#' @export rmm
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#' @details 
#' 
#' \bold{General formula structure}  
#' 
#' \code{Y ~ 1 + X.L2 + mm(id(l1id, l2id), mmc(X.L1), mmw(w ~ 1 / X.W, constraint=1))} 
#' \itemize{
#'   \item Dependent variable: Y 
#'   \item Vector of level-2 predictors: X.L2, being something like 1 + X1 + ... + XN 
#'   \item Micro-macro object: mm()
#' }
#' \bold{Currently supported dependent variables / link functions}  
#' 
#' \itemize{
#'   \item Gaussian continuous variable \bold{Y}
#'   \item Weibull survival time: \bold{Surv(survivaltime, event)}
#' }
#' 
#' \bold{Vector of level-2 predictors}
#' 
#' An intercept can be added by including a `1` at the beginning. Interaction terms have to be included as separate variables. 
#' Currently no support for nonlinear relationships.
#' 
#' \bold{Micro-macro object}
#' 
#' \itemize{
#'   \item \code{id()} to indicate level-1 and level-2 ids
#'   \item \code{mmc()} to specify level-1 predictors. No intercept allowed. Interaction terms have to be included as separate variables. Currently no support for nonlinear relationships.
#'   \item \code{mmw()} to specify the weight function (i.e. aggregation). Function can be nonlinear and contain variables. For instance, \code{w ~ 1/N} specifies mean-aggregation (with N being a variables that indicates the number of level-1 units per level-2 entity). If no mmw() is specified, w ~ 1/N is assumed. Moreover, two different identification restrictions are provided:
#'   \itemize{
#'     \item \code{mmw(w ~ 1/N, constraint=1)}: \code{constraint=1} restricts the weights to sum to 1 for each level-2 entity. 
#'     \item \code{mmw(w ~ 1/N, constraint=2)}: \code{constraint=2} restricts the weights to sum to the total number of level-2 entities over the whole dataset, allowing some level-2 entities to have weights smaller/larger than 1.
#'   }
#' }

rmm <- function(formula, family="Gaussian", priors, iter=1000, burnin=100, chains=3, seed=NULL, data=NULL) {

  stopifnot(!is.null(data))
  #devtools::load_all()
  
  # 1. Dissect formula ========================================================================== #
  
  getTerms <- function(e, x = list()) {
    if (is.symbol(e)) c(e, x)
    else if (identical(e[[1]], as.name("+"))) Recall(e[[2]], c(e[[3]], x))
    else c(e, x)
  }
  
  # Left-hand side
  Y <- all.vars(formula[[2]])
  
  # Right-hand side
  rhs <- gsub("[[:space:]]", "", getTerms(formula[3][[1]]))
  
  # mm()
  mmi <- which(startsWith(rhs, "mm("))
  mmstring <- rhs[mmi]
  
  # id()
  ids <- el(stringr::str_split(gsub(".*id\\((.*?)\\).*", "\\1", mmstring), ","))
  
  # mmc()
  l1vars <- el(stringr::str_split(gsub(".*mmc\\((.*?)\\).*", "\\1", mmstring), "\\+"))

  # mmlw(), function
  mmwstring <- el(stringr::str_split(gsub(".*mmw\\((.*?\\))\\).*", "\\1", mmstring), ","))
  mmwfunction <- mmwstring[1]
  lwvars <- all.vars(as.formula(stringr::str_remove(mmwfunction, fixed("b*",))))[-1]
  
  for(i in 1:length(lwvars)) {
    mmwfunction <- stringr::str_replace(mmwfunction, fixed(lwvars[i]), paste0("X.w[i,", i, "]"))
    mmwfunction <- stringr::str_replace(mmwfunction, fixed(paste0("b*X.w[i,", i, "]")), paste0("b.w[", i, "]*X.w[i,", i, "]"))
  }
  mmwfunction <- stringr::str_replace(mmwfunction, "w~", "uw[i] <- ")
  
  # mmlw(), constraint
  mmwconstraint <- ifelse(length(mmwstring)>1, as.numeric(gsub("\\D", "", mmwstring[2])), 1)

  if(mmwconstraint == 1) {
    mmwconstraint <- "w[i] <- uw[i] / sum(uw[l1i1.l1[i]:l1i2.l1[i]])" # Constraint 1
  } else if (mmwconstraint == 2) {
    mmwconstraint <- "w[i] <- uw[i] * n.l2/sum(uw[])" # Constraint 2
  } else {
    stop("Constraint must be either 1 or 2")
  }
  
  # level 2
  l2vars <- sub("^1$", "Intercept", rhs[-mmi])
  
  # 3. Edit modelstring ========================================================================= #
  
  # Load modelstring
  if(family=="Gaussian") {
    modelstring <- Model_Gaussian
  } else if(family=="Weibull") {
    modelstring <- Model_Weibull
  }
  
  # Priors
  vars <- names(priors)
  prior <- priors
  
  for(i in 1:length(vars)) {
    if(vars[i] %in% l1vars) {
      body(modelstring)[[5]][[4]][[2]] <- el(parse(text = paste0("b.l1[b] ~ ", prior[i])))
    } else if(vars[i] %in% l2vars) {
      body(modelstring)[[6]][[4]][[2]] <- el(parse(text = paste0("b.l2[b] ~ ", prior[i])))
    } else if(vars[i] %in% lwvars) {
      body(modelstring)[[7]][[4]][[2]] <- el(parse(text = paste0("b.lw[b] ~ ", prior[i])))
    } else if(vars[i] == "tau.l1") {
      body(modelstring)[[8]] <- el(parse(text = paste0("tau.l1 ~ ", prior[i])))
    } else if(vars[i] == "sigma.l2") {
      body(modelstring)[[10]]  <- el(parse(text = paste0("sigma.l2 ~ ", prior[i])))
    }
  }
  
  # Weight function
  body(modelstring)[[3]][[4]][[3]] <- el(parse(text = mmwfunction))
  body(modelstring)[[3]][[4]][[4]] <- el(parse(text = mmwconstraint))

  # 2. Disentangle level1 and level2 data ======================================================= #
  
  data <- data %>% dplyr::rename(l1id = !!ids[1], l2id = !!ids[2]) %>% ungroup() 
  
  # Level 1
  level1 <-
    data %>% 
    arrange(l2id) %>%
    dplyr::select(l1id, l2id, !!l1vars) %>%
    ungroup()
  
  # Level 2
  level2 <- 
    data %>% 
    arrange(l2id) %>%
    group_by(l2id) %>%
    add_count(name="l1n") %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    arrange(l2id) %>%
    mutate(l1i2=cumsum(l1n), l1i1=lag(l1i2)+1) %>%
    mutate(l1i1 = ifelse(row_number()==1, 1, l1i1)) %>%
    mutate(Intercept = 1) %>%
    dplyr::select(l1i1, l1i2, l2id, l1n, !!Y, !!l2vars) %>%
    ungroup()

  # 3. Prepare data for JAGS ==================================================================== #
  
  # IDs
  l1id <- level1 %>% .$l1id 
  l2id <- level2 %>% .$l2id 
  
  l1n <- level2 %>% .$l1n # number of l1-members per l2-unit
  
  l1i1 <- level2 %>% .$l1i1 # first index of l1-members per l2-unit (@ level2)
  l1i1.l1 <- rep(l1i1, l1n) # first index of l1-members per l2-unit (@ level1)
  
  l1i2 <- level2 %>% .$l1i2 # last index of l1-members per l2-unit (@ level2)
  l1i2.l1 <- rep(l1i2, l1n) # last index of l1-members per l2-unit (@ level1)
  
  # Y
  if(family=="Gaussian") {
    Y <- level2 %>% dplyr::rename(Y = !!Y) %>% .$Y
  } else if(family=="Weibull") {
    t <- level2 %>% dplyr::rename(t=Y[1], ev=Y[2]) %>% mutate(t=case_when(ev==0 ~ NA_real_, TRUE ~ t)) %>% .$t
    event <- level2 %>% dplyr::rename(ev=Y[2]) %>% .$ev
    t.cen <- level2 %>% dplyr::rename(t=Y[1], ev=Y[2]) %>% mutate(t.cens = t + ev) %>% .$t.cens
    censored <- 1-event
  }

  # Xs
  X.l1 <- as.matrix(level1 %>% dplyr::select(-l1id, -l2id))
  X.l2 <- as.matrix(level2 %>% select(-Y, -l1i1, -l1i2, -l2id, -l1n))
  X.w <- as.matrix(data %>% select(!!lwvars)) 
  
  # Ns
  n.l1 <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2 <- length(l2id)
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  n.Xlw <- dim(X.w)[2]
  
  if(family=="Gaussian") {
    
    # Monitored parameters
    b.w <- c()
    if(str_count(mmwfunction, fixed("b*"))>0) b.w <- "b.w"
    jags.params <- c("b.l1", "b.l2", b.w, "ppp.bl1", "ppp.bl2", "sigma.l1", "sigma.l2")
    
    # Data
    jags.data   <- c("l1id", "l1i1", "l1i1.l1", "l1i2", "l1i2.l1", "Y", "X.l1", "X.l2", "X.w", "l1n", "n.l1", "n.ul1", "n.l2", "n.Xl1", "n.Xl2", "n.Xlw") 
    
    # Initial values
    if(is.null(seed)) seed <- runif(1, 0,1000)
    jags.inits <- list( 
      list(".RNG.seed" = seed),
      list(".RNG.seed" = seed+1),
      list(".RNG.seed" = seed+2))
    
  } else if(family=="Weibull") {
    
    # Monitored parameters
    b.w <- c()
    if(str_count(mmwfunction, fixed("b*"))>0) b.w <- "b.w"
    jags.params <- c("b.l1", "b.l2", b.w, "ppp.bl1", "ppp.bl2", "sigma.l1", "sigma.l2")
    
    # Data
    jags.data   <- c("l1id", "l1i1", "l1i1.l1", "l1i2", "l1i2.l1", "t", "t.cen", "censored", "X.l1", "X.l2", "X.w", "l1n", "n.l1", "n.ul1", "n.l2", "n.Xl1", "n.Xl2", "n.Xlw") 
    
    # Initial values
    t.0 <- t
    t.0[censored==0] <- NA
    t.0[censored==1] <- t.cen[censored==1] + 1 
    
    jags.inits <- list( 
      list(".RNG.seed" = seed,   t=t.0, sigma.l2=1.0),
      list(".RNG.seed" = seed+1, t=t.0, sigma.l2=1.1),
      list(".RNG.seed" = seed+2, t=t.0, sigma.l2=1.3))
  }
  
  jags(data=jags.data, inits = jags.inits, parameters.to.save = jags.params, n.chains = 3, n.iter = iter, n.burnin = burnin, model.file = modelstring)

}
