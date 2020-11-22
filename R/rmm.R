#' @title Bayesian multiple membership multilevel models with endogenized weights using 'JAGS'
#'
#' @description The \strong{rmm} package provides an interface to fit Bayesian multiple membership 
#' multilevel models with endogenized weights using \href{http://mcmc-jags.sourceforge.net/}{JAGS}.
#' 
#' @details 
#' The main function of \strong{rmm} is \code{rmm}, which uses formula syntax to specify a multiple membership 
#' multilevel model with endogenized weights. Based on the supplied formulas, data, and additional information, 
#' it writes code to fit the model in JAGS. Subsequently, the JAGS output is processed to ease the
#' interpretation of the model results.
#'  
#' In order to fit the models, JAGS must be \href{https://sourceforge.net/projects/mcmc-jags/files/}{installed}.
#' 
#' The package \strong{rmm} estimates models with a complex nonstandard multilevel structure, known 
#' as multiple membership multilevel structure. The difference of this package to other packages and programs 
#' to estimate multiple membership multilevel models, such as \code{\link[brms]{brms}} or \href{http://www.bristol.ac.uk/cmm/software/mlwin/}{MLwiN}, 
#' is that \strong{rmm} allows to endogenize the membership weights with a weight function. In doing so, rmm allows 
#' to examine the process by which the effects of lower-level units aggregate to a higher level (micro-macro link).
#' 
#' Accessible introductions to multiple membership models are given by the report by Fielding and Goldstein (2006) 
#' and the book chapter by Beretvas (2010). More advanced treatments of multiple membership models are  
#' provided in the multilevel textbook by Goldstein (2011, Chapter 13), the book chapters on 
#' multiple membership models by Rasbash and Browne (2001, 2008), the paper by Browne et al. (2001), and the report by Leckie (2013).
#' 
#' \bold{General formula structure}  
#' 
#' \code{Y ~ 1 + mm(id(l1id, l2id), mmc(X.L1), mmw(w ~ 1 / offset(N), constraint=1)) + X.L2 + X.L3 + hm(id=l3id, name=l3name, type=FE, showFE=F)} 
#' \itemize{
#'   \item Dependent variable: \strong{Y}
#'   \item Multiple membership object: \strong{mm()} to analyze how the effects of level-1 predictors from multiple constituting members aggregate to level 2
#'   \item Level-2 predictors: \strong{X.L2}, being something like X1 + ... + XN 
#'   \item Level-3 predictors: \strong{X.L3}, being something like X1 + ... + XN 
#'   \item Hierarchical membership object: \strong{hm()} to recognize that level 2 units embedded in a third level
#' }
#' \bold{Currently supported dependent variables / link functions}  
#' 
#' \itemize{
#'   \item Gaussian continuous variable \code{Y}
#'   \item Binomial outcome for logistic regression \code{Y}
#'   \item Conditional logistic outcomes \code{???}
#'   \item Weibull survival time: \code{Surv(survivaltime, event)}
#'   \item Cox survival time: \code{Surv(survivaltime, event)}
#' }
#' 
#' \bold{Vector of level-2 predictors}
#' 
#' An intercept at the main level 2 is added whether or not a \code{1} is specified in the beginning. Interaction terms have to be included as separate variables. 
#' 
#' \bold{Multiple membership object mm()}
#' 
#' \itemize{
#'   \item \code{id()} to indicate level-1 and level-2 ids
#'   \item \code{mmc()} to specify level-1 predictors. No intercept allowed. Interaction terms have to be included as separate variables. 
#'   \item \code{mmw()} to specify the weight function (micro-macro link). The function can be nonlinear and contain variables but needs to be identifiable. 
#'         To give a few examples: \code{w ~ 1/offset(N)} specifies mean-aggregation, with N being a variables that indicates the number of level-1 units per level-2 entity.
#'         If no mmw() is specified, \code{w ~ 1/offset(N)} is assumed. In Rosche (2021), I propose to use \code{w ~ 1/offset(N)^exp(-(X.W))} as the general form for weight functions. 
#'         This function ensures that weights are limited by 0 and 1. Specifying variables as \code{offset(X.W)} will not estimate a parameter for this variable. 
#'         
#'         Two different identification restrictions are provided:
#'   \itemize{
#'     \item \code{mmw(w ~ 1/N, constraint=1)}: \code{constraint=1} restricts the weights to sum to 1 for each level-2 entity. (default)
#'     \item \code{mmw(w ~ 1/N, constraint=2)}: \code{constraint=2} restricts the weights to sum to the total number of level-2 entities over the whole dataset, allowing some level-2 entities to have weights smaller/larger than 1.
#'   }
#' }
#' 
#' \bold{Hierachical membership object hm()}
#' 
#' \itemize{
#'   \item \code{id=l3id} to indicate level-3 id
#'   \item \code{name=l3name} to specify value labels for level 3 units. 
#'   \item \code{type=RE} (default) or \code{type=FE} to choose between random- or fixed effect estimation. 
#'         If RE is chosen, level 3 predictors can be added. If FE is chosen, each level 3 unit has its own intercept and level 3 predictors are removed. 
#'         If \code{showFE=TRUE} the fixed effects are reported, otherwise omitted (default).
#' }
#' 
#' \bold{More details on the weight function}
#' 
#' Add variables as \code{offset(X)} to not estimate their effect (but assume Beta=1). 
#' The general weight function as proposed in Rosche (XXXX: XX) is implemented by \code{w ~ 1/offset(N)^exp(-(X.W))}
#' 
#' \bold{More details on adding a third level}
#' 
#' ....
#'
#' \bold{More details on constructing the data}
#' ...
#' 
#' @param formula A symbolic description of the model in form of an R formula. More details below.
#' @param family Character vector. Currently supported are "Gaussian", "Logit", "Condlogit", "Weibull", or "Cox".
#' @param priors A list with parameter or variable names as tags and their prior specification as values. More details below.
#' @param iter Total number of iterations.
#' @param burnin Number of iterations that will be discarded .
#' @param chains Number of chains.
#' @param seed A random number.
#' @param run A logical value (True or False) indicating whether JAGS should estimate the model.
#' @param modelfile Path to modelfile. If \code{NULL}, the modelfile is created and will be located at \code{'rmm/temp/modelstring.txt'}.
#' @param monitor A logical value (True or False). If \code{True}, weights and random effects are monitored and saved in the global environment.
#' @param hdi Numeric or False. If confidence level \code{x} is specified (default \code{x=0.95}), \code{mode} and \code{(x*100)%} HDI estimates are given. If \code{False} is specified, \code{mean} and \code{95% CI} are given.
#' @param r Numeric. Rounding value. Default is 3.
#' @param data Dataframe object. The dataset must have level 1 as unit of analysis. More details below.
#'
#' @return JAGS output. More details the output ...
#'
#' @examples data(coalgov)
#' rmm(Y ~ 1 + mm(id(l1id, l2id), mmc(X11 + X12), mmw(w ~ 1/offset(N)^exp(-(X1+X4)), constraint=1)) + X21 + X31 + hm(id=l3id, name=l3name, type=FE, showFE=F),
#'     family="Gaussian",
#'     priors=list(X11="dnorm(0,0.0001)", tau.l1="dscaled.gamma(25, 1)"),
#'     data=coalgov)
#'
#' @export rmm
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#' @references 
#' Rosche (XXXX): The multilevel structure of coalition government outcomes

rmm <- function(formula, family="Gaussian", priors=NULL, iter=1000, burnin=100, chains=3, seed=NULL, run=T, modelfile=NULL, monitor=F, hdi=0.95, r=3, data=NULL) {

  # formula <- sim_y ~ 1 + mwc + investiture + hetero + mm(id(pid, gid), mmc(ipd +fdep), mmw(w ~ 1/offset(n), c=1)) + hm(id="cid", type=RE); family <- "Weibull"; priors=NULL; iter=1000; burnin=100; chains <- 3; seed <- NULL; run <- T; modelfile <- NULL; monitor <- F; hdi=0.95; r=3; data <- coalgov;
  
  if(is.null(data)) stop("No data supplied.")
  DIR <- system.file(package = "rmm")
  #devtools::load_all()
  
  # ----------------------------------------------------------------------------------------------#
  # 1. Dissect formula 
  # --------------------------------------------------------------------------------------------- #
  
  getTerms <- function(e, x = list()) {
    if (is.symbol(e)) c(e, x)
    else if (identical(e[[1]], as.name("+"))) Recall(e[[2]], c(e[[3]], x))
    else c(e, x)
  }
  
  # Left-hand side
  lhs <- all.vars(formula[[2]])
  if(family=="Gaussian" & length(lhs)>1) stop("family=\"Gaussian\" takes only one variable on the left-hand side of the formula.")
  if(family=="Weibull" & length(lhs)!=2) stop("family=\"Weibull\" takes two variables on the left-hand side: 'Surv(survtime, event)'")
  
  # Right-hand side
  rhs <- gsub("[[:space:]]", "", getTerms(formula[3][[1]]))
  
  # mm() object ---------------------------------------------------------------------------------- #
  
  # mm()
  mmi <- which(startsWith(rhs, "mm("))
  mmstring <- rhs[mmi]
  
  if(length(mmstring)==0) stop("Multiple membership construct 'mm()'  missing. Please use the lme4 package to estimate conventional multilevel models.")
  
  # id()
  l12ids <- if(stringr::str_detect(mmstring, "id\\(.*\\)")) unlist(stringr::str_split(stringr::str_replace(mmstring, ".*id\\((.*?)\\).*", "\\1"), ",")) else stop("'id(l1id,l2id)' missing within 'mm()'")
    
  if(all(l12ids %in% names(data))==FALSE) stop("Either l1id or l2id could not be found in the dataset.")
  
  # mmc()
  l1vars <- if(stringr::str_detect(mmstring, "mmc\\(.*\\)")) unlist(stringr::str_split(stringr::str_replace(mmstring, ".*mmc\\((.*?)\\).*", "\\1"), "\\+")) else stop("l1vars 'mmc()' missing within 'mm()'")
  if(all(nzchar(l1vars)) == FALSE) l1vars <- c()
  if(all(l1vars %in% names(data))==FALSE) stop("l1vars within 'mmc()' could not be found in the dataset.")
  
  # mmlw(), function
  mmwstring <- if(stringr::str_detect(mmstring, "mmw\\(.*\\)")) unlist(stringr::str_split(stringr::str_replace(mmstring, ".*mmw\\((.*)\\)(.*)\\)", "\\1 \\2"), ",")) else stop("Weight function 'mmw()' missing within 'mm()'")
                                                                                              
  mmwfunction <- mmwstring[1] # weight function
  
  if(!startsWith(mmwfunction, "w")) stop("Weight function within 'mmw()' missing")
  
  lwvars <- tryCatch(all.vars(as.formula(mmwfunction))[-1], error=function(e) { print("Weight variables within 'mmw()' missing") }) # all variables but 'w'
  
  if(all(lwvars %in% names(data))==FALSE | length(lwvars)==0) stop("Weight variables within 'mmw()' could not be found in the dataset.")
  if(length(lwvars)<length(all.vars(as.formula(mmwfunction),unique=F)[-1])) stop("Weight function must not include the same variable multiple times.") # must be at this position
  
  offsetvars <- stringr::str_replace(stringr::str_extract_all(mmwfunction, "offset\\((.+?)\\)", simplify = T), "offset\\((.+?)\\)", "\\1")
  
  # Translate weight function into JAGS format
  mmwcoefstring <- c()
  for(i in 1:length(lwvars)) { 
    if(lwvars[i] %in% offsetvars) {
      mmwfunction <- stringr::str_replace(mmwfunction, fixed(paste0("offset(", lwvars[i], ")")), paste0("(X.w[i,", i, "])")) # without coefficient
    } else {
      mmwfunction <- stringr::str_replace(mmwfunction, fixed(lwvars[i]), paste0("(b.w[", i, "]*X.w[i,", i, "])")) # with coefficient
      mmwcoefstring  <- append(mmwcoefstring, paste0("b.w[", i, "] ~ dnorm(0,0.0001)\r\n  ppp.b.w[", i, "] <- step(b.w[", i, "])\r\n  "))
    }
  }
  mmwfunction <- stringr::str_replace(mmwfunction, "w~", "uw[i] <- ")
  
  # mmlw(), constraint
  if(length(mmwstring)>1) {
    mmwconstraint <- as.numeric(str_extract(paste0(mmwstring[-1], collapse = ", "), "(?!constraint\\=)[0-9]"))
    if(is.na(mmwconstraint) & length(lwvars)>length(offsetvars)) mmwconstraint <- 2
    if(is.na(mmwconstraint) & length(lwvars)<=length(offsetvars)) mmwconstraint <- 1
  } else {
    if(length(lwvars)>length(offsetvars)) mmwconstraint <- 2
    if(length(lwvars)<=length(offsetvars)) mmwconstraint <- 1
  }
  
  if(!mmwconstraint %in% c(1,2)) stop("Constraint must be either 1 or 2")
  
  # mmlw(), ar
  if(length(mmwstring)>1) {
    mmwar <- as.logical(str_extract(paste0(mmwstring[-1], collapse = ", "), "(?!ar\\=)(T|F|TRUE|FALSE)"))
    if(is.na(mmwar)) mmwar <- FALSE
  } else {
    mmwar <- FALSE
  }
  if(!mmwar %in% c(T,F)) stop("AR must be either True or False")
    
  # hm() object ---------------------------------------------------------------------------------- #
  
  hmi <- which(startsWith(rhs, "hm("))
  hmstring <- rhs[hmi]
  
  if(length(hmi)>0) {

    l3id    <- stringr::str_remove(stringr::str_extract(hmstring, "(?<=id\\=)(.*?)(?=[:punct:]?,|\\))"), "[:punct:]*") # level-3 id
    if(!l3id %in% names(data)) stop("'l3id' in 'hm()' could not be found in the dataset.")
    
    l3name  <- stringr::str_remove(stringr::str_extract(hmstring, "(?<=name\\=)(.*?)(?=[:punct:]?,|\\))"), "[:punct:]*") # level-3 names
    if(!is.na(l3name) & !l3name %in% names(data)) stop("'l3name' in 'hm()' is specified but could not be found in the dataset.")
    if(is.na(l3name)) l3name <- c()
    
    l3type  <- stringr::str_remove(stringr::str_extract(hmstring, "(?<=type\\=)(.*?)(?=[:punct:]?,|\\))"), "[:punct:]*") # FE|RE
    if(!is.na(l3type) & !l3type %in% c("RE","FE")) stop("'type' in 'hm()' must be either 'RE', 'FE', or left unspecified")
    if(is.na(l3type)) l3type <- c("RE")
    
    showFE  <- stringr::str_remove(stringr::str_extract(hmstring, "(?<=showFE\\=)(.*?)(?=[:punct:]?,|\\))"), "[:punct:]*") # Should FE be monitored?
    if(!is.na(showFE) & !showFE %in% c("T", "F", "TRUE", "FALSE")) stop("'show' in 'hm()' must be either 'TRUE', 'FALSE', or left unspecified")
    showFE <- if(is.na(showFE)) FALSE else as.logical(showFE)
    
  } else {
    l3id <- "one_vec" # an existing variable
  }
  
  # Level 2 + 3 vars------------------------------------------------------------------------------ #
  
  l23vars <- stringr::str_replace(rhs[c(-mmi, -hmi)], "^1$", "X0") # Intercept
  if(!"X0" %in% l23vars)  l23vars <- append(l23vars, "X0", after=0)
  
  # ---------------------------------------------------------------------------------------------- #
  # 2. Disentangle vars and data into l1-3
  # ---------------------------------------------------------------------------------------------- #
  
  data <- # rename and regroup ids and sort
    data %>% 
    dplyr::rename(l1id = !!l12ids[1], l2id = !!l12ids[2]) %>%
    dplyr::mutate(one_vec = 1) %>% # must be here
    dplyr::mutate(l3id = !!rlang::sym(l3id)) %>% 
    dplyr::group_by(l1id) %>% dplyr::mutate(l1id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(l2id) %>% dplyr::mutate(l2id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(l3id) %>% dplyr::mutate(l3id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::arrange(l2id, l1id) 
    
  # Level 1
  level1 <-
    data %>% 
    dplyr::arrange(l2id, l1id) %>% # important
    dplyr::select(l1id, l2id, !!l1vars) %>%
    dplyr::mutate_at(l1vars, function(x) { if(is.numeric(x) & dim(table(x))>2) x-mean(x) else x }) # center continuous vars 

  # Level 3
  if(length(hmi)>0) {
    
    l3vars <- 
      data %>% 
      dplyr::select(l3id, l23vars[-1]) %>%
      dplyr::group_by(l3id) %>%
      dplyr::mutate_at(l23vars[-1], ~var(., na.rm = TRUE)) %>% # select variables that do not vary within levels
      dplyr::select_if(~sum(.)==0) %>%  
      dplyr::ungroup() %>% 
      dplyr::select(-l3id) %>%
      colnames() 
    
    l2vars <- l23vars[!l23vars %in% l3vars] # must be at this position to be able to overwrite l3vars
    
    if(l3type=="FE") { # FE
      
      level3 <-
        data %>% 
        dplyr::rename(l3name=!!l3name) %>%
        dplyr::select(l3id, l3name) %>%
        group_by(l3id) %>%
        filter(row_number()==1) %>%
        ungroup() %>%
        mutate(rn = row_number(), val = 1) %>%
        pivot_wider(names_from = l3id, names_prefix="l3id", values_from = val, values_fill = list(val = 0)) %>%
        dplyr::rename(l3id=rn) 
      
      l3vars <- paste0("l3id", 2:dim(level3)[1]) # leave out first country
      
    } else { # RE
      level3 <-
        data %>%
        dplyr::select(l3id, all_of(l3vars)) %>%
        dplyr::group_by(l3id) %>%
        dplyr::summarise_all(.funs = first) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(l3id) %>% 
        dplyr::mutate_at(l3vars, function(x) { if(is.numeric(x) & dim(table(x))>2) x-mean(x) else x }) 
    }
    
  } else {
    l3vars <- c()
    level3 <- NULL
  }

  # Level 2
  level2 <- 
    data %>% 
    dplyr::arrange(l2id) %>%
    dplyr::group_by(l2id) %>%
    dplyr::add_count(name="l1n") %>%
    dplyr::filter(row_number()==1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(l2id) %>%
    dplyr::mutate(l1i2=cumsum(l1n), l1i1=lag(l1i2)+1) %>%
    dplyr::mutate(l1i1 = ifelse(row_number()==1, 1, l1i1)) %>%
    dplyr::mutate(X0 = 1) %>%
    dplyr::select(l2id, l1i1, l1i2, l1n, !!lhs, !!l2vars) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_at(l2vars, function(x) { if(is.numeric(x) & dim(table(x))>2) x-mean(x) else x }) 
  
  # ---------------------------------------------------------------------------------------------- #
  # 3. Transform data to JAGS format
  # ---------------------------------------------------------------------------------------------- #
  
  # IDs
  l1id <- level1 %>% .$l1id # length = rows @ level1, must be sorted by l2id
  l2id <- level2 %>% .$l2id # length = rows @ level2
  l3id <- if(length(hmi)>0) data %>% dplyr::group_by(l2id) %>% dplyr::filter(row_number()==1) %>% .$l3id else c() # length = rows @ level2
  
  l1n <- level2 %>% .$l1n # number of l1-members per l2-unit
  
  l1i1 <- level2 %>% .$l1i1 # first index of l1-members per l2-unit (@ level2)
  l1i1.l1 <- rep(l1i1, l1n) # first index of l1-members per l2-unit (@ level1)
  
  l1i2 <- level2 %>% .$l1i2 # last index of l1-members per l2-unit (@ level2)
  l1i2.l1 <- rep(l1i2, l1n) # last index of l1-members per l2-unit (@ level1)
  
  # Outcome
  if(family=="Gaussian") {
    Y <- level2 %>% dplyr::rename(Y = !!lhs) %>% .$Y
  } else if(family=="Weibull") {
    t <- level2 %>% dplyr::rename(t=lhs[1], ev=lhs[2]) %>% dplyr::mutate(t=case_when(ev==0 ~ NA_real_, TRUE ~ t)) %>% .$t
    event <- level2 %>% dplyr::rename(ev=lhs[2]) %>% .$ev
    t.cen <- level2 %>% dplyr::rename(t=lhs[1], ev=lhs[2]) %>% dplyr::mutate(t.cens = t + ev) %>% .$t.cens
    censored <- 1-event
  }
  
  # Xs
  X.l1 <- if(length(l1vars)>0) as.matrix(level1 %>% dplyr::select(!!l1vars)) else c()
  X.l2 <- if(length(l2vars)>0) as.matrix(level2 %>% dplyr::select(!!l2vars)) else c()
  X.l3 <- if(length(l3vars)>0) as.matrix(level3 %>% dplyr::select(!!l3vars)) else c()
  X.w  <- if(length(lwvars)>0) as.matrix(data %>% dplyr::select(!!lwvars)) else c()
  
  # Ns
  n.l1 <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2 <- length(l2id)
  n.l3 <- length(unique(l3id))
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  n.Xl3 <- dim(X.l3)[2]
  n.Xw <- dim(X.w)[2]
  
  n.GPN  <- level1 %>% group_by(l1id) %>% count() %>% .$n %>% as.numeric() %>% max() # max number of gov participations
  n.GPNi <- level1 %>% arrange(l1id) %>% group_by(l1id) %>% count() %>% .$n %>% as.numeric() # number of gov partipations per party, sorted l1id=1,2,3,...
  n.GPn  <- level1 %>% group_by(l1id) %>% dplyr::mutate(n=row_number()) %>% .$n %>% as.numeric() # participation index, sorted l1id=2,6,2,...
  
  # ---------------------------------------------------------------------------------------------- #
  # 3. Edit modelstring 
  # ---------------------------------------------------------------------------------------------- #

  # Model family
  if(family=="Gaussian") {
    modelstring <- if(length(hmi)>0) readr::read_file(paste0(DIR, "/JAGS/Gaussian_l123.txt")) else readr::read_file(paste0(DIR, "/JAGS/Gaussian_l12.txt")) # levels
  } else if(family=="Weibull") {
    modelstring <- if(length(hmi)>0) readr::read_file(paste0(DIR, "/JAGS/Weibull_l123.txt")) else readr::read_file(paste0(DIR, "/JAGS/Weibull_l12.txt")) # levels
  }
  
  # No covariates at level 3?
  if(length(hmi)>0 & is.null(n.Xl3)) { 
    modelstring <- stringr::str_remove(modelstring,  fixed("for(x in 1:n.Xl3) {\r\n    b.l3[x] ~ dnorm(0,0.0001)\r\n    ppp.b.l3[x] <- step(b.l3[x])\r\n  }\r\n  \r\n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l3[k] <- inprod(X.l3[k,], b.l3) + re.l3[k]"), "l3[k] <- re.l3[k]")
  }
  
  # No covariates at level 1?
  if(length(l1vars)==0) { 
    modelstring <- stringr::str_remove(modelstring,  fixed("  for(x in 1:n.Xl1) {\r\n    b.l1[x] ~ dnorm(0,0.0001)\r\n    ppp.b.l1[x] <- step(b.l1[x])\r\n  }\r\n  \r\n"))
    modelstring <- stringr::str_replace(modelstring, fixed("l1[i] <- inprod(X.l1[i,], b.l1) + re.l1[l1id[i]]"), "l1[i] <- re.l1[l1id[i]]")
  }
  
  # Auto-regressive l1 effect?
  if(mmwar==T) { 
    modelstring <- stringr::str_replace(modelstring, fixed("re.l1[l1id[i]]"), "re.l1[l1id[i], n.GPn[i]]")
    modelstring <- stringr::str_replace(modelstring, fixed("re.l1[i] ~ dnorm(0, tau.l1)\r\n  "), "re.l1[i,1] ~ dnorm(0, tau.l1)\r\n    for(j in 2:n.GPN) {\r\n      re.l1[i,j] ~ dnorm(re.l1[i,j-1], tau.l1)\r\n    }\r\n  ")
  }
  
  # Weight function
  if(mmwfunction != "uw[i] <- 1/X.w[i,1]") modelstring <- stringr::str_replace(modelstring, fixed("uw[i] <- 1/X.w[i,1]"), mmwfunction)
  if(mmwconstraint == 2) modelstring <- stringr::str_replace(modelstring, fixed("w[i] <- uw[i] / sum(uw[l1i1.l1[i]:l1i2.l1[i]])"), "w[i] <- uw[i] * n.l2/sum(uw[]) # rescale to sum up to n.l1 overall")
  if(stringr::str_detect(mmwfunction, "b.w\\[.\\]")) { 
    modelstring <- stringr::str_replace(modelstring, fixed("b.w\r\n"), paste0(mmwcoefstring, collapse = "")) 
  } else {
    modelstring <- stringr::str_replace(modelstring, fixed("b.w\r\n"), "") 
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
      } else if(vars[i] == "tau.l2") { # Gaussian
        modelstring <- stringr::str_replace(modelstring, fixed("tau.l2 ~ dscaled.gamma(25, 1)"), paste0("tau.l2 ~ ", prior[i]))
      } else if(vars[i] == "tau.l3") {
        modelstring <- stringr::str_replace(modelstring, fixed("tau.l3 ~ dscaled.gamma(25, 1)"), paste0("tau.l3 ~ ", prior[i]))
      } else if(vars[i] == "sigma.l2") { # Weibull
        modelstring <- stringr::str_replace(modelstring, fixed("sigma.l2 ~ dexp(0.0001)"), paste0("sigma.l2 ~ ", prior[i]))
      }
    }
  }

  if(is.null(modelfile)) readr::write_file(modelstring, paste0(DIR, "/temp/modelstring.txt"))
  
  # ---------------------------------------------------------------------------------------------- #
  # 4. Set up JAGS
  # ---------------------------------------------------------------------------------------------- #
  
  # Level 1
  if(mmwar==T) { # AR==T
    re.param <- c()
    for(j in 1:n.ul1) {
      re.param <- append(re.param, paste0("re.l1[", j, "," , seq(1,n.GPNi[j]), "]"))
    }
    l1.param <- if(monitor==T) c("sigma.l1", re.param) else c("sigma.l1")
    l1.data  <- c("l1id", "l1i1", "l1i2", "n.l1", "n.ul1", "n.GPn", "n.GPN")
  } else { # AR==F
    l1.param <- if(monitor==T) c("sigma.l1", "re.l1") else c("sigma.l1")
    l1.data  <- c("l1id", "l1i1", "l1i2", "n.l1", "n.ul1")
  }
  if(!is.null(n.Xl1)) { l1.data <- append(l1.data, c("X.l1", "n.Xl1")); l1.param <- append(l1.param, c("b.l1", "ppp.b.l1")) }

  # Level 2
  l2.param <- c("sigma.l2")
  l2.data  <- c("n.l2")
  if(!is.null(n.Xl2)) { l2.data <- append(l2.data, c("X.l2", "n.Xl2")); l2.param <- append(l2.param, c("b.l2", "ppp.b.l2")) }
  
  # Level 3
  l3.param <- if(length(hmi)>0 & monitor==T) c("sigma.l3", "re.l3") else if(length(hmi)>0 & monitor==F) c("sigma.l3") else c()
  l3.data  <- if(length(hmi)>0) c("l3id", "n.l3") else c()
  if(!is.null(n.Xl3)) l3.data <- append(l3.data, c("X.l3", "n.Xl3"))
  if(!is.null(n.Xl3) & (l3type=="RE" | (l3type=="FE" & showFE==T))) l3.param <- append(l3.param, c("b.l3", "ppp.b.l3")) 
  
  # Weight function
  lw.param <- if(length(lwvars)>length(offsetvars) & monitor == T) c("b.w", "ppp.b.w", "w") else if(length(lwvars)>length(offsetvars) & monitor == F) c("b.w", "ppp.b.w") else if(length(lwvars)>0 & monitor==T) c("w") else c()
  lw.data  <- c("X.w")
  if(mmwconstraint==1) lw.data <- append(lw.data, c("l1i1.l1", "l1i2.l1"))
  
  jags.params <- c(l1.param, l2.param, l3.param, lw.param)
  jags.data   <- c(l1.data, l2.data, l3.data, lw.data)
  
  # Initial values
  if(is.null(seed)) seed <- runif(1, 0, 1000)

  # Model-specific specifications ---------------------------------------------------------------- #
  
  if(family=="Gaussian") {

    # Dependent variable
    jags.data <- append(jags.data, "Y")
    
    jags.inits <- list( 
      list(".RNG.seed" = seed+1, tau.l1=1.5, tau.l2=0.5),
      list(".RNG.seed" = seed+2, tau.l1=1.0, tau.l2=1.0),
      list(".RNG.seed" = seed+3, tau.l1=0.5, tau.l2=1.5))
    
    if(length(hmi)>0) jags.inits <- lapply(jags.inits, FUN=function(x) { append(x, list(tau.l3=1.0)) }) 
    # 2do: change so that adding inits is easier and maybe also a parameter of rmm()
    
  } else if(family=="Weibull") {
    
    # Dependent variable
    jags.data   <- append(jags.data, c("t", "t.cen", "censored"))
    
    # Initial values
    t.init <- t
    t.init[] <- NA
    t.init[censored==1] <- t.cen[censored==1] + 1 
    
    jags.inits <- list( 
      list(".RNG.seed" = seed+1, t=t.init, shape=1.3, tau.l1=0.5),
      list(".RNG.seed" = seed+2, t=t.init, shape=1.1, tau.l1=1.0),
      list(".RNG.seed" = seed+3, t=t.init, shape=0.9, tau.l1=1.5))
    
  }
  
  if(is.null(modelfile)) modelfile <- paste0(DIR, "/temp/modelstring.txt")
  
  # ---------------------------------------------------------------------------------------------- #
  # Run JAGS 
  # ---------------------------------------------------------------------------------------------- #
  
  if(run==T) {
    
    jags.out <<- jags(data=jags.data, inits = jags.inits, parameters.to.save = jags.params, n.chains = 3, n.iter = iter, n.burnin = burnin, model.file = modelfile)
    
    # -------------------------------------------------------------------------------------------- #
    # Format results
    # -------------------------------------------------------------------------------------------- #
    
    reg.table <- 
      as.data.frame(jags.out$BUGSoutput$summary[, c(1, 2, 3, 7)]) %>% 
      tibble::rownames_to_column(var="name") %>% 
      dplyr::rename(estimate=2, sd=3, lb=4, ub=5)
    
    pppvalues <- # extract ppp here because they should remain a mean estimate 
      reg.table %>%  
      dplyr::filter(startsWith(name, "ppp.")) %>% 
      dplyr::mutate_at(.vars=c("name"), .funs=function(x) str_remove(x, "ppp.")) %>%
      dplyr::select(-sd, -lb, -ub)
    
    # 2DO: ppp values for variance terms
    
    # Mean and 95% CI or Mode and HDI ------------------------------------------------------------ #
    
    if(hdi != FALSE) {
      
      mcmc.out <- MCMCvis::MCMCchains(mcmcplots::as.mcmc.rjags(jags.out)) 
      
      reg.table.hdi <- data.frame(name=NA, estimate=NA, sd=NA, lb=NA, ub=NA) 
      
      for(i in 1:dim(mcmc.out)[[2]]) {
        dres <- density(mcmc.out[,i]) # density
        mode <- dres$x[which.max(dres$y)] # get mode from density
        hdiinterval <- HDInterval::hdi(mcmc.out[,i], credMass=hdi) # get HDI
        sd <- sd(mcmc.out[,i])
        reg.table.hdi[i,] <- c(colnames(mcmc.out)[i], mode, sd, hdiinterval[1] , hdiinterval[2])
      }
      
      reg.table <- reg.table.hdi # replace
    } 
    
    # Remove random effects, weights, and ppp from reg.table and save separately ----------------- #
    
    if(monitor==T) {
      # Level-1 RE
      
      re.l1 <- reg.table %>% dplyr::filter(startsWith(name, "re.l1")) %>% dplyr::select(-sd, -lb, -ub) %>% dplyr::mutate(estimate = round(as.numeric(estimate), r))
      
      if(mmwar==T) { # AR == T
        
        re.l1 <- 
          re.l1 %>%
          tidyr::separate(name, c("i", "j"), ",", remove = F) %>%
          dplyr::mutate(i=as.numeric(str_remove(i, "re.l1\\[")), j=as.numeric(str_remove(j, "]"))) %>%
          dplyr::arrange(i,j) # sort by l1 unit and random walk (important)
        
        remat <- matrix(NA, nrow=n.ul1, ncol=n.GPN) # rows: random walks / columns: parties
        rownames(remat) <- paste0("L1 unit ", seq(1,n.ul1))
        colnames(remat) <- paste0("Random walk ", seq(1,n.GPN))
        
        for(i in 1:dim(re.l1)[1]) { # populate matrix
          remat[re.l1[i,"i"], re.l1[i,"j"]] <- re.l1[i, "estimate"] 
        }
        
        re.l1 <<- remat
        
      } else { # AR == F
        re.l1 <<- re.l1 %>% .$estimate
      }
      
      reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "re.l1["))
      
      # Level-3 RE
      re.l3 <<- reg.table %>% dplyr::filter(startsWith(name, "re.l3")) %>% dplyr::mutate(estimate = round(as.numeric(estimate), r)) %>% .$estimate
      reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "re.l3[")) 
      
      # Weights
      wmat <- matrix(NA, nrow = n.l2, ncol = max(l1n))
      w <- reg.table %>% dplyr::filter(startsWith(name, "w")) %>% dplyr::mutate(estimate = round(as.numeric(estimate), r)) %>% .$estimate
      
      id1 <- cumsum(l1n)-l1n+1
      id2 <- cumsum(l1n)
      
      for(i in 1:n.l2) {
        wmat[i,1:l1n[i]] <- w[id1[i]:id2[i]]
      }
      
      w <<- wmat
      
      reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "w")) 
    }
    
    # PPP values (add as separate column, must be after HDI) 
    reg.table <- 
      reg.table %>% 
      dplyr::left_join(pppvalues %>% dplyr::select(name, estimate) %>% dplyr::rename(ppp=estimate), by=c("name")) %>%
      dplyr::mutate(ppp=case_when(estimate>0 ~ 1-ppp, # for positive estimates, the correct p-value = 1-p
                                  TRUE ~ ppp)) %>%
      dplyr::filter(!startsWith(name, "ppp.")) 
    
    # Rename ------------------------------------------------------------------------------------- #
    
    newnames <- reg.table %>% .$name
    newnames[stringr::str_detect(newnames, "b.l1")] <- l1vars
    newnames[stringr::str_detect(newnames, "b.l2")] <- l2vars
    newnames[stringr::str_detect(newnames, "b.l3")] <- if(!is.null(l3name) & l3type=="FE") {level3 %>% .$l3name}[-1] else l3vars
    newnames[stringr::str_detect(newnames, "b.w")]  <- lwvars[!lwvars %in% offsetvars]
    
    reg.table <- 
      reg.table %>% 
      rename(coefficients=estimate) %>% # compatibility with lm() etc
      dplyr::mutate(variable=newnames) %>% relocate(variable, .before = coefficients) %>%
      filter(!variable=="deviance") %>%
      rbind(c("DIC", "DIC", jags.out$BUGSoutput$DIC, NA, NA, NA, NA)) %>%
      dplyr::mutate_at(3:7, list(~round(as.numeric(.), r))) %>%
      tibble::column_to_rownames(var = "name") 
    
  }
  
  # -------------------------------------------------------------------------------------------- #
  # Return
  # -------------------------------------------------------------------------------------------- #
  
  return(reg.table)

}
