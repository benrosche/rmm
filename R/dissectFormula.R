# ================================================================================================ #
# Function dissectFormula 
# ================================================================================================ #

dissectFormula <- function(formula, family) {
  
  # This function has two arguments:
  # ...
  
  getTerms <- function(e, x = list()) {
    if (is.symbol(e)) c(e, x)
    else if (identical(e[[1]], as.name("+"))) Recall(e[[2]], c(e[[3]], x))
    else c(e, x)
  }
  
  # Left-hand side ------------------------------------------------------------------------------- #
  
  lhs <- all.vars(formula[[2]])
  if(family=="Gaussian" & length(lhs)>1) stop("family=\"Gaussian\" takes only one variable on the left-hand side of the formula.")
  if(family=="Weibull" & length(lhs)!=2) stop("family=\"Weibull\" takes two variables on the left-hand side: 'Surv(survtime, event)'")
  
  # Right-hand side ------------------------------------------------------------------------------ #
  
  rhs <- gsub("[[:space:]]", "", getTerms(formula[3][[1]]))
  
  # ---------------------------------------------------------------------------------------------- #
  # mm() object 
  # ---------------------------------------------------------------------------------------------- #
  
  # mm() ----------------------------------------------------------------------------------------- #
  
  mmi <- which(startsWith(rhs, "mm("))
  mmstring <- rhs[mmi]
  
  if(length(mmstring)==0) stop("Multiple membership construct 'mm()'  missing. Please use the lme4 package to estimate conventional multilevel models.")
  
  # id() ----------------------------------------------------------------------------------------- #
  
  l12ids <- if(stringr::str_detect(mmstring, "id\\(.*\\)")) unlist(stringr::str_split(stringr::str_replace(mmstring, ".*id\\((.*?)\\).*", "\\1"), ",")) else stop("'id(l1id,l2id)' missing within 'mm()'")
  
  if(all(l12ids %in% names(data))==FALSE) stop("Either l1id or l2id could not be found in the dataset.")
  
  # mmc() ---------------------------------------------------------------------------------------- #
  l1vars <- if(stringr::str_detect(mmstring, "mmc\\(.*\\)")) unlist(stringr::str_split(stringr::str_replace(mmstring, ".*mmc\\((.*?)\\).*", "\\1"), "\\+")) else stop("l1vars 'mmc()' missing within 'mm()'")
  if(all(nzchar(l1vars)) == FALSE) l1vars <- c()
  if(all(l1vars %in% names(data))==FALSE) stop("l1vars within 'mmc()' could not be found in the dataset.")
  
  # mmlw(), function  ---------------------------------------------------------------------------- #
  
  mmwstring <- if(stringr::str_detect(mmstring, "mmw\\(.*\\)")) unlist(stringr::str_split(stringr::str_replace(mmstring, ".*mmw\\((.*)\\)(.*)\\)", "\\1 \\2"), ",")) else stop("Weight function 'mmw()' missing within 'mm()'")
  
  mmwfunction <- mmwstring[1] # weight function
  
  if(!startsWith(mmwfunction, "w")) stop("Weight function within 'mmw()' missing")
  
  lwvars <- tryCatch(all.vars(as.formula(mmwfunction))[-1], error=function(e) { print("Weight variables within 'mmw()' missing") }) # all variables but 'w'
  
  if(all(lwvars %in% names(data))==FALSE | length(lwvars)==0) stop("Weight variables within 'mmw()' could not be found in the dataset.")
  if(length(lwvars)<length(all.vars(as.formula(mmwfunction),unique=F)[-1])) stop("Weight function must not include the same variable multiple times.") # must be at this position
  
  offsetvars <- stringr::str_replace(stringr::str_extract_all(mmwfunction, "offset\\((.+?)\\)", simplify = T), "offset\\((.+?)\\)", "\\1")
  
  # Translate weight function into JAGS format  -------------------------------------------------- #
  
  mmwcoefstring <- c()
  for(i in 1:length(lwvars)) { 
    if(lwvars[i] %in% offsetvars) {
      mmwfunction <- stringr::str_replace(mmwfunction, fixed(paste0("offset(", lwvars[i], ")")), paste0("(X.w[i,", i, "])")) # without coefficient
    } else {
      mmwfunction <- stringr::str_replace(mmwfunction, fixed(lwvars[i]), paste0("(b.w[", i, "]*X.w[i,", i, "])")) # with coefficient
      mmwcoefstring  <- append(mmwcoefstring, paste0("b.w[", i, "] ~ dnorm(0,0.0001)\n  ppp.b.w[", i, "] <- step(b.w[", i, "])\n  "))
    }
  }
  mmwfunction <- stringr::str_replace(mmwfunction, "w~", "uw[i] <- ")
  
  # mmlw(), constraint --------------------------------------------------------------------------- #
  
  if(length(mmwstring)>1) {
    mmwconstraint <- as.numeric(str_extract(paste0(mmwstring[-1], collapse = ", "), "(?!constraint\\=)[0-9]"))
    if(is.na(mmwconstraint) & length(lwvars)>length(offsetvars)) mmwconstraint <- 2
    if(is.na(mmwconstraint) & length(lwvars)<=length(offsetvars)) mmwconstraint <- 1
  } else {
    if(length(lwvars)>length(offsetvars)) mmwconstraint <- 2
    if(length(lwvars)<=length(offsetvars)) mmwconstraint <- 1
  }
  
  if(!mmwconstraint %in% c(1,2)) stop("Constraint must be either 1 or 2")
  
  # mmlw(), ar ----------------------------------------------------------------------------------- #
  
  if(length(mmwstring)>1) {
    mmwar <- as.logical(str_extract(paste0(mmwstring[-1], collapse = ", "), "(?!ar\\=)(T|F|TRUE|FALSE)"))
    if(is.na(mmwar)) mmwar <- FALSE
  } else {
    mmwar <- FALSE
  }
  if(!mmwar %in% c(T,F)) stop("AR must be either True or False")
  
  # ---------------------------------------------------------------------------------------------- #
  # hm() object 
  # ---------------------------------------------------------------------------------------------- #
  
  hmi <- which(startsWith(rhs, "hm("))
  hmstring <- rhs[hmi]
  hm <- if(length(hmi)>0) TRUE else FALSE
  
  if(hm) {
    
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
    l3name <- c()
    l3type <- FALSE
    showFE <- FALSE
    
  }
  
  # Level 2 + 3 vars ----------------------------------------------------------------------------- #
  
  l23vars <- stringr::str_replace(rhs[c(-mmi, -hmi)], "^1$", "X0") # Intercept
  if(!"X0" %in% l23vars)  l23vars <- append(l23vars, "X0", after=0)
  
  # Collect terms -------------------------------------------------------------------------------- #
  
  return(list("ids"=list("l12ids"=l12ids, "l3id"=l3id),
              "vars"=list("lhs"=lhs, "l1vars"=l1vars, "l23vars"=l23vars, "lwvars"=lwvars, "offsetvars"=offsetvars),
              "l3"=list("hm"=hm, "l3name"=l3name, "l3type"=l3type, "showFE"=showFE),
              "mm"=list("mmwfunction"=mmwfunction, "mmwcoefstring"=mmwcoefstring, "mmwconstraint"=mmwconstraint, "mmwar"=mmwar)))
}
