# ================================================================================================ #
# Function dissectFormula 
# ================================================================================================ #

dissectFormula <- function(formula, family, data) {
  
  # This function takes the formula object and turns it into variables for further processing.
  # Arguments:
  # - data (df) data.frame to check whether specified variables are in the dataset
  # - family (str) model family
  # - formula (formula) formula object
  # Return:
  # - Returns a list with four elements: ids, vars, l1, l3.
  
  # Extract left- and right-hand side from formula ----------------------------------------------- #
  
  lhs <- formula[[2]] %>% all.vars()
  rhs <- terms(formula) %>% attr(., "term.labels")
  
  if(family=="Gaussian" & length(lhs)>1) stop("family=\"Gaussian\" takes only one variable on the left-hand side of the formula.")
  if(family=="Weibull" & length(lhs)!=2) stop("family=\"Weibull\" takes two variables on the left-hand side: 'Surv(survtime, event)'")
  
  # ---------------------------------------------------------------------------------------------- #
  # mm() object 
  # ---------------------------------------------------------------------------------------------- #
  
  mmstring <- {rhs[str_starts(rhs, fixed("mm("))] %>% stringr::str_match(., "mm\\((.*)\\)")}[2]
  
  mm <- !is.na(mmstring)
  
  if(mm) {
    
    # id() --------------------------------------------------------------------------------------- #

    # Extract l1id and l2id
    l12ids <- stringr::str_extract(mmstring, "id\\([^)]+\\)") 
    if (!is.na(l12ids)) {
      l12ids <- {stringr::str_match(l12ids, "id\\(([^)]+)\\)")[[2]]} %>% {stringr::str_split(., ", ")[[1]]} %>% trimws()
      if(isFALSE(all(l12ids %in% names(data)))) stop("Either l1id or l2id could not be found in the dataset.")
    } else {
      stop("'id(l1id,l2id)' missing within 'mm()'")
    }

    # mmc() -------------------------------------------------------------------------------------- #
    
    # Extract l1vars if any
    l1vars <- stringr::str_extract(mmstring, "mmc\\([^)]*\\)")
    if (!is.na(l1vars)) {
      l1vars <- {stringr::str_match(l1vars, "mmc\\(([^)]+)\\)")[[2]]} %>% {stringr::str_split(., "\\+")[[1]]} %>% trimws()
      if(length(l1vars) == 1 && is.na(l1vars)) l1vars <- c() # no l1vars specified
      if(isFALSE(all(l1vars %in% names(data)))) stop("l1vars within 'mmc()' could not be found in the dataset.")
    } else {
      stop("l1vars 'mmc()' missing within 'mm()'")
    }

    # mmw() -------------------------------------------------------------------------------------- #
    
    mmwstring <- stringr::str_match(mmstring, "mmw\\((.*)\\)")[2]
    
    mmw <- !is.na(mmwstring)
    
    if(mmw) { # is mmw() specified?
      
      # Extract weight function
      mmwfunction <- str_extract(mmwstring, "^[^,]+") %>% trimws() # extract w ~ ...
      if(is.na(mmwfunction)) stop("Weight function: w ~ ... not specified.")
      
      # Extract weight variables
      lwvars <- mmwfunction %>% as.formula() %>% all.vars()
      lwvars <- lwvars[!str_detect(lwvars, "^b[0-9]+$|^w$")] %>% trimws() # extract rhs of weight function
      if(isFALSE(all(lwvars %in% names(data)))) stop("Weight variables in weight function w ~ ... could not be found in the dataset.")
      if(length(lwvars)==0) stop("No weight variables in weight function w ~ ... specified.")
      
      lwparams <- str_extract_all(mmwfunction, "\\bb\\d+\\b")[[1]] %>% unique()
      
      # Determine constraint
      mmwconstraint <- stringr::str_match(mmwstring, "(constraint|c)\\s*=\\s*([^,)]+)") %>% trimws()
      if(!is.na(mmwconstraint[1])) {
        mmwconstraint <- mmwconstraint[3] %>% as.numeric()
        if(!mmwconstraint %in% c(1,2)) stop("Constraint must be either 1 or 2")
      } else {
        mmwconstraint <- 1 # default
      }
      
      # Determine whether auto-regressive random effect is TRUE or FALSE
      mmwar <- stringr::str_match(mmstring, "ar\\s*=\\s*(T|F)") %>% trimws()
      if (!is.na(mmwar[1])) {
        mmwar <- mmwar[2]
      } else {
        mmwar <- F
      }
      
    } else {
      stop("Weight function 'mmw()' missing within 'mm()'")
    }
    
  } else {
    
    l12ids <- c("l1id", "l2id") # must be here
    l1vars <- c()
    lwvars <- c()
    lwparams <- c()
    mmwfunction <- c()
    mmwconstraint <- F
    mmwar <- F
    
  }
  
  # ---------------------------------------------------------------------------------------------- #
  # hm() object 
  # ---------------------------------------------------------------------------------------------- #
  
  hmstring <- {rhs[str_starts(rhs, fixed("hm("))] %>% stringr::str_match(., "hm\\((.*)\\)")}[2]
  
  hm <- !is.na(hmstring)
  
  if(hm) {
    
    # Extract l3id
    l3id <- str_match(hmstring, "id\\s*=\\s*([^,)]+)")[2] %>% trimws()
    if(isFALSE(l3id %in% names(data))) stop("'l3id' in 'hm()' could not be found in the dataset.")
    
    # Extract labels for l3 effects
    l3name  <- str_match(hmstring, "name\\s*=\\s*([^,)]+)")[2] %>% trimws()
    if(!is.na(l3name)) {
      if(isFALSE(l3name %in% names(data))) stop("'l3name' in 'hm()' is specified but could not be found in the dataset.")
    } else {
      l3name <- c() 
    }
    
    # Determine whether fixed or random effects should be estimated
    l3type <- str_match(hmstring, "type\\s*=\\s*([^,)]+)")[2] %>% trimws()
    if(!is.na(l3type)) {
      if(!l3type %in% c("RE","FE")) stop("'type' in 'hm()' must be either 'RE', 'FE', or left unspecified")
    } else {
      l3type <- c("RE")
    }
  
    # Determine whether fixed effects should be monitored
    showFE  <- str_match(hmstring, "showFE\\s*=\\s*([^,)]+)")[2] %>% trimws() 
    if(!is.na(showFE) & !showFE %in% c("T", "F", "TRUE", "FALSE")) stop("'show' in 'hm()' must be either 'TRUE', 'FALSE', or left unspecified")
    showFE <- if(is.na(showFE)) FALSE else as.logical(showFE)

  } else {
    
    l3id <- "l3id" 
    l3name <- c()
    l3type <- FALSE
    showFE <- FALSE
    
  }
  
  # ---------------------------------------------------------------------------------------------- #
  # Collect terms and return
  # ---------------------------------------------------------------------------------------------- #
  
  l23vars <- rhs[!sapply(rhs, \(x) any(str_starts(x, fixed(c("hm(", "mm(", "1", "X0")))))] %>% append(., "X0", after=0)
  
  return(
    list(
      "ids"=
        list(
          "l12ids"=l12ids, 
          "l3id"=l3id
        ),
      "vars"=
        list(
          "lhs"=lhs, 
          "l1vars"=l1vars, 
          "l23vars"=l23vars, 
          "lwvars"=lwvars,
          "lwparams"=lwparams
        ),
      "l1"=
        list(
          "mm"=mm, 
          "mmwfunction"=mmwfunction, 
          "mmwconstraint"=mmwconstraint, 
          "mmwar"=mmwar
        ),
      "l3"=
        list(
          "hm"=hm, 
          "l3name"=l3name, 
          "l3type"=l3type, 
          "showFE"=showFE
        )
    )
  )
}
