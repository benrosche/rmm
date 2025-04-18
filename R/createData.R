# ================================================================================================ #
# Function createData
# ================================================================================================ #

createData <- function(data, ids, vars, l1, l3) {
  
  # Unpack lists --------------------------------------------------------------------------------- #
  
  l12ids <- ids[["l12ids"]]
  l3id <- ids[["l3id"]]
  
  lhs <- vars[["lhs"]]
  l1vars <- vars[["l1vars"]]
  l23vars <- vars[["l23vars"]]
  wvars <- vars[["wvars"]]
  wvars_p <- vars[["wvars_p"]]
  wparams <- vars[["wparams"]]
  
  mm <- l1[["mm"]]
  
  hm <- l3[["hm"]]
  l3type <- l3[["l3type"]]
  l3name <- l3[["l3name"]]
  
  # Rename and regroup ids and sort -------------------------------------------------------------- #
  
  if(isFALSE(mm)) data <- data %>% dplyr::mutate(l1id = 1, l2id = 1)
  if(isFALSE(hm)) data <- data %>% dplyr::mutate(l3id = 1)
  
  data <-
    data %>% 
    dplyr::rename(l1id = all_of(l12ids[1]), l2id = all_of(l12ids[2]), l3id = all_of(l3id)) %>%
    dplyr::group_by(l1id) %>% dplyr::mutate(l1id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(l2id) %>% dplyr::mutate(l2id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(l3id) %>% dplyr::mutate(l3id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::arrange(l2id, l1id) 
  
  # Level 1 -------------------------------------------------------------------------------------- #
  
  if(mm) {
    
    l1dat <-
      data %>% 
      dplyr::arrange(l2id, l1id) %>% # important
      dplyr::select(l1id, l2id, all_of(l1vars))
    
    wdat <- 
      data %>%
      dplyr::add_count(l2id, name="n") %>% 
      dplyr::select(l1id, l2id, all_of(wvars)) 
    
  } else { # no l1
    
    l1dat <- NULL
    wdat  <- NULL
    
  }
  
  # Level 3 -------------------------------------------------------------------------------------- #
  
  if(hm) {
    
    l3vars <- 
      data %>% 
      dplyr::select(l3id, all_of(l23vars[-1])) %>%
      dplyr::group_by(l3id) %>%
      dplyr::mutate(across(all_of(l23vars[-1]), ~var(., na.rm = TRUE))) %>% # select variables that do not vary within levels
      dplyr::ungroup() %>% 
      dplyr::select_if(~sum(.)==0) %>%  
      names() 
    
    l2vars <- l23vars[!l23vars %in% l3vars] # must be at this position to be able to overwrite l3vars
    
    if(l3type=="FE") { # FE
      
      l3dat <-
        data %>% 
        dplyr::rename(l3name=all_of(l3name)) %>%
        dplyr::select(l3id, any_of("l3name")) %>%
        dplyr::group_by(l3id) %>%
        dplyr::filter(row_number()==1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(rn = row_number(), val = 1) %>%
        tidyr::pivot_wider(names_from = l3id, names_prefix="l3id", values_from = val, values_fill = list(val = 0)) %>%
        dplyr::rename(l3id=rn) 
      
      l3vars <- paste0("l3id", 2:dim(l3dat)[1]) # leave out first country
      
    } else { # RE
      
      l3dat <-
        data %>%
        dplyr::select(l3id, all_of(l3vars)) %>%
        dplyr::group_by(l3id) %>%
        dplyr::filter(row_number()==1) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(l3id) 
      
    }
    
  } else { # no l3
    
    l2vars <- l23vars
    l3vars <- c()
    l3dat <- NULL
    
  }
  
  # Level 2 -------------------------------------------------------------------------------------- #
  
  if(isFALSE(mm)) data <- data %>% dplyr::mutate(l2id = row_number())
  
  l2dat <- 
    data %>% 
    dplyr::arrange(l2id) %>%
    dplyr::group_by(l2id) %>%
    dplyr::add_count(l2id, name="l1n") %>%
    dplyr::filter(row_number()==1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(l2id) %>%
    dplyr::mutate(l1i2=cumsum(l1n), l1i1=lag(l1i2)+1) %>%
    dplyr::mutate(l1i1 = ifelse(row_number()==1, 1, l1i1)) %>%
    dplyr::mutate(X0 = 1) %>%
    dplyr::select(l2id, l3id, l1i1, l1i2, l1n, all_of(lhs), all_of(l2vars)) 
  
  # Collect return ------------------------------------------------------------------------------- #
  
  return(
    list(
      "data"=data, 
      "level1"=list("dat"=l1dat, "vars"=l1vars), 
      "level2"=list("dat"=l2dat, "vars"=l2vars, "lhs"=lhs), 
      "level3"=list("dat"=l3dat, "vars"=l3vars), 
      "weight"=list("dat"=wdat,  "vars"=wvars, "vars_p"=wvars_p, "params"=wparams)
    )
  )
  
}
