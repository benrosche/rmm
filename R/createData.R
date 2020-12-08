# ================================================================================================ #
# Function createData
# ================================================================================================ #

createData <- function(data, ids, vars, l3, transform) {
  
  # Unpack lists --------------------------------------------------------------------------------- #
  
  l12ids <- ids$l12ids
  l3id   <- ids$l3id
  
  lhs <- vars$lhs
  l1vars <- vars$l1vars
  l23vars <- vars$l23vars
  lwvars <- vars$lwvars
  offsetvars <- vars$offsetvars
  
  hm <- l3$hm
  l3type <- l3$l3type
  
  # Center or Standardize ------------------------------------------------------------------------ #
  
  cen_std <- function(x) { 
    if(transform=="center") {
      if(is.numeric(x) & dim(table(x))>2) x-mean(x) else x 
    } else if(transform=="std") {
      if(is.numeric(x) & dim(table(x))>2) (x-mean(x))/sqrt(var(x)) else x 
    } else {
      x
    }
  }
  
  # Rename and regroup ids and sort -------------------------------------------------------------- #
  
  data <-
    data %>% 
    dplyr::rename(l1id = !!l12ids[1], l2id = !!l12ids[2]) %>%
    dplyr::mutate(one_vec = 1) %>% # must be here
    dplyr::mutate(l3id = !!rlang::sym(l3id)) %>% 
    dplyr::group_by(l1id) %>% dplyr::mutate(l1id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(l2id) %>% dplyr::mutate(l2id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::group_by(l3id) %>% dplyr::mutate(l3id = cur_group_id()) %>% dplyr::ungroup() %>%
    dplyr::arrange(l2id, l1id) 
  
  # Level 1 -------------------------------------------------------------------------------------- #
  
  level1 <-
    data %>% 
    dplyr::arrange(l2id, l1id) %>% # important
    dplyr::select(l1id, l2id, !!l1vars) %>%
    dplyr::mutate_at(l1vars, cen_std) # center continuous vars 
  
  # Level 3 -------------------------------------------------------------------------------------- #
  
  if(hm) {
    
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
        dplyr::mutate_at(l3vars, cen_std) 
    }
    
  } else {
    l2vars <- l23vars
    l3vars <- c()
    level3 <- NULL
  }
  
  # Level 2 -------------------------------------------------------------------------------------- #
  
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
    dplyr::mutate_at(l2vars, cen_std) 
  
  # Collect return ------------------------------------------------------------------------------- #
  
  return(list("data"=data, "level1"=list("dat"=level1, "vars"=l1vars), "level2"=list("dat"=level2, "vars"=l2vars), "level3"=list("dat"=level3, "vars"=l3vars)))
  
}
