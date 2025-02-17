# ================================================================================================ #
# Function formatJags
# ================================================================================================ #

formatJags <- function(jags.out, monitor, Ns, l1, l3, level1, level2, level3, weightf) {

  # Unpack lists --------------------------------------------------------------------------------- #

  l1vars <- level1$vars
  l2vars <- level2$vars
  l3vars <- level3$vars
  
  l3dat  <- level3$dat
  
  hm <- l3$hm
  l3name <- l3$l3name
  l3type <- l3$l3type
  
  mm <- l1$mm
  mmwar <- l1$mmwar
  lwvars <- weightf$vars
  
  n.ul1 <- Ns$n.ul1
  l1n  <- Ns$l1n
  n.l2 <- Ns$n.l2
  n.GPN <- Ns$n.GPN
  
  # Create reg.table from JAGS output ------------------------------------------------------------ #
  
  reg.table <- 
    as_tibble(jags.out$BUGSoutput$summary[, c(1, 2, 3, 7)], rownames="name") %>% 
    dplyr::rename(coefficients=2, sd=3, lb=4, ub=5) # compatibility with lm() 
  
  # Remove random effects and weights from reg.table and save separately ------------------------- #
  
  if(monitor) {
    
    # Level-1 RE #
    
    if(mm) {
      
      re.l1 <- reg.table %>% dplyr::filter(startsWith(name, "re.l1")) %>% dplyr::select(-sd, -lb, -ub)
      
      if(mm & mmwar) { # AR
        
        re.l1 <- 
          re.l1 %>%
          tidyr::separate(name, c("i", "j"), ",", remove = F) %>%
          dplyr::mutate(i=as.numeric(str_remove(i, "re.l1\\[")), j=as.numeric(str_remove(j, "]"))) %>%
          dplyr::arrange(i,j) # sort by l1 unit and random walk (important)
        
        remat <- matrix(NA, nrow=n.ul1, ncol=n.GPN) # rows: random walks / columns: parties
        rownames(remat) <- paste0("L1 unit ", seq(1,n.ul1))
        colnames(remat) <- paste0("Random walk ", seq(1,n.GPN))
        
        for(i in 1:dim(re.l1)[1]) { # populate matrix
          remat[re.l1[i,"i"], re.l1[i,"j"]] <- re.l1[i, "coefficients"] 
        }
        
        re.l1 <- remat
        
      } else { # AR == F
        re.l1 <- re.l1 %>% pull(coefficients)
      }
      
    } else {
      
      re.l1 <- c()
      
    }

    reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "re.l1["))
    
    # Level-3 RE #
    
    re.l3 <- if(hm) reg.table %>% dplyr::filter(startsWith(name, "re.l3")) %>% pull(coefficients) else c()
    reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "re.l3[")) 
    
    # Weights #
    
    if(mm) {
      
      w <- reg.table %>% dplyr::filter(startsWith(name, "w")) %>% pull(coefficients)
      
      id1 <- cumsum(l1n)-l1n+1
      id2 <- cumsum(l1n)
      
      wmat <- matrix(NA, nrow = n.l2, ncol = max(l1n))
      rownames(wmat) <- paste0("L2 unit ", seq(1,n.l2))
      colnames(wmat) <- paste0("W", seq(1,max(l1n)))
      for(i in 1:n.l2) {
        wmat[i,1:l1n[i]] <- w[id1[i]:id2[i]]
      }
      
      w <- wmat
      reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "w")) 
      
    } else {
      w <- c()
    }
    
  }
  
  # Predicted values #
  pred <- reg.table %>% dplyr::filter(startsWith(name, "pred")) %>% dplyr::select(-sd, -lb, -ub) %>% pull(coefficients)
  reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "pred"))
  
  # Rename --------------------------------------------------------------------------------------- #
  
  newnames <- reg.table %>% pull(name)
  newnames[stringr::str_detect(newnames, "b.l1")] <- l1vars
  newnames[stringr::str_detect(newnames, "b.l2")] <- l2vars
  newnames[stringr::str_detect(newnames, "b.l3")] <- if(!is.null(l3name) & l3type=="FE") {l3dat %>% pull(l3name)}[-1] else l3vars
  
  reg.table <- 
    reg.table %>% 
    dplyr::mutate(variable=newnames) %>% relocate(variable, .before = coefficients) %>%
    filter(!variable=="deviance") %>%
    rbind(data.frame(name="DIC", variable="DIC", coefficients=as.numeric(jags.out$BUGSoutput$DIC), sd=NA_real_, lb=NA_real_, ub=NA_real_)) %>%
    tibble::column_to_rownames(var = "name")
  
  # Return --------------------------------------------------------------------------------------- #
  
  if(monitor==T) return(list("reg.table"=reg.table, "w"=w, "re.l1"=re.l1, "re.l3"=re.l3, "pred"=pred)) else return(list("reg.table"=reg.table, "w"=c(), "re.l1"=c(), "re.l3"=c(), "pred"=c()))
  
}
