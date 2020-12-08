# ================================================================================================ #
# Function formatJags
# ================================================================================================ #

formatJags <- function(jags.out, hdi, r, monitor, vars, Ns, mm, l3, level3) {

  # Unpack lists --------------------------------------------------------------------------------- #

  l1vars <- vars$l1vars
  l2vars <- vars$l2vars
  l3vars <- level3$vars
  l3dat  <- level3$dat
  
  l3name <- l3$l3name
  l3type <- l3$l3type
  
  n.ul1 <- Ns$n.ul1
  l1n  <- Ns$l1n
  n.l2 <- Ns$n.l2
  n.GPN <- Ns$n.GPN
  
  mmwar <- mm$mmwar
  lwvars <- vars$lwvars
  offsetvars <- vars$offsetvars
  
  # Create reg.table from JAGS output ------------------------------------------------------------ #
  
  reg.table <- 
    as.data.frame(jags.out$BUGSoutput$summary[, c(1, 2, 3, 7)]) %>% 
    tibble::rownames_to_column(var="name") %>% 
    dplyr::rename(estimate=2, sd=3, lb=4, ub=5)
  
  pppvalues <- # extract ppp here because they should remain a mean estimate 
    reg.table %>%  
    dplyr::filter(startsWith(name, "ppp.")) %>% 
    dplyr::mutate_at(.vars=c("name"), .funs=function(x) str_remove(x, "ppp.")) %>%
    dplyr::select(-sd, -lb, -ub)
  
  # 2do: ppp values for variance terms
  
  # Mean and 95% CI or Mode and HDI -------------------------------------------------------------- #
  
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
  
  # Remove random effects, weights, and ppp from reg.table and save separately ------------------- #
  
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
      
      re.l1 <- remat
      
    } else { # AR == F
      re.l1 <- re.l1 %>% .$estimate
    }
    
    reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "re.l1["))
    
    # Level-3 RE
    re.l3 <- reg.table %>% dplyr::filter(startsWith(name, "re.l3")) %>% dplyr::mutate(estimate = round(as.numeric(estimate), r)) %>% .$estimate
    reg.table <- reg.table %>% dplyr::filter(!startsWith(name, "re.l3[")) 
    
    # Weights
    w <- reg.table %>% dplyr::filter(startsWith(name, "w")) %>% dplyr::mutate(estimate = round(as.numeric(estimate), r)) %>% .$estimate
    
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
    
  }
  
  # PPP values (add as separate column, must be after HDI) 
  reg.table <- 
    reg.table %>% 
    dplyr::left_join(pppvalues %>% dplyr::select(name, estimate) %>% dplyr::rename(ppp=estimate), by=c("name")) %>%
    dplyr::mutate(ppp=case_when(estimate>0 ~ 1-ppp, # for positive estimates, the correct p-value = 1-p
                                TRUE ~ ppp)) %>%
    dplyr::filter(!startsWith(name, "ppp.")) 
  
  # Rename --------------------------------------------------------------------------------------- #
  
  newnames <- reg.table %>% .$name
  newnames[stringr::str_detect(newnames, "b.l1")] <- l1vars
  newnames[stringr::str_detect(newnames, "b.l2")] <- l2vars
  newnames[stringr::str_detect(newnames, "b.l3")] <- if(!is.null(l3name) & l3type=="FE") {l3dat %>% .$l3name}[-1] else l3vars
  newnames[stringr::str_detect(newnames, "b.w")]  <- lwvars[!lwvars %in% offsetvars]
  
  reg.table <- 
    reg.table %>% 
    rename(coefficients=estimate) %>% # compatibility with lm() 
    dplyr::mutate(variable=newnames) %>% relocate(variable, .before = coefficients) %>%
    filter(!variable=="deviance") %>%
    rbind(c("DIC", "DIC", jags.out$BUGSoutput$DIC, NA, NA, NA, NA)) %>%
    dplyr::mutate_at(3:7, list(~round(as.numeric(.), r))) %>%
    tibble::column_to_rownames(var = "name")
  
  # Return --------------------------------------------------------------------------------------- #
  
  if(monitor==T) return(list("reg.table"=reg.table, "w"=w, "re.l1"=re.l1, "re.l3"=re.l3)) else return(list("reg.table"=reg.table, "w"=c(), "re.l1"=c(), "re.l3"=c()))
  
}
