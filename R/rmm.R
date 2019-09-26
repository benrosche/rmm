# Rmm
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

rm(list=ls())

library(tidyr)
library(dplyr)
library(R2jags)
library(stringr)

# create data =================================================================================== #

level1adat <- data.frame(l1id=seq(1:10), l3id=1, x=rnorm(10, 0, 1))
level1bdat <- data.frame(l1id=c(10,1,2,3,4,5,6,7,8,9,10,1,9,4,6,8,3,9,9,2,4), l2id=rep(1:10, times=c(2, 3, 1, 1, 2, 2, 2, 3, 2, 3)), l3id=1, x=rnorm(21, 0, 1))
level2dat <- data.frame(l1id=c("1,10", "2,3,4", "5", "6", "7,8", "9, 10", "1,9", "4,6,8", "3,9", "9,2,4"), l2id=seq(1:10), l3id=1, x=rnorm(10, 0, 1), y=rnorm(10, 0, 1))
level3dat <- data.frame(l3id=1, x=rnorm(1, 0, 1))

build_data <- function(level1dat=NA, level2dat=NA, level3dat=NA) {
  
  # 1. L12=MM, L23=H ============================================================================ #
  
  # Create level-3 data #
  level3 <<- level3dat %>% arrange(l3id)
  
  # Create level-2 data #
  X.l2 <- level2dat %>% dplyr::select(-l1id, -l2id, -l3id)
  
  level2 <<- level2dat %>%
    dplyr::select(l1id, l2id, l3id) %>%
    arrange(l2id) %>%
    mutate(nl1pl2=lengths(gregexpr("[0-9]+", l1id))) %>%
    mutate(l1i2=cumsum(nl1pl2), l1i1=lag(l1i2)+1) %>%
    mutate(l1i1 = ifelse(row_number()==1, 1, l1i1)) %>%
    dplyr::select(l1i1, l1i2, l2id, l3id, nl1pl2) %>%
    bind_cols(X.l2)
  
  # Test whether level1 units have stable features #
  if(all(level1dat %>% group_by(l1id) %>% mutate(n=n()) %>% .$n == 1) == TRUE) {
    
    # If so, create level1 in l1/l2 format #
    
    level1 <<- level2dat %>%
      arrange(l2id) %>%
      separate_rows(l1id) %>%
      mutate(l1id=as.numeric(l1id)) %>% 
      arrange(l2id, l1id) %>% 
      dplyr::select(l1id, l2id) %>%
      left_join(level1dat, by="l1id")
    
  } else {
    
    # Data is already in l1/l2 format #
    # Test whether number of level1 units provided in level1dat and level2dat are the same
    test_order <- 
      level2dat %>% dplyr::select(l2id, l1id) %>% separate_rows(l1id) %>% mutate(l1id=as.numeric(l1id)) %>% arrange(l2id, l1id) == 
      level1dat %>% dplyr::select(l2id, l1id) %>% arrange(l2id, l1id) 
    order_false <- which(test_order[,2]==FALSE)
    if(length(order_false>0)) stop("The number of level1 units provided in level1dat and level2dat are not the same. Rows in level1dat affected: ", paste(order_false, collapse=", "))
    
    level1 <<- level1dat %>% arrange(l2id, l1id)
    
  }

}

build_data(level1dat=level1bdat, level2dat=level2dat, level3dat=level3dat)

rm(level1adat, level1bdat, level2dat, level3dat)

data <- level1 %>% rename(pid=l1id, gid=l2id) %>% mutate(survival=rnorm(21, 0, 1), var1=rnorm(21, 0, 1), var2=rnorm(21, 0, 1), var3=rnorm(21, 0, 1), var4=rnorm(21, 0, 1))


# rmm ============================================================================================ #

rmm <- function(formula, data, family, seed=1, burnin=100, iter=1000) {
  
  formula <- formula(survival ~ var1 + var2 + var3 + mm(id(pid, gid), mmc(var3, var4), mmw(pupils^exp(teacher*b))))
  
  # Dissect formula 
  f.chr <- sort(el(strsplit(as.character(formula)[3], " \\+ ")))
  
  id <- str_trim(str_split(gsub(".*id\\((.*?)\\).*", "\\1", f.chr[1]), ",")[[1]])
  mmc <- str_trim(str_split(gsub(".*mmc\\((.*?)\\).*", "\\1", f.chr[1]), ",")[[1]])
  mmw <- gsub(".*mmw\\((.*?\\))\\).*", "\\1", f.chr[1])
  
  # Rename IDs
  data <- data %>% rename(l1id = !!id[1], l2id = !!id[2]) 
  
  # Disentangle level1 and level2 data
  
  # Level 1
  level1 <-
    data %>% 
    arrange(l2id) %>%
    dplyr::select(l1id, l2id, !!mmc)
  
  # Level 2
  level2 <- 
    data %>% 
    group_by(l2id) %>%
    add_count(name="l1n") %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    arrange(l2id) %>%
    mutate(l1i2=cumsum(l1n), l1i1=lag(l1i2)+1) %>%
    mutate(l1i1 = ifelse(row_number()==1, 1, l1i1)) %>%
    rename(Y = !!formula[[2]]) %>%
    dplyr::select(l1i1, l1i2, l2id, l1n, Y, !!f.chr[-1])
  
  # Prepare data for JAGS
  
  # IDs
  l1id <- level1 %>% .$l1id
  l2id <- level2 %>% .$l2id

  l1i1 <- level2 %>% .$l1i1
  l1i2 <- level2 %>% .$l1i2
  
  # Y
  Y <- level2 %>%  .$Y
  
  # Xs
  X.l1 <- as.matrix(level1 %>% select(-l1id, -l2id))
  X.l2 <- as.matrix(level2 %>% select(-l1i1, -l1i2, -l2id, -l1n))
  l1n <- level2 %>% .$l1n
  
  # Ns
  n.l1 <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2 <- length(l2id)
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  
  # JAGS
  jags.data   <- c("l1id", "l1i1", "l1i2", "Y", "X.l1", "X.l2", "l1n", "n.l1", "n.ul1", "n.l2", "n.Xl1", "n.Xl2") 
  jags.params <- c("b.l1", "b.l2", "ppp.bl1", "ppp.bl2", "sigma.l1", "sigma.l2")
  jags.inits <- list( 
    list(".RNG.seed" = seed),
    list(".RNG.seed" = seed+1),
    list(".RNG.seed" = seed+2))
  
  jags(data=jags.data, inits = jags.inits, parameters.to.save = jags.params, n.chains = 3, n.iter = iter, n.burnin = burnin, model.file = modelstring)

}






display_result <- function() {
  
}
