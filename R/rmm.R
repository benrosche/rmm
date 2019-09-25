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

level1 <- level1 %>% rename(pid=l1id, gid=l2id) %>% mutate(survival=rnorm(21, 0, 1))



rmm <- function(formula, data, family) {
  
  # Formula
  formula <- formula(survival ~ var1 + mm(id(pid, gid), mmc(var3, var4), mmw(pupils^exp(teacher*b))))
  f.chr <- sort(el(strsplit(as.character(formula)[3], " \\+ ")))[1]
  
  mmc <- gsub(".*mmc\\((.*?)\\).*", "\\1", f.chr)
  mmw <- gsub(".*mmw\\((.*?\\))\\).*", "\\1", f.chr)

  # Data
  data <- level1
  
  formula.vars <- all.vars(formula)
  
  # IDs
  id <- str_trim(str_split(gsub(".*id\\((.*?)\\).*", "\\1", f.chr), ",")[[1]])
  l1id <- data %>% rename(l1id = !!id[1]) %>% .$l1id
  l2id <- data %>% rename(l2id = !!id[2]) %>% .$l2id
  
  # Y
  Y <- data %>% rename(Y = !!formula.vars[1]) %>% .$Y
  
  # X
  
  data %>% group_by(l2id)
  
  CONTINUE HERE
  
  level2 <<- level2dat %>%
    dplyr::select(l1id, l2id, l3id) %>%
    arrange(l2id) %>%
    mutate(nl1pl2=lengths(gregexpr("[0-9]+", l1id))) %>%
    mutate(l1i2=cumsum(nl1pl2), l1i1=lag(l1i2)+1) %>%
    mutate(l1i1 = ifelse(row_number()==1, 1, l1i1)) %>%
    dplyr::select(l1i1, l1i2, l2id, l3id, nl1pl2) %>%
    bind_cols(X.l2)
  

  
  Y <- quote(survival)
  X <- quote(majority)
  mm <- formula(pid)
  mmc <- quote(aipd)
  mmw <- quote(1/N)


}




estimate_model <- function(burnin=100, iter=1000, seed=1) {
  
  # IDs
  l1id <- level1 %>% .$l1id
  l2id <- level2 %>% .$l2id
  l3id <- level2 %>% .$l3id
  
  l1i1 <- level2 %>% .$l1i1
  l1i2 <- level2 %>% .$l1i2
  
  # Y
  y.l1 <- level1 %>% .$y
  y.l2 <- level2 %>% .$y
  y.l3 <- level3 %>% .$y
  
  if (!is.null(y.l1)) {
    
    Y <- level1 %>% .$y
    level1 <- level1 %>% dplyr::select(-y)
    
  } else if (!is.null(y.l2)) {
    
    Y <- level2 %>% .$y
    level2 <- level2 %>% dplyr::select(-y)
    
  } else if(!is.null(y.l3)) {
    
    Y <- level3 %>% .$y
    level3 <- level3 %>% dplyr::select(-y)
    
  }
  
  # Xs
  X.l1 <- as.matrix(level1 %>% select(-l1id, -l2id, -l3id))
  X.l2 <- as.matrix(level2 %>% select(-l1i1, -l1i2, -l2id, -l3id, -nl1pl2))
  X.l3 <- as.matrix(level3 %>% select(-l3id))
  
  nl1pl2 <- level2 %>% .$nl1pl2
  
  # Ns
  n.l1 <- length(l1id)
  n.ul1 <- length(unique(l1id))
  n.l2 <- length(l2id)
  n.l3 <- length(l3id)
  n.ul3 <- length(unique(l3id))
  n.Xl1 <- dim(X.l1)[2]
  n.Xl2 <- dim(X.l2)[2]
  n.Xl3 <- dim(X.l3)[2]
  
  # JAGS
  jags.data   <- c("l1id", "l1i1", "l1i2", "l2id", "l3id", "Y", "X.l1", "X.l2", "X.l3", "nl1pl2", "n.l1", "n.ul1", "n.l2", "n.l3", "n.ul3", "n.Xl1", "n.Xl2", "n.Xl3") 
  jags.params <- c("b.l1", "b.l2", "b.l3", "ppp.bl1", "ppp.bl2", "ppp.bl3", "sigma.l1", "sigma.l2", "sigma.l3")
  jags.inits <- list( 
    list(".RNG.seed" = seed),
    list(".RNG.seed" = seed+1),
    list(".RNG.seed" = seed+2))
  
  jags(data=jags.data, inits = jags.inits, parameters.to.save = jags.params, n.chains = 3, n.iter = iter, n.burnin = burnin, model.file = modelstring)
}

display_result <- function() {
  
}
