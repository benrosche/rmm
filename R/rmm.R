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

level1adat <- data.frame(l1id=seq(1:11), l3id=1, x=rnorm(11, 0, 1))
level1bdat <- data.frame(l1id=c(10,1,2,3,4,5,6,7,8,9,10,1,9,4,6,8,3,9,9,2,4), l2id=rep(1:10, times=c(2, 3, 1, 1, 2, 2, 2, 3, 2, 3)), l3id=1, x=rnorm(21, 0, 1))
level2dat <- data.frame(l1id=c("1,10", "2,3,4", "5", "6", "7,8", "9, 10", "1,9", "4,6,8", "3,9", "9,2,4"), l2id=seq(1:10), l3id=1, x=rnorm(10, 0, 1))
level3dat <- data.frame(l3id=NA)

build_data <- function(level1dat=NA, level2dat=NA, level3dat=NA) {
  
  # 1. L12=MM, L23=H ============================================================================ #
  
  # Create level-3 data #
  level3 <<- level3dat %>% arrange(l3id)
  
  # Create level-2 data #
  level2 <<- level2dat %>%
    arrange(l2id) %>%
    mutate(nl1=lengths(gregexpr(",",l1id)) + 1) %>%
    mutate(n2=cumsum(nl1), n1=lag(n2)+1) %>%
    mutate(n1 = ifelse(row_number()==1, 1, n1)) %>%
    dplyr::select(l1id, l2id, l3id, nl1, n1, n2, x)
  
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


build_model <- function(level12=NA, level23=NA) {
  
}

estimate_model <- function() {
  
}

display_result <- function() {
  
}
