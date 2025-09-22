#' @title summary() method for an rmm object
#'
#' @description summary() method for an rmm object
#' 
#' @return Returns a table with regression results.
#'
#' @examples data(coalgov)
#' m1 <- rmm(Surv(govdur, earlyterm, govmaxdur) ~ 1 + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/n, constraint=1)) + majority + hm(id=cid, name=cname, type=RE, showFE=F),
#'           family="Weibull", monitor=T, data=coalgov)
#' summary(m1)
#'
#' @exportS3Method summary rmm
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

summary.rmm <- function(rmm, r=3) {
  
  return( rmm$reg.table %>% dplyr::mutate(across(where(is.numeric), ~round(.,r))) )
  
}




