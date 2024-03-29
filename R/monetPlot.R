#' @title Create monetPlots for your rmm results.
#'
#' @description The function \strong{rmmmonetPlot} creates a density plot of the posterior 
#' distribution of your model parameters and the traceplot that led to this density.
#' 
#' @param rmm A rmm object. rmm has to be run with monitor=T
#' @param parameter A string with the parameter name. The internal name has to be used, which are the rownames in the rmm reg.table output.
#' @param lab String to describe the parameter on the graph's x-axis. Optional. If not specified, the internal parameter name is used.
#' @param r Specify number of decimal places. Default equals 3.
#' @param sav TRUE or FALSE (default). If \code{TRUE}, the graph is saved to the current working directory as .png
#'
#' @return Returns a plot. The solid vertical is at 0 and the dashed vertical line is the mode of the posterior distributions.
#'
#' @examples data(coalgov)
#' m1 <- rmm(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/offset(n), constraint=1)) + majority + hm(id=cid, name=cname, type=RE, showFE=F),
#'           family="Weibull", monitor=T, data=coalgov)
#' monetPlot(m1, parameter="b.l1")
#'
#' @export monetPlot
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

monetPlot <- function(rmm, parameter, lab=F, r=3, sav=F) {
  
  library(ggplot2)
  library(cowplot)
  library(ggmcmc)
  
  # Retrieve mcmc list --------------------------------------------------------------------------- #
  
  if(is.null(rmm$jags.out)) stop("JAGS output could not be retrieved. monitor=T must be specified when running rmm.")
  mcmclist <- mcmcplots::as.mcmc.rjags(rmm$jags.out)
  
  # Differentiate between b.l1 and b.l1[1] ------------------------------------------------------- #
  
  parameter <- unlist(strsplit(parameter, fixed=T, split="["))
  
  if(length(parameter)>1) {
    parameter[2] <- gsub("[^0-9]", "", parameter[2]) 
    pobj <- ggs(mcmclist, family=paste0("^",parameter[1], "\\[[", parameter[2], "]\\]"))
  } else {
    pobj <- ggs(mcmclist, family=parameter) %>% filter(Parameter == parameter)
  }
  
  # Parameter statistics ------------------------------------------------------------------------- #
  
  # getmode <- function(v) {
  #   uniqv <- unique(v)
  #   uniqv[which.max(tabulate(match(v, uniqv)))]
  # }
  
  pmean <- round(mean(pobj$value), r)
  p05 <- round(min(pobj$value), r)
  p95 <- round(max(pobj$value), r)
  
  # Modify theme --------------------------------------------------------------------------------- #
  
  theme_light2 <- 
    theme_light() +
    theme(legend.position="none",
          plot.title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.text = element_blank()
          ) 
  theme_set(theme_light2)
  
  lbl <- ifelse(lab==F, parameter, lab)

  # Create plots --------------------------------------------------------------------------------- #
  
  # Density plot
  p1 <- ggs_density(pobj) + 
    labs(y = "Density") +
    geom_vline(xintercept=0) +
    geom_vline(xintercept=pmean, linetype="dashed") +
    theme_light2 +
    theme(axis.ticks=element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
   
  # Traceplot
  p2 <- ggs_traceplot(pobj, original_burnin = F) + coord_flip() +
    labs(x = paste0("Scans (", max(pobj$Chain), " chains)"), y = paste0("Parameter: ", lbl)) +
    scale_y_continuous(breaks = c(p05, 0, pmean, p95)) +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=pmean, linetype="dashed") +
    theme_light2
  
  plot <- plot_grid(p1, p2, align = 'hv', nrow = 2)
  # 2do: Try to model as facet_grid so that people can add to it as they please. plot_grid does not
  #      allow for anymore changes in form of + ...
  
  # Save?
  if(sav==T) {
    ggsave(
      paste0(gsub("[[:space:]]", "", str_squish(lbl)), ".png"),
      plot,
      width = 11,
      height = 6, 
      units='in', 
      dpi=300
    )
  }
  
  return(plot)
}
