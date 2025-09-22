#' @title Create monetPlots for your rmm results.
#'
#' @description The function \strong{monetPlot} creates a density plot of the posterior 
#' distribution of your model parameters and the traceplot that led to this density.
#' 
#' @param rmm A rmm object. rmm has to be run with monitor=T
#' @param parameter A string with the parameter name. The internal name has to be used, which are the rownames in the rmm reg.table output.
#' @param centrality A string specifying one of the following options: "median", "mean", "MAP", or "mode".
#' @param lab String to describe the parameter on the graph's x-axis. Optional. If not specified, the internal parameter name is used.
#' @param r Specify number of decimal places. Default equals 3.
#' @param sav TRUE or FALSE (default). If \code{TRUE}, the graph is saved to the current working directory as .png
#'
#' @return Returns a plot. The solid vertical is at 0 and the dashed vertical line is the mode of the posterior distributions.
#'
#' @examples data(coalgov)
#' m1 <- rmm(Surv(govdur, earlyterm) ~ 1 + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/n, constraint=1)) + majority + hm(id=cid, name=cname, type=RE, showFE=F),
#'           family="Weibull", monitor=T, data=coalgov)
#' monetPlot(m1, parameter="b.l1")
#'
#' @export monetPlot
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

monetPlot <- function(rmm, parameter, centrality="median", lab=F, r=3, sav=F) {
  
  library(ggplot2)
  library(ggmcmc)
  library(coda)
  library(patchwork)
  library(bayestestR)
  
  # Get mcmclist and posterior stats ------------------------------------------------------------- #
  
  if(is.null(rmm$jags.out)) stop("JAGS output could not be retrieved. monitor=T must be specified when running rmm.")
  
  mcmclist <- 
    rmm$jags.out %>% 
    as.mcmc() %>% 
    .[, parameter, drop = FALSE]
  
  param <- describe_posterior(mcmclist, centrality=centrality, ci=1) 
  
  # Modify theme --------------------------------------------------------------------------------- #
  
  theme_light2 <- 
    theme_light() +
    theme(
      text = element_text(size = 16),
      legend.position="none",
      plot.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text = element_blank()
    ) 
  theme_set(theme_light2)
  
  # Create plots --------------------------------------------------------------------------------- #
  
  mcmc_ggs <- ggs(mcmclist)
  
  # Density plot
  p1 <- 
    ggs_density(mcmc_ggs) + 
    geom_vline(xintercept=0) +
    geom_vline(xintercept=round(param[,2], r), linetype="dashed") +
    theme_light2 +
    theme(
      axis.ticks=element_blank(), 
      axis.text.x=element_blank(), 
      axis.text.y=element_blank()
    ) +
    labs(x=NULL, y = "Density", subtitle = paste0("Parameter: ", parameter, ifelse(!isFALSE(lab), paste0(" - ", lab), "")))
  
  # Traceplot
  p2 <- 
    ggs_traceplot(mcmc_ggs, original_burnin = F) + 
    geom_hline(yintercept=0) +
    geom_hline(yintercept=round(param[,2], r), linetype="dashed") +
    coord_flip() +
    scale_y_continuous(breaks = c(0, round(param[,2], r), round(param[,4], r), round(param[,5], r)) %>% sort()) +
    theme_light2 +
    labs(x = paste0("Iteration (", max(mcmc_ggs$Chain), " chains)"), y = NULL) 
  
  p <- 
    (p1 / p2) +
    plot_layout(heights = c(1, 1)) +  
    plot_annotation(theme = theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) 
  
  
  # Save?
  if(sav==T) {
    ggsave(
      paste0(gsub("[[:space:]]", "", str_squish(parameter)), ".png"),
      p,
      width = 11,
      height = 6, 
      units='in', 
      dpi=300
    )
  }
  
  return(p)
}

# NEW VERSION (INTEGRATE!)
# monetPlot <- function(g.ergm, param, paramlabel = NULL, breaks = waiver(), yaxis = T) {
#   
#   if (is.null(paramlabel)) paramlabel = param
#   
#   if (is.null(g.ergm$sample)) stop("Not estimated by MCMC ...")
#   mcmclist <- g.ergm$sample %>% as.mcmc.list()
#   
#   ggs.obj <- ggs(mcmclist, family = paste0("A", escapeRegex(param), "$"), 
#                  par_labels = data.frame(Parameter = param, Label = paramlabel))
#   
#   p.median <- ggs.obj %>% pull(value) %>% median()
#   
#   p1 <-
#     ggs_density(ggs.obj, hpd = TRUE) +
#     geom_vline(xintercept = 0) +
#     geom_vline(xintercept = p.median, linetype = "dashed") +
#     scale_color_viridis_d(direction = -1) +
#     scale_fill_viridis_d(direction = -1) +
#     labs(y = if (isTRUE(yaxis)) "Density" else "") +
#     theme_minimal() +
#     theme(
#       plot.margin = unit(c(0, 0, 0, 0), "cm"),
#       legend.position = "none",
#       panel.grid = element_blank(),
#       axis.ticks = element_blank(),
#       axis.title.x = element_blank(),
#       axis.text = element_blank(),
#       strip.text = element_text(size=15, face-"bold", color="black"),
#       strip.background = element_blank()
#     )
#   
#   p2 <-
#     ggs_traceplot(ggs.obj, original_burnin = FALSE) +
#     geom_hline(yintercept = 0) +
#     geom_hline(yintercept = p.median, linetype = "dashed") +
#     coord_flip() +
#     scale_y_continuous(breaks = breaks) +
#     scale_color_viridis_d(direction = -1) +
#     labs(
#       x = if (isTRUE(yaxis)) paste0("Scans (", ggs.obj %>% pull(Chain) %>% max(), " chains)") else "",
#       y = ""
#     ) +
#     theme_minimal() +
#     theme(
#       plot.margin = unit(c(-0.6, 0, 0, 0), "cm"),
#       legend.position = "none",
#       panel.grid = element_blank(),
#       axis.text.y = element_blank(),
#       strip.text = element_blank()
#     )
#   
#   return(plot_grid(p1, p2, align = "v", nrow = 2))
# }

