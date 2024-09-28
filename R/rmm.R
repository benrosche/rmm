#' @title Bayesian multiple membership multilevel models with endogenized weights using 'JAGS'
#'
#' @description The \strong{rmm} package provides an interface to fit Bayesian multiple membership 
#' multilevel models with endogenized weights using \href{http://mcmc-jags.sourceforge.net/}{JAGS}.
#' 
#' @details 
#' The main function of \strong{rmm} is \code{rmm}, which uses formula syntax to specify a multiple membership 
#' multilevel model with endogenized weights. Based on the supplied formulas, data, and additional information, 
#' it writes code to fit the model in JAGS. Subsequently, the JAGS output is processed to ease the
#' interpretation of the model results.
#'  
#' In order to fit the models, JAGS must be \href{https://sourceforge.net/projects/mcmc-jags/files/}{installed}.
#' 
#' The package \strong{rmm} estimates models with a complex nonstandard multilevel structure, known 
#' as multiple membership multilevel structure. The difference of this package to other packages and programs 
#' to estimate multiple membership multilevel models, such as \code{\link[brms]{brms}} or \href{http://www.bristol.ac.uk/cmm/software/mlwin/}{MLwiN}, 
#' is that \strong{rmm} allows to endogenize the membership weights with a weight function. In doing so, rmm allows 
#' to examine the process by which the effects of lower-level units aggregate to a higher level (micro-macro link).
#' 
#' Accessible introductions to multiple membership models are given by the report by Fielding and Goldstein (2006) 
#' and the book chapter by Beretvas (2010). More advanced treatments of multiple membership models are  
#' provided in the multilevel textbook by Goldstein (2011, Chapter 13), the book chapters on 
#' multiple membership models by Rasbash and Browne (2001, 2008), the paper by Browne et al. (2001), and the report by Leckie (2013).
#' 
#' \bold{General formula structure}  
#' 
#' \code{Y ~ 1 + mm(id(l1id, l2id), mmc(X.L1), mmw(w ~ 1 / offset(N), constraint=1)) + X.L2 + X.L3 + hm(id=l3id, name=l3name, type=FE, showFE=F)} 
#' \itemize{
#'   \item Dependent variable: \strong{Y}
#'   \item Multiple membership object: \strong{mm()} to analyze how the effects of level-1 predictors from multiple constituting members aggregate to level 2
#'   \item Level-2 predictors: \strong{X.L2}, being something like X1 + ... + XN 
#'   \item Level-3 predictors: \strong{X.L3}, being something like X1 + ... + XN 
#'   \item Hierarchical membership object: \strong{hm()} to recognize that level 2 units embedded in a third level
#' }
#' \bold{Currently supported dependent variables / link functions}  
#' 
#' \itemize{
#'   \item Gaussian continuous variable \code{Y}
#'   \item Binomial outcome for logistic regression \code{Y}
#'   \item Conditional logistic outcomes \code{???}
#'   \item Weibull survival time: \code{Surv(survivaltime, event)} or \code{Surv(survivaltime, event, upperlimit)} (where upperlimit is the upper limit for predictions(!))
#'   \item Cox survival time: \code{Surv(survivaltime, event)}
#' }
#' 
#' \bold{Vector of level-2 predictors}
#' 
#' An intercept at the main level 2 is added whether or not a \code{1} is specified in the beginning. Interaction terms have to be included as separate variables. 
#' 
#' \bold{Multiple membership object mm()}
#' 
#' \itemize{
#'   \item \code{id()} to indicate level-1 and level-2 ids
#'   \item \code{mmc()} to specify level-1 predictors. No intercept allowed. Interaction terms have to be included as separate variables. 
#'   \item \code{mmw()} to specify the weight function (micro-macro link). The function can be nonlinear and contain variables but needs to be identifiable. 
#'         To give a few examples: \code{w ~ 1/offset(N)} specifies mean-aggregation, with N being a variables that indicates the number of level-1 units per level-2 entity.
#'         If no mmw() is specified, \code{w ~ 1/offset(N)} is assumed. In Rosche (2021), I propose to use \code{w ~ 1/offset(N)^exp(-(X.W))} as the general form for weight functions. 
#'         This function ensures that weights are limited by 0 and 1. Specifying variables as \code{offset(X.W)} will not estimate a parameter for this variable. 
#'         
#'         Two different identification restrictions are provided:
#'   \itemize{
#'     \item \code{mmw(w ~ ..., constraint=1)}: \code{constraint=1} restricts the weights to sum to 1 for each level-2 entity. (default)
#'     \item \code{mmw(w ~ ..., constraint=2)}: \code{constraint=2} restricts the weights to sum to the total number of level-2 entities over the whole dataset, allowing some level-2 entities to have weights smaller/larger than 1.
#'   }
#'   
#'   \itemize{
#'     \item \code{mmw(w ~ ..., ar=TRUE)}: Allows random effects of level-1 units to change across memberships in level-2 entities.
#'     \item \code{mmw(w ~ ..., ar=FALSE)}: Assumes all level-1 units to have one random effect
#'   }
#' }
#' 
#' \bold{Hierachical membership object hm()}
#' 
#' \itemize{
#'   \item \code{id=l3id} to indicate level-3 id
#'   \item \code{name=l3name} to specify value labels for level 3 units. 
#'   \item \code{type=RE} (default) or \code{type=FE} to choose between random- or fixed effect estimation. 
#'         If RE is chosen, level 3 predictors can be added. If FE is chosen, each level 3 unit has its own intercept and level 3 predictors are removed. 
#'         If \code{showFE=TRUE} the fixed effects are reported, otherwise omitted (default). The first l3id is the base.
#' }
#' 
#' \bold{More details on changing priors}
#' 
#' Priors of the following parameters may be changed: \code{b.l1, b.l2, b.l3, b.w, tau.l1, tau.l2, tau.l3}. 
#' The priors are specified as a list with parameter names as tags and their prior specification as values:
#' \code{priors=list("b.l1"="dnorm(0,0.01)")}. In this example, the priors of all level-1 regression coefficients 
#' are changed to a more informative prior that has a smaller variance than the default (\code{dnorm(0,0.0001)}). 
#' I refer to the \href{https://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf}{JAGS manual} 
#' for more details on possible prior specifications. 
#' 
#' \bold{More details on the weight function}
#' 
#' ...
#'
#' \bold{More details on constructing the data}
#' 
#' ...
#' 
#' \bold{Tips}
#' 
#' \itemize{
#' \item \code{Error in update.jags(model, n.iter, ...) : Error in node w[1285] Invalid parent values} The weight function must be designed such that the distribution of weights is in 
#' line with the priors for all other parameters. This error could, for instance, be caused if weights can be negative but negative weights cause the distribution of other parameters to be outside of the distribution of their priors.
#' Carefully designing the weight function so that it is properly bounded may therefore help. Specifying \code{transform="std"} may help as well. 
#' \item Including weight regressors demands a lot from your data It is therefore a good idea to start with slightly more informative priors. 
#' I suggest starting with \code{priors = list("b.w"="dnorm(0,0.1)"} and then increasing the variance step by step.
#' }
#' ...
#' @param formula A symbolic description of the model in form of an R formula. More details below.
#' @param family Character vector. Currently supported are "Gaussian", "Binomial", "Weibull", and "Cox". Not yet implemented: "CondLogit"
#' @param priors A list with parameter names as tags and their prior specification as values. More details below.
#' @param inits A list with parameter as tags and their initial values as values. This list will be used in all chains. If NULL, JAGS and rmm select appropriate inits.
#' @param n.iter Total number of iterations.
#' @param n.burnin Number of iterations that will be discarded.
#' @param n.thin Thinning rate.
#' @param chains Number of chains.
#' @param seed A random number.
#' @param run A logical value (True or False) indicating whether JAGS should estimate the model.
#' @param monitor A logical value (True or False). If \code{True}, weights, random effects, predictions, and JAGS output is saved as well.
#' @param transform Character vector or FALSE. Specifying \code{center} or \code{std} to center or standardize continuous predictors before estimation. Specifying \code{std2} will divide by two times the standard deviation, so that regression coefficients are comparable to those of binary predictors (Gelman 2008). 
#' @param modelfile Character vector or TRUE|False. If TRUE, the JAGS model is saved in rmm/temp/modelstring.txt. If a file path is supplied as string, rmm will just create the data structure and use the provided modelfile. Run \code{.libPaths()} to see where R packages are stored.
#' @param data Dataframe object. The dataset must have level 1 as unit of analysis. More details below.
#'
#' @return A list with 7 elements: reg.table, w, re.l1, re.l3, pred, input, jags.out. If monitor=F, only the 
#'         regression table is returned. If monitor=T, the predicted weights, level-1 random effects (if specified in the model),
#'         level-3 random effects (if specified in the model), predicted values of the dependent variable, and the internally created variables are returned. 
#'         The last element of the list is the unformatted Jags output. 
#' @examples data(coalgov)
#' m1 <- rmm(Surv(govdur, earlyterm, govmaxdur) ~ 1 + mm(id(pid, gid), mmc(fdep), mmw(w ~ 1/offset(n), constraint=1)) + majority + hm(id=cid, name=cname, type=RE, showFE=F),
#'           family="Weibull", monitor=T, data=coalgov)
#' m1$reg.table # the regression output
#' m1$w         # the estimated weights
#' m1$re.l1     # the level-1 random effects
#' m1$re.l3     # the level-3 random effects
#' m1$pred      # posterior predictions of the dependent variable (linear predictor for \code{family="Gaussian"}, survival time for \code{family="Weibull"})
#' m1$input     # internal variables
#' jags.out <- m1$jags.out # JAGS output
#' monetPlot(m1, "b.l1") # monetPlot to inspect the posterior distribution of the model parameters
#'
#' @export rmm
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#' @references 
#' Rosche, B. (2021). A multilevel model for coalition governments: Uncovering dependencies within and across governments due to parties. https://doi.org/10.31235/osf.io/4bafr

rmm <- function(formula, family="Gaussian", priors=NULL, inits=NULL, n.iter = 1000, n.burnin = 500, n.thin = max(1, floor((n.iter - n.burnin) / 1000)), chains=3, seed=NULL, run=T, parallel=F, monitor=T, transform="center", modelfile=F, data=NULL) {

  # formula = Surv(govdur, earlyterm) ~ 1 + majority + mwc; family = "Weibull"; priors=NULL; inits=NULL; n.iter=100; n.burnin=10; n.thin = max(1, floor((n.iter - n.burnin) / 1000)); chains = 3; seed = 123; run = T; parallel = F; monitor = T; transform = "center"; modelfile = T; data = coalgov
  # source("./R/dissectFormula.R"); source("./R/createData.R"); source("./R/editModelstring.R"); source("./R/createJagsVars.R"); source("./R/formatJags.R"); 
  
  # ---------------------------------------------------------------------------------------------- #
  # 0. Checks
  # ---------------------------------------------------------------------------------------------- #
  
  if(is.null(data)) stop("No data supplied.")

  # ---------------------------------------------------------------------------------------------- #
  # 1. Dissect formula 
  # ---------------------------------------------------------------------------------------------- #
  
  DIR <- system.file(package = "rmm")
  
  c(ids, vars, l1, l3) %<-% dissectFormula(data, family, formula)
  
  # ---------------------------------------------------------------------------------------------- #
  # 2. Disentangle vars and data into l1-3
  # ---------------------------------------------------------------------------------------------- #
  
  c(data, level1, level2, level3, weightf) %<-% createData(data, ids, vars, l1, l3, transform)
  
  # Remove varlist
  rm(vars)

  # ---------------------------------------------------------------------------------------------- #
  # 3. Edit modelstring 
  # ---------------------------------------------------------------------------------------------- #
  
  modelstring <- editModelstring(family, priors, l1, l3, level1, level2, level3, DIR, monitor, modelfile)

  # ---------------------------------------------------------------------------------------------- #
  # 4. Transform data into JAGS format
  # ---------------------------------------------------------------------------------------------- #
  
  c(ids, Ns, Xs, Ys, jags.params, jags.inits, jags.data) %<-% createJagsVars(data, family, level1, level2, level3, weightf, ids, l1, l3, monitor, modelfile, chains, inits)
  
  list2env(c(ids, Ns, Xs, Ys), envir=environment())
  
  # ---------------------------------------------------------------------------------------------- #
  # 5. Run JAGS 
  # ---------------------------------------------------------------------------------------------- #
  
  if(run==T) {
    
    # Get seed
    if(is.null(seed)) seed <- round(runif(1, 0, 1000)) 
    
    if(parallel) {
      
      # Run parallel ----------------------------------------------------------------------------- #
      
      readr::write_file(modelstring, paste0(DIR, "/temp/jags-parallel.txt"))
      jags.out <- do.call(jags.parallel, list(data=jags.data, inits = jags.inits[1], n.chains = chains, parameters.to.save = jags.params, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, jags.seed = seed, model.file = paste0(DIR, "/temp/jags-parallel.txt")))
      file.remove(paste0(DIR, "/temp/jags-parallel.txt"))
      # Three peculiarities about jags.parallel:
      # - It cannot read the model from textConnection(modelstring)
      # - It cannot read variables from the global environment - do.call needs to be used 
      # - There seems to be a bug in that it wants just one list element of inits instead of n.chains number of list elements
      
    } else {
      
      # Run sequentially ------------------------------------------------------------------------- #
      
      set.seed(seed) 
      jags.out <- jags(data=jags.data, inits = jags.inits, n.chains = chains, parameters.to.save = jags.params, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = textConnection(modelstring)) 
      
    } 
   
    # Format JAGS output ------------------------------------------------------------------------- #
    
    c(reg.table, w, re.l1, re.l3, pred) %<-% formatJags(jags.out, monitor, Ns, l1, l3, level1, level2, level3, weightf) 
    
    # Prepare return ----------------------------------------------------------------------------- #
    
    # Save info on transformed vars
    isTransformed <- function(df) {
      tvars <- if(dim(df)[2]>0) sapply(names(df), function(x) mean(df[[x]])<0.0001) else F
      return(names(tvars)[tvars==T])
    }
  
    transformedVars <- c(if(!is.null(level1[["dat"]])) isTransformed(level1 %>% .$dat %>% dplyr::select(!!level1[["vars"]])),
                         if(!is.null(level2[["dat"]])) isTransformed(level2 %>% .$dat %>% dplyr::select(!!level2[["vars"]])),
                         if(!is.null(level3[["dat"]])) isTransformed(level3 %>% .$dat %>% dplyr::select(!!level3[["vars"]])),
                         if(!is.null(weightf[["dat"]])) isTransformed(weightf %>% .$dat %>% dplyr::select(!!weightf[["vars"]])))
    
    # Save info on input
    input <- 
      if(isTRUE(monitor)) {
        append(
          list(
            "family"=family, "priors"=priors, "inits"=inits, 
            "n.iter"=n.iter, "n.burnin"=n.burnin, "n.thin"=n.thin, "chains"=chains, "seed"=seed, "run"=run, "parallel"=parallel, 
            "monitor"=monitor, "transform"=transform, "modelfile"=modelfile,
            "lhs" = level2$lhs, "l1vars"=level1$vars, "l2vars"=level2$vars, "l3vars"=level3$vars, "transformedVars"=transformedVars,
            "n.ul1"=Ns$n.ul1, "n.l1"=Ns$n.l1, "n.l2"=Ns$n.l2, "n.l3"=Ns$n.l3
          ), 
          c(l1, l3)
        )
      } else c()  
    
    # Return ------------------------------------------------------------------------------------- #
    
    out <- list("reg.table"=reg.table, "w"=w, "re.l1"=re.l1, "re.l3"=re.l3, "pred"=pred, "input"=input, "jags.out"=if(isTRUE(monitor)) jags.out else c())
    
    class(out) <- "rmm"
    
    return(out)
    
  } else {
    
    message("Data and model have been created without any errors.")
    
  }
  
}
