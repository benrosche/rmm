rm(list=ls())

library(dplyr)
library(haven)
library(labelled)
library(MASS)
library(lqmm)
library(rmm)

# ================================================================================================ #
# Function to create Weibull data
# ================================================================================================ #

simWeib <- function(N, lambda, rho, LP, rateC) {
    
    # Bender et al 2005: Generating survival times to simulate Cox proportional hazards models
    # https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2059
    # https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring
    
    # Function arguments:
    # N = sample size    
    # lambda = scale parameter in h0()
    # rho = shape parameter in h0()
    # LP = linear predictor
    # rateC = rate parameter of the exponential distribution of C
    
    # Weibull latent event times
    v <- runif(N, 0, 1)
    Tlat <- (- log(v) / (lambda * exp(LP)))^(1 / rho)
    
    # censoring times
    C <- rexp(n=N, rate=rateC)
    
    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)
    
    return(cbind(survtime=time, event=status))
}

# ================================================================================================ #
# Function to create data
# ================================================================================================ #

crDat <- function(party=c(0,0,0), gov=c(1,3,3,0), country=c(3,0), weight=c(0,0,0), Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), level=1, seed=NULL) {
  
  if(!level %in% c(1,2)) stop("Level can only be 1 or 2")
  
  # Create X ------------------------------------------------------------------------------------- #
  
  # Sigma = covariance matrix of party, gov, and country random effects 
  
  if(is.null(seed)) seed <- runif(1, 0, 1000)
  set.seed(seed)
  
  # Load data 
  dat.party   <- read_dta("C:/Users/benja/OneDrive - Cornell University/GitHub/govsurvival/data/Rosche2017-party.dta")
  dat.gov     <- read_dta("C:/Users/benja/OneDrive - Cornell University/GitHub/govsurvival/data/Rosche2017-government.dta")
  dat.country <- read_dta("C:/Users/benja/OneDrive - Cornell University/GitHub/govsurvival/data/Rosche2017-country.dta")
  
  # Variable selection
  dat <- 
    dat.party %>% 
    dplyr::select(pid, gid, cid, prime, pipd, finance1, pseatrel) %>% dplyr::rename(ipd=pipd, fdep=finance1) %>%
    dplyr::inner_join(dat.gov %>% 
                        dplyr::select(gid, cid, country, gstart, gend, nPG, majority, mwc, rilegov2) %>% 
                        dplyr::rename(n=nPG, cname=country, hetero=rilegov2), by=c("gid", "cid")) %>% # add gov vars
    dplyr::inner_join(dat.country %>% dplyr::select(cid, investiture, pmpower), by=c("cid")) %>% # add country vars
    dplyr::relocate(c(cname, gstart, gend, n), .after=cid) 
  
  # Create RE
  
  if(det(Sigma)<0) {
    Sigma <- make.positive.definite(Sigma)
    message("Sigma changed")
  }
  
  re.party   <- mvrnorm(n = length(unique(dat %>% .$pid)), mu=c(0,0,0), Sigma)[,1]
  re.gov     <- mvrnorm(n = length(unique(dat %>% .$gid)), mu=c(0,0,0), Sigma)[,2]
  re.country <- mvrnorm(n = length(unique(dat %>% .$cid)), mu=c(0,0,0), Sigma)[,3]
  
  # Add RE
  dat <- 
    dat %>%
    dplyr::group_by(pid) %>% dplyr::mutate(l1id = cur_group_id()) %>% dplyr::ungroup() %>%  dplyr::mutate(re.party=re.party[l1id]) %>%
    dplyr::group_by(gid) %>% dplyr::mutate(l2id = cur_group_id()) %>% dplyr::ungroup() %>%  dplyr::mutate(re.gov=re.gov[l2id]) %>% 
    dplyr::group_by(cid) %>% dplyr::mutate(l3id = cur_group_id()) %>% dplyr::ungroup() %>%  dplyr::mutate(re.country=re.country[l3id]) %>% 
    dplyr::select(-l1id, -l2id, -l3id)
  
  # Create Y ------------------------------------------------------------------------------------- #
  
  # Create LP
  crLP <- 
    dat %>%
    dplyr::inner_join(dat.gov %>% 
                        dplyr::select(gid, cid, event_wkb, dur_wkb) %>% 
                        dplyr::rename(earlyterm=event_wkb, govdur=dur_wkb), by=c("gid", "cid")) %>%
    dplyr::mutate(partyeffect=party[1]*ipd+party[2]*fdep+party[3]*re.party) %>%
    dplyr::mutate(w=1/n^exp(-(weight[1]*pseatrel+weight[2]*hetero+weight[3]*pmpower)), ng=max(gid), w=w*ng/sum(w)) %>%
    dplyr::group_by(gid) %>%
    dplyr::mutate(aggpartyeffect=sum(w*partyeffect)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(LP=aggpartyeffect+gov[1]+gov[2]*majority+gov[3]*mwc+gov[4]*re.gov+country[1]*investiture+country[2]*re.country) %>%
    dplyr::select(-partyeffect, -aggpartyeffect, -re.party, -re.gov, -re.country) %>%
    dplyr::rename(sim_w=w)
  
  # Create DVs
  crY <-
    crLP %>%
    dplyr::mutate(sim_y=LP) %>%
    cbind(simWeib(N=dim(crLP)[1], lambda=0.01, rho=1, LP=crLP %>% .$LP, rateC=0.001)) %>% 
    dplyr::rename(sim_st=survtime, sim_e=event) %>%
    dplyr::select(-LP, -ng)

  
  # Aggregate to second level 
  if(level==2) {
    crY <- 
      crY %>%
      dplyr::group_by(gid) %>%
      dplyr::mutate(ipd=mean(ipd), fdep=mean(fdep)) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-pid, -sim_w)
  }
  
  return(crY)
}



# ================================================================================================ #
# Save data for rmm() package
# ================================================================================================ #

# coalgov <- crDat(party=c(3,3,1), gov=c(1,3,3,1), country=c(3,1), weight=c(0,0,0), Sigma=matrix(c(1,0,0, 0,1,0, 0,0,1),3,3), level=1, seed=1)
# var_label(coalgov) <- c("Unique party ID", "Unique government ID", "Unique country ID", "Country name", "Government start date", "Government end date", "# government parties", "Prime minister party", "Intra-party democracy", "Financial dependency", "Party's relative seat share within coalition", "Majority government", "Minimal winning coalition", "SD(rile) of goverment / SD(rile) of parliament", "Investiture vote", "Prime ministerial powers", "Discretionary early termination ", "Government duration", "Simulated weights", "Simulated linear outcome", "Simulated survival time", "Simulated event status")
# save(coalgov, file="data/coalgov.RData")

