# items to load including packages, R scripts, data files etc. 

# Packages ----

if(!require(SimInf)){
  install.packages("SimInf")
  library(SimInf)
}


if(!require(lhs)){
  install.packages("lhs")
  library(lhs)
}

if(!require(parallel)){
  install.packages("parallel")
  library(parallel)
}

if(!require(pbmcapply)){
  install.packages("pbmcapply")
  library(pbmcapply)
}

if(!require(igraph)){
  install.packages("igraph")
  library(igraph)
}

if(!require(actuar)){
  install.packages("actuar")
  library(actuar)
}

if(!require(dirmult)){
  install.packages("dirmult")
  library(dirmult)
}

#library(igraph)
#library(SimInf)
#library(parallel)
#library(scales)
#library(RColorBrewer)
require(tidyverse)
#library(cowplot)

# select species ----
sp <- c("animals", "cattle", "caprine")[2]

# load functions ----

# functions to simulate movement between markets and wards
source('2. MarketWardMove2025.R')

# functions to split the wards into matrices (subnodes) and run an SIR in each subnode. 
source("nodeSIR2025SimpleMod.R")

# Ward data files import and ward selection
source("4. ImportData2025.R")

# movement parameters for the simulation of the hurdle - look at measures of ward-level network connectedness
source("5. MovementPars2025.R")

# load the parameter sets for the model 
source("6. GlobalSettings2025.R")

# load interventions 
source("7. Interventions2025.R")

# load events 
source("8. Events2025.R")


u0 <- u0.outer
u0$Sh <- u0$Sh+10
# make the vector popululation relative to vector population for now
u0$Sm <- u0$Sh*100

# split nodes into subnodes 
u0 <- as.data.frame(lapply(u0, as.integer))
mod.input <- split.nodes.sir(n.subnodes = n.subnodes, u0 = u0)

# create N matrix, which determines how internal and external transfer events
# shift indiviudals between compartments. here the only event that uses shift 
# is vaccination
N.mat <- cbind(mod.input$Ev[, 1])
# susceptibles move to recovered/removed
N.mat[substr(rownames(N.mat), 1, 2) == "Sh", ] <- 
  which(substr(rownames(N.mat), 1, 2) == "Rh") - which(substr(rownames(N.mat), 1, 2) == "Sh")
# infectious stay infectious (vax too late), and cumulative incidence should never move
N.mat[substr(rownames(N.mat), 1, 2) == "Ih", ] <- 0
# recovereds shouldn't move, as they're already effectively vaxed
N.mat[substr(rownames(N.mat), 1, 2) == "Rh", ] <- 0

# events #

events <- input.list[[1]]


networkSimpMod <- mparse(
  transitions = mod.input$transitions
  , compartments = mod.input$compartments
  , gdata = c(
    recovery  = inf.pars$recovery
    , muH  = inf.pars$muH  
    , muBirth = inf.pars$muBirth 
    , mu   = inf.pars$mu  
    , p_hv = inf.pars$p_hv 
    , p_vh = inf.pars$p_vh 
    , laying = inf.pars$laying 
    , extrinInc = inf.pars$extrinInc 
    , K = inf.pars$K 
    , Kh = inf.pars$Kh
    , development = inf.pars$development 
    , biteRate = inf.pars$biteRate  
    , firstBiteDelay =  inf.pars$firstBiteDelay 
    , pupsMort =  inf.pars$pupsMort 
    , coupling = coupling
  )
  ,u0 = mod.input$u0
  ,tspan = 1:n.days
  ,events = events
  ,E = mod.input$Ev
  ,N = N.mat
)


start <- Sys.time()
resultNetworkSimpMod <- run(model = networkSimpMod)
(trajresultNetworkSimpMod <- trajectory(resultNetworkSimpMod))
end <- Sys.time()
(start-end)
save.image("networkModSimpMod.RData")