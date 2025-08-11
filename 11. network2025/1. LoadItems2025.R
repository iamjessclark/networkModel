# items to load including packages, R scripts, data files etc. 

# Packages ----

if(!require(SimInf)){
  install.packages("SimInf")
  library(SimInf) <- 
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
source("3. NodeSIRS2025.R")

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

