
# movement multiplier ----
# optionally prevent any movement, either via markets or across ward borders
# set to 0 to prevent all movement, set to 1 for usual simulation behaviour
movt.multiplier <- 1

# how many subnodes to split the wards into, to simulate heterogeneity within wards ----
# subnodes form a square grid 
n.subnodes <- 2^2 # 1 means don't split # 3^2 means 3x3 grid, etc

# strength of coupling between subnodes (used if n.subnodes > 1) ----
coupling <- 0.02   
#coupling <- 0.0001

# SEIR model parameters ----
inf.pars <- 
  list( recovery  = 0.2
        , muH  = 0.000456621
        , muBirth = 0.001369863
        , mu   = 0.03743286
        , p_hv = 0.7082972
        , p_vh = 0.2618384
        , laying = 20000
        , extrinInc = 0.08136993
        , K = 80000000
        , Kh = 100000
        , development = 0.1
        , biteRate = 0.515714
        , firstBiteDelay = 0.3840011
        , pupsMort = 0.2079803
       ,cov = 0.7
       ,i0 = 0
  )
# exclude seed ward from intervention ----
exclude.seed <- FALSE   

# reporting threshold ----
rep.thresh <- 0.005

# number of animals moving between contiguous wards ----
contig.move.n <- 
  0.25/n.subnodes * coupling * mean(wards[, paste0("pop.", sp)]) * n.days/12 * inf.pars$recovery * movt.multiplier

# allow mosquito reproduction to vary with NDVI ----
# in the previous iteration, paul used 1-ndvi but for mosquito this doesn't make sense 
# higher NDVI relates to lush green which is where more mosquitos would be found
# as opposed to closer to zero which is bare earth - no vegetation, no water = no mosqitos? 

 # if(inf.pars$var.beta > 0) {
 #   wards$betaV <- inf.pars$betaV * (1 - wards[ward.names, "ndvi"]) / mean(1 - wards[ward.names, "ndvi"])
 # } else wards$betaV <- inf.pars$betaV

 # if(inf.pars$var.beta > 0) {
 #   wards$betaV <- inf.pars$betaV * (1 - wards[ward.names, "ndvi"]) / mean(1 - wards[ward.names, "ndvi"])
 # } else wards$betaV <- inf.pars$betaV

# how many vaccines are available? Inf means unlimited ----
n.vax.dose <- c(Inf, round(sum(wards[, paste0("pop.", sp)]) / 10))[2]

# how many replicates of each scenario to run? ----
nrep <- nrow(not.mkts.rvf.risk.tab) * 1

# initial node states ----
# these inital states are common to all simulations (outside the loop)
n <- length(ward.names)

u0.outer <- 
  data.frame(
    Sh = rep(0, n)
    ,Ih = rep(0, n)
    ,Rh = rep(0, n) 
    ,Pm = rep(20000, n) 
    ,Jm = rep(1000, n) 
    ,Sm = rep(1000, n)
    ,Em = rep(0, n) 
    ,Im = rep(0, n)
    ,Ic = rep(0, n) 
    )

rownames(u0.outer) <- ward.names

# add susceptible animals ----
u0.outer$Sh[!is.mkt] <- wards[ward.names[!is.mkt], paste0("pop.", sp)]

# add cumulative count of infecteds and vax column ----
u0.outer$Ic <- u0.outer$Ih

# inspect u0.outer
head(u0.outer)
dim(u0.outer)

