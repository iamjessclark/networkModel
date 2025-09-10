# Requirements ####
source("1. LoadItems2025.R")

if(!require(SimInf)){
  install.packages("SimInf")
  library(SimInf)
}


compartments <- c("Sh", "Ih", "Rh", "Pm", "Jm", "Sm","Em","Im", "Ic")
transitions <- c("@ -> muBirth*(Sh+Rh)*(1-(Sh+Ih+Rh)/Kh) -> Sh"
                 ,"Sh -> ((Sh+Ih+Rh) > 0) ? biteRate*p_vh*Im*Sh/(Sh+Ih+Rh) : 0 -> Ih+Ic"
                 ,"Ih ->  recovery*Ih -> Rh"
                 ,"Sh -> muH*Sh -> @"
                 ,"Ih -> muH*Ih -> @"
                 ,"Rh -> muH*Rh -> @"
                 ,"@ ->  laying*biteRate*(1-(Sm+Em+Im)/K) -> Pm"
                 ,"Pm -> development*Pm -> Jm"
                 ,"Jm -> firstBiteDelay*Jm -> Sm"
                 ,"Sm -> ((Sh+Ih+Rh) > 0) ? (p_hv*biteRate*Ih/(Sh+Ih+Rh))*Sm : 0 -> Em"
                 ,"Em -> extrinInc*Em -> Im"
                 ,"Pm -> pupsMort*Pm -> @"
                 ,"Sm -> mu*Sm -> @"
                 ,"Em -> mu*Em -> @"
                 ,"Im -> mu*Im -> @"
)


# Model Inputs ####

# n nodes
n <- 398

# time
years <- 1
maxTime <- 365*years

# inits
u0 <- u0.outer
u0$Sh <- u0$Sh+10
u0$Sm <- u0$Sh*100
events <- input.list[[1]]

unique(events$select)

events$select[which(events$select==6)] <- 3

# events matrix
n.compartments <- length(compartments)-1

Ev <- matrix(rep(0, n.compartments * 3),
             nrow = n.compartments,
             dimnames = list(compartments[1:n.compartments]))

Ev[substr(rownames(Ev), 1, 2) == "Sh", 1] <- 1 # just susceptible cattle
Ev[substr(rownames(Ev), 1, 2) == "Sh"|
     substr(rownames(Ev), 1, 2) == "Ih"|
     substr(rownames(Ev), 1, 2) == "Rh", 2] <- 1 # all cattle 
Ev[substr(rownames(Ev), 1, 2) == "Ih", 3] <- 1 # just infected hosts

Ev <- rbind(Ev, Ic = 0)

# # model
# VBDmodMovement <- mparse(
#                    transitions = transitions
#                  , compartments = compartments
#                  , gdata = c(
#                       recovery  = 0.2 
#                      , muH  = 0.000456621 
#                      , muBirth = 0.001369863
#                      , mu   = 0.03058104
#                      , p_hv = 0.25
#                      , p_vh = 0.04
#                      , laying = 20000
#                      , extrinInc = 0.1094668
#                      , K = 80000000
#                      , Kh = 100000
#                      , development = 0.1
#                      , biteRate = 0.3#0.68
#                      , firstBiteDelay = 0.14
#                      , pupsMort = 0.8
#                  )
#                  , u0 = u0
#                  , tspan = seq(1, maxTime)
#                  , events = events
#                  , E = Ev
# )
# 
# start <- Sys.time()
# resultMovement <- run(model = VBDmodMovement)
# (trajMovement <- trajectory(resultMovement))
# end <- Sys.time()
# (start-end)
#save.image("test.RData")


# trajMovement %>%
#   mutate(propIh = Ih/(Sh+Ih+Rh), 
#          propRh = Rh/(Sh+Ih+Rh), 
#          VtH = (Sm+Em+Im)/(Sh+Ih+Rh)) %>%
#   ggplot()+
#   geom_line(aes(x = time, y = VtH, group = node))

#### Sensitivity ####

muRange = c(0.03058104, 0.1428571) 
p_hvRange = c(0.001,1)
p_vhRange = c(0.001,1)
#layingRange = c(300,30000)
extrinIncRange = c(0.07142857, 0.1428571) 
#KRange = c(150000, 80000000)
biteRateRange = c(0.2, 0.68)
pupsMortRange = c(0.2,0.9)
firstBiteDelayRange = c(0.14, 0.5)

# put into vectors of min and max
min <- c(muRange[1]
         ,p_hvRange[1]
         ,p_vhRange[1]
         #,layingRange[1]
         ,extrinIncRange[1]
         #,KRange[1]
         ,biteRateRange[1]
         ,pupsMortRange[1]
         ,firstBiteDelayRange[1]
)

max <- c(muRange[2]
         ,p_hvRange[2]
         ,p_vhRange[2]
         #,layingRange[2]
         ,extrinIncRange[2]
         #,KRange[2]
         ,biteRateRange[2]
         ,pupsMortRange[2]
         ,firstBiteDelayRange[2]
)

# make a dataframe with ranges and names
params <- c(
  "mu"
  ,"p_hv"
  ,"p_vh"
  #,"laying"
  ,"extrinInc"
  #,"K"
  ,"biteRate"
  ,"pupsMort"
  ,"firstBiteDelay"
)

params <- cbind.data.frame(params,min,max)
nLHS <- 1000

# select random sets of parameter values within parameter value ranges given above
r <- randomLHS(nLHS,length(params[,1]))
# randomly sample from those random values
parmVals <- lapply(1:length(params[,1]),function(x){
  
  temp <- params[x,]
  randomSample <- runif(r[,x],min=temp$min,max=temp$max)
  
})
# list to df
parmVals <- do.call(cbind.data.frame,parmVals)
names(parmVals) <- params$params
parmVals$run <- 1:nLHS

(start <- Sys.time())
VBDmodMovementSA <- mclapply(parmVals$run, function(x){
  VBDmod <- mparse(transitions = transitions,
                   , compartments = compartments
                   , gdata = c(
                     recovery  = 0.2 
                     ,muH  = 0.000456621 
                     ,muBirth = 0.001369863
                     ,mu   = parmVals[x, "mu"]
                     ,p_hv = parmVals[x, "p_hv"] 
                     ,p_vh = parmVals[x, "p_vh"] 
                     ,laying = 300
                     ,extrinInc = parmVals[x, "extrinInc"]
                     ,K = 80000000
                     ,Kh = 500000
                     ,development = 0.1
                     ,biteRate = parmVals[x, "biteRate"]
                     ,firstBiteDelay = parmVals[x, "firstBiteDelay"]
                     ,pupsMort = parmVals[x, "pupsMort"]
                   )
                   , u0 = u0
                   , tspan = seq(1, maxTime)
                   , events = events
                   , E = Ev
  )
  
  
  result <- run(model = VBDmod)
  traj <- trajectory(result)
  
}  , mc.cores = detectCores()-2
) 
(end <- Sys.time())
start-end

save.image("simplenetworkSA.RDdata")

load("simplenetworkSA.RDdata")
testList <- list()
for(i in 1:length(VBDmodMovementSA)){
  VBDmodMovementSA[[i]]$run <- i
  
  testList[[i]] <- 
    VBDmodMovementSA[[i]] %>%
    mutate(propIh = Ih/(Sh+Ih+Rh), 
           propRh = Rh/(Sh+Ih+Rh), 
           VtH = (Sm+Em+Im)/(Sh+Ih+Rh)) %>%
    group_by(time) %>%
    summarise(meanIh = mean(propIh), 
              meanRh = mean(propRh), 
              meanVtH = mean(VtH)) %>%
    filter(time == max(time))
}

tempdf <- 
  do.call(rbind,testList)
tempdf$run <- 1:nrow(tempdf)
tempdf %>%
  filter(meanIh > 0) %>%
  ggplot()+
  geom_point(aes(x = meanIh, y = meanRh))


tempdf %>% 
  filter(meanRh > 0.1 & meanRh < 0.2)
  
parmVals[995,]
  
maxTime = 365*1 
# model
VBDmodMovementUpdatePars <- mparse(
                   transitions = transitions
                 , compartments = compartments
                 , gdata = c(
                      recovery  = 0.2
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
                 )
                 , u0 = u0
                 , tspan = seq(1, maxTime)
                 , events = events
                 , E = Ev
)

(start <- Sys.time())
resultMovementUpdatePars <- run(model = VBDmodMovementUpdatePars)
(trajMovementUpdatePars <- trajectory(resultMovementUpdatePars))
end <- Sys.time()
(start-end)
save.image("VBDmodMovementUpdatePars855.RData")

trajMovementUpdatePars %>%
  mutate(propIh = Ih/(Sh+Ih+Rh), 
         propRh = Rh/(Sh+Ih+Rh), 
         VtH = (Sm+Em+Im)/(Sh+Ih+Rh)) %>%
  ggplot() + 
  #geom_line(aes(x = time, y = propIh, group = node, colour = node)) +
  geom_line(aes(x = time, y = propRh, group = node, colour = node))
