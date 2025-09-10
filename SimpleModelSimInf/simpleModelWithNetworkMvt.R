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

load("simplenetworkSA.RDdata")

# model
# pars from row 855

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
