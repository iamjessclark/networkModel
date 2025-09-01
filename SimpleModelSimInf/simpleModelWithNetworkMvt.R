# Requirements ####

if(!require(SimInf)){
  install.packages("SimInf")
  library(SimInf)
}

#,"Sh -> biteRate*p_vh*Im*Sh/(Sh+Ih+Rh) -> Ih+Ic"

compartments <- c("Sh", "Ih", "Rh", "Pm", "Jm", "Sm","Em","Im", "Ic")
transitions <- c("@ -> muBirth*(Sh+Rh)*(1-(Sh+Ih+Rh)/Kh) -> Sh"
                 ,"Sh -> (Sh > 0) ? biteRate*p_vh*Im*Sh/(Sh+Ih+Rh) : 0 -> Ih+Ic"
                 ,"Ih ->  recovery*Ih -> Rh"
                 ,"Sh -> muH*Sh -> @"
                 ,"Ih -> muH*Ih -> @"
                 ,"Rh -> muH*Rh -> @"
                 ,"@ ->  laying*biteRate*(1-(Sm+Em+Im)/K) -> Pm"
                 ,"Pm -> development*Pm -> Jm"
                 ,"Jm -> firstBiteDelay*Jm -> Sm"
                 ,"Sm -> (Ih > 0) ? (p_hv*biteRate*Ih/(Sh+Ih+Rh))*Sm : 0 -> Em"
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
     substr(rownames(Ev), 1, 2) == "Rh", 2] <- 1 # all cattle apart from the dead ones
Ev[substr(rownames(Ev), 1, 2) == "Ih", 3] <- 1 # just infected hosts

Ev <- rbind(Ev, Ic = 0)

# model
VBDmodMovement <- mparse(
                   transitions = transitions
                 , compartments = compartments
                 , gdata = c(
                      recovery  = 0.2 
                     , muH  = 0.000456621 
                     , muBirth = 0.001369863
                     , mu   = 0.03058104
                     , p_hv = 0.25
                     , p_vh = 0.04
                     , laying = 20000
                     , extrinInc = 0.1094668
                     , K = 80000000
                     , Kh = 100000
                     , development = 0.1
                     , biteRate = 0.3#0.68
                     , firstBiteDelay = 0.14
                     , pupsMort = 0.8
                 )
                 , u0 = u0
                 , tspan = seq(1, maxTime)
                 , events = events
                 , E = Ev
)

start <- Sys.time()
resultMovement <- run(model = VBDmodMovement)
(trajMovement <- trajectory(resultMovement))
end <- Sys.time()
(start-end)
#save.image("test.RData")


trajMovement %>%
  mutate(propIh = Ih/(Sh+Ih+Rh), 
         propRh = Rh/(Sh+Ih+Rh), 
         VtH = (Sm+Em+Im)/(Sh+Ih+Rh)) %>%
  ggplot()+
  geom_line(aes(x = time, y = VtH, group = node))

#### Sensitivity ####


