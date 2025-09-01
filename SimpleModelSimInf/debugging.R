is.error <- function(x) inherits(x, "try-error")
succeeded <- !vapply(stochSApups, is.error, logical(1))
good <- which(succeeded == T)
fails <- which(succeeded == F)

stochSApups2 <- stochSApups

for(i in 1:length(good)){
  stochSApups2[[good[i]]]$run <- NA
  stochSApups2[[good[i]]][,ncol(stochSApups2[[good[i]]])] <- good[i]
}

# make df
stochSAout2 <- do.call(rbind, stochSA)# Requirements ####

if(!require(SimInf)){
  install.packages("SimInf")
  library(SimInf)
}

# require(tidyverse)
# 
# require(DescTools)

#set.seed(123)

# Model Framework ####
compartments <- c("Sh", "Ih", "Rh", "Ic", "Nc","Dh", "Pm","Sm","Em","Im")
transitions <- c("@ -> muBirth*(Sh+Rh)*(1-(Sh+Ih+Rh)/Kh) -> Sh + Nc"
                 ,"Sh -> biterate*p_vh*Im*Sh/(Sh+Ih+Rh) -> Ih+Ic"
                 ,"Ih ->  recovery*Ih -> Rh"
                 ,"Sh -> muH*Sh -> Dh"
                 ,"Ih -> muH*Ih -> Dh"
                 ,"Rh -> muH*Rh -> Dh"
                 ,"@ ->  laying*biterate*(1-(Sm+Em+Im)/K) -> Pm"
                 ,"Pm -> er*Pm -> Sm"
                 ,"Sm -> (p_hv*biterate*Ih/(Sh+Ih+Rh))*Sm -> Em"
                 ,"Em -> epsilon*Em -> Im"
                 ,"Sm -> mu*Sm -> @"
                 ,"Em -> mu*Em -> @"
                 ,"Im -> mu*Im -> @"
)


# Model Inputs ####

# n nodes
n <- 1

# time
years <- 70
maxTime <- 365*years

# inits
u0 <- data.frame(Sh = rep(99, n), Ih = rep(1, n), Rh = rep(0, n), Dh = rep(0,n)
                 ,Pm=rep(1000,n),Sm = rep(1000,n),Em=rep(0,n),Im=rep(0,n))

u0$Ic <- u0$Ih
u0$Nc <- u0$Sh + u0$Ih + u0$Rh

ldata <- data.frame(biterate = 0.68 # unknown
)

VBDmod1 <- mparse(transitions = transitions
                  , compartments = compartments
                  , gdata = c(recovery  = 0.2 # host recovery
                              ,muH  = 0.000456621 #0.005479452 gives a steady infection but 100% sero+ when not differentiating birth and death rate
                              ,muBirth = 0.001369863
                              ,mu = 0.03058104 # https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-023-05792-3#:~:text=The%20average%20adult%20lifespan%20for,Culex%20species%20(Table%202)
                              ,p_hv = 0.25 # unknwown
                              ,p_vh = 0.04 # unknown
                              ,epsilon  = 0.1094668 # unknown extrinsic incubation - some sources say as little as 1 day in high temps
                              ,laying = 20000 # rafts of 300 eggs
                              ,er   = 0.1  # eggs take 10 days to develop into adults
                              ,K = 110000
                              ,Kh = 1200
                  )
                  
                  , ldata = ldata
                  , u0 = u0
                  , tspan = seq(1, maxTime)
)

result <- run(model = VBDmod1)
(traj <- trajectory(result))

# plot recovered and infected pop
stochSAout2 %>%
  # pivot_longer(cols = c(Sh:Im), names_to = "agent", values_to = "N") %>%
  # group_by(time, run, agent) %>%
  # summarise(N = sum(N)) %>%
  # pivot_wider(names_from = agent, values_from = N) %>%
  mutate(Shprop = Sh/(Sh+Ih+Rh),
         Ihprop = Ih/(Sh+Ih+Rh),
         Rhprop = Rh/(Sh+Ih+Rh)
  ) %>%
  select(time, run, Shprop:Rhprop) %>%
  pivot_longer(cols = Shprop:Rhprop, names_to = "agent", values_to = "Prop") %>%
  #group_by(agent, run, time) %>%
  #summarise(quantiles = quantile(Prop, c(0.025, 0.5, 0.975)), q = c("min", "mid", "upper")) %>%
  #pivot_wider(names_from = q, values_from = quantiles) %>%
  filter(agent == "Shprop") %>% 
  ggplot()+
  geom_line(aes(x = time/365, y = Prop, group = run)) 

endPoint <- 
  stochSAout2 %>%
  mutate(Shprop = Sh/(Sh+Ih+Rh),
         Ihprop = Ih/(Sh+Ih+Rh),
         Rhprop = Rh/(Sh+Ih+Rh)
  ) %>%
  filter(time == max(time)) 

endPointSH <- 
  endPoint %>%
  dplyr::select(run, Shprop)

endPointSH <- 
  parmVals %>%
  full_join(endPointSH, by = "run") 

#endPointSH <- endPointSH[-which(is.na(endPointSH$Shprop)),]

endPointSH <- 
  endPointSH %>%
  pivot_longer(-Shprop)

endPointSH %>%
  ggplot() +
  geom_point(aes(x = value, y = Shprop)) +
  facet_wrap(.~name, scales = "free")

# end point Rh
endPointRH <- 
  endPoint %>%
  dplyr::select(run, Rhprop)

endPointRH <- 
  parmVals %>%
  full_join(endPointRH, by = "run") 

#endPointRH <- endPointRH[-which(is.na(endPointRH$Rhprop)),]

endPointRH <- 
  endPointRH %>%
  pivot_longer(-Rhprop)

endPointRH %>%
  ggplot() +
  geom_point(aes(x = value, y = Rhprop)) +
  facet_wrap(.~name, scales = "free")

# end point Ih
endPointIH <- 
  endPoint %>%
  dplyr::select(run, Ihprop)

endPointIH <- 
  parmVals %>%
  full_join(endPointIH, by = "run") 

#endPointIH <- endPointIH[-which(is.na(endPointIH$Ihprop)),]

endPointIH <- 
  endPointIH %>%
  pivot_longer(-Ihprop)

endPointIH %>%
  ggplot() +
  geom_point(aes(x = value, y = Ihprop)) +
  facet_wrap(.~name, scales = "free")

# which pops have ~ 85% susceptible individuals? 
head(endPoint)
potentialsSh <- 
  endPoint %>%
  filter(Shprop > 0.5 & Shprop < 0.95) %>%
  dplyr::select(run, Shprop) %>%
  full_join(parmVals, by = "run") %>%
  filter(!is.na(Shprop))

goodRuns <- potentialsSh$run
# extract runs that have a sensible Sh at the end 
forPlotRunsSh <- subset(stochSAout2, run %in% goodRuns)
forPlotRunsSh %>%
  mutate(Shprop = Sh/(Sh+Ih+Rh),
         Ihprop = Ih/(Sh+Ih+Rh),
         Rhprop = Rh/(Sh+Ih+Rh)
  ) %>%
  ggplot() + 
  geom_line(aes( x = time, y = Shprop, group = run))

# what about infected dynamics
potentialsIh <- 
  endPoint %>%
  filter(Ihprop >0 & Ihprop < 0.05 ) %>%
  dplyr::select(run, Ihprop) %>%
  full_join(parmVals, by = "run") %>%
  filter(!is.na(Ihprop))

goodRunsIh <- potentialsIh$run
# extract runs that have a sensible Ih at the end 
forPlotRunsIh <- subset(stochSAout2, run %in% goodRunsIh)
forPlotRunsIh %>%
  mutate(Shprop = Sh/(Sh+Ih+Rh),
         Ihprop = Ih/(Sh+Ih+Rh),
         Rhprop = Rh/(Sh+Ih+Rh)
  ) %>%
  ggplot() + 
  geom_line(aes( x = time, y = Ihprop))

forPlotRunsIh %>%
  mutate(Shprop = Sh/(Sh+Ih+Rh),
         Ihprop = Ih/(Sh+Ih+Rh),
         Rhprop = Rh/(Sh+Ih+Rh)
  ) %>%
  ggplot() + 
  geom_line(aes( x = time, y = Shprop))

forPlotRunsIh %>%
  mutate(Shprop = Sh/(Sh+Ih+Rh),
         Ihprop = Ih/(Sh+Ih+Rh),
         Rhprop = Rh/(Sh+Ih+Rh)
  ) %>%
  ggplot() + 
  geom_line(aes( x = time, y = Rhprop))
