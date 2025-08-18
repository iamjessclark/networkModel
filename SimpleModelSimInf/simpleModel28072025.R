# Requirements ####

if(!require(SimInf)){
  install.packages("SimInf")
  library(SimInf)
}

compartments <- c("Sh", "Ih", "Rh", "Ic", "Pm", "Jm", "Sm","Em","Im")
transitions <- c("@ -> muBirth*(Sh+Rh)*(1-(Sh+Ih+Rh)/Kh) -> Sh"
                 ,"Sh -> biteRate*p_vh*Im*Sh/(Sh+Ih+Rh) -> Ih+Ic"
                 ,"Ih ->  recovery*Ih -> Rh"
                 ,"Sh -> muH*Sh -> @"
                 ,"Ih -> muH*Ih -> @"
                 ,"Rh -> muH*Rh -> @"
                 ,"@ ->  laying*biteRate*(1-(Sm+Em+Im)/K) -> Pm"
                 ,"Pm -> development*Pm -> Jm"
                 ,"Jm -> firstBiteDelay*Jm -> Sm"
                 ,"Sm -> (p_hv*biteRate*Ih/(Sh+Ih+Rh))*Sm -> Em"
                 ,"Em -> extrinInc*Em -> Im"
                 ,"Pm -> pupsMort*Pm -> @"
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
# inits
u0 <- data.frame(Sh = rep(99, n), Ih = rep(1, n), Rh = rep(0, n), Dh = rep(0,n)
                 ,Pm=rep(1000,n), Jm = rep(0, n), Sm = rep(1000,n), Em=rep(0,n), Im=rep(0,n))

u0$Ic <- u0$Ih

VBDmod <- mparse(transitions = transitions
                   , compartments = compartments
                   , gdata = c(
                     recovery  = 0.2 
                     ,muH  = 0.000456621 
                     ,muBirth = 0.001369863
                     ,mu   = 0.03058104
                     ,p_hv = 0.25
                     ,p_vh = 0.04
                     ,laying = 20000
                     ,extrinInc = 0.1094668
                     ,K = 110000
                     ,Kh = 1200
                     ,development = 0.1
                     ,biteRate = 0.68
                     ,firstBiteDelay = 0.14
                     ,pupsMort = 0.2
                   )
                   , u0 = u0
                   , tspan = seq(1, maxTime)
                 )

start <- Sys.time()
result <- run(model = VBDmod)
(traj <- trajectory(result))
end <- Sys.time()
(start-end)
save.image("test.RData")