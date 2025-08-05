# Requirements ####

source("1. LoadItems2025.R")

set.seed(123)


interventionsSA <- intervention[1,]

muRange = c(0.03058104, 0.1428571) 
p_hvRange = c(0.001,1)
p_vhRange = c(0.001,1)
layingRange = c(300,30000)
extrinIncRange = c(0.07142857, 0.1428571) 
KRange = c(50000, 500000)
biteRateRange = c(0.02, 0.68)
pupsMortRange = c(0.001,1)
firstBiteDelayRange = c(0.14, 0.5)
#developmentRange = c(0.1, round((1/1.5),3))

# put into vectors of min and max
min <- c(muRange[1]
         ,p_hvRange[1]
         ,p_vhRange[1]
         ,layingRange[1]
         ,extrinIncRange[1]
         ,KRange[1]
         ,biteRateRange[1]
         ,pupsMortRange[1]
         ,firstBiteDelayRange[1]
#         ,developmentRange[1]
         )

max <- c(muRange[2]
         ,p_hvRange[2]
         ,p_vhRange[2]
         ,layingRange[2]
         ,extrinIncRange[2]
         ,KRange[2]
         ,biteRateRange[2]
         ,pupsMortRange[2]
         ,firstBiteDelayRange[2]
#         ,developmentRange[2]
         )

# make a dataframe with ranges and names
params <- c(
  "mu"
  ,"p_hv"
  ,"p_vh"
  ,"laying"
  ,"extrinInc"
  ,"K"
  ,"biteRate"
  ,"pupsMort"
  ,"firstBiteDelay"
#  ,"development"
)

params <- cbind.data.frame(params,min,max)
nLHS <- 10

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

# run SEIR model----
### seed ward here
# which ward(s) to put infected animals into
seed.ward <- 
  if(substr(interventionsSA$seed.ward, 1, 4) == "rand") {
    sample(not.mkts, gsub("rand", "", interventionsSA$seed.ward))
  } else
    if(interventionsSA$seed.ward == "balanced") {
      not.mkts.rvf.risk.tab[interventionsSA$rep %% nrow(not.mkts.rvf.risk.tab) + 1, ]
    } else
      if(substr(interventionsSA$seed.ward, 1, 8) == "rvf.risk") {
        sample(not.mkts, gsub("rvf.risk", "", interventionsSA$seed.ward), 
               prob = wards[not.mkts, "rvf.risk"]) 
      } else
        interventionsSA$seed.ward


# add infected animals to u0
u0 <- u0.outer
u0[seed.ward, "Ih"] <- u0[seed.ward, "Ic"] <- inf.pars$i0
u0$Sh <- u0$Sh - u0$Ih

# get vaccination plan for scenario k 
vax.df <- vax.list[[interventionsSA$vax]]
vax.per.ward <- 0
if(all(vax.df$day == 0) & !(interventionsSA$vax %in% c("none"))) {
  #ward.names.ns <- ward.names[!ward.names %in% seed.ward]
  ward.names.ns <- ward.names
  vax.df <- vax.df[vax.df$ward %in% ward.names.ns, ]
  vax.df <- vax.df[1:min(interventionsSA$n.vax, nrow(vax.df)), ]
  if(interventionsSA$vax == "rand") vax.df <- vax.df[sample(nrow(vax.df)), ]
  if(!(interventionsSA$vax %in% c("none", "all.lo", "all.hi"))) {
    vax.df <- vax.df[cumsum(u0[vax.df$ward, "Sh"] * vax.df$cov) < n.vax.dose, ]
  }
  vax.per.ward <- round(u0[vax.df$ward, "Sh"] * vax.df$cov)
  u0[vax.df$ward, "Rh"] <- u0[vax.df$ward, "Rh"] + vax.per.ward
  u0[vax.df$ward, "Sh"] <- u0[vax.df$ward, "Sh"] - vax.per.ward
}
n.ward.vax <- nrow(vax.df[vax.df$ward != "none", ])
vax.used <- sum(vax.per.ward)
if(is.na(vax.used)) vax.used <- 0
print(interventionsSA$vax)
print(paste(n.ward.vax, "wards vaccinated.", vax.used, "doses used."))

## set up SEIR model ##

# make sure starting values for states are integers 
u0$Ic <- 0
u0 <- as.data.frame(lapply(u0,as.integer))
rownames(u0) <- ward.names

# split nodes into subnodes 
mod.input <- split.nodes.sir(n.subnodes = n.subnodes, u0 = u0)
#mod.input <- split.nodes.sir(n.subnodes = 1, u0 = u0, seir = TRUE)

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



(start <- Sys.time())
#### Run model ####
x = 1
#networkSA2025 <- pbmcapply::pbmclapply(parmVals$run, function(x){
#networkSA2025 <- mclapply(parmVals$run, function(x...){
  VBDmod <- mparse(transitions = mod.input$transitions,
                   , compartments = as.character(mod.input$compartments)
                   , gdata = c(
                     recovery  = 0.2 
                     ,muH  = 0.000456621 
                     ,muBirth = 0.001369863
                     ,mu   = parmVals[x, "mu"]
                     ,p_hv = parmVals[x, "p_hv"] 
                     ,p_vh = parmVals[x, "p_vh"] 
                     ,laying = parmVals[x, "laying"] 
                     ,extrinInc = parmVals[x, "extrinInc"]
                     ,K = parmVals[x, "K"]
                     ,Kh = 500000
                     ,development = 0.1
                     ,biteRate = parmVals[x, "biteRate"]
                     ,firstBiteDelay = parmVals[x, "firstBiteDelay"]
                     ,pupsMort = parmVals[x, "pupsMort"]
                     ,coupling = coupling
                   )
                   ,u0 = mod.input$u0
                   ,tspan = c(1, n.days*5:n.days)
                   ,events = as.data.frame(input.list[1])
                   ,E = mod.input$Ev
                   ,N = N.mat
                   )

  
  result <- run(model = VBDmod)
  traj <- trajectory(result)
  save.image("SAnetwork2025.RData")
  #out <- 
  #  out %>%
  #  filter(time > 365*5)
  #out
  
#  }  , mc.cores = detectCores()/2
#) 
#(end <- Sys.time())
#start-end



#### # extract the total cumulative incidence across all wards (except markets) ####
# # and the total number of animals
# use.wards <- match(not.mkts, rownames(u0))
# sir.comp <- mod.input$compartments[mod.input$compartments != "Ic"]
# i.comp <- sir.comp[substring(sir.comp, 1, 2) %in% c("Eh", "Ih")]
# cumInc.full <- sapply(use.wards, function(w) trajectory(result, index = w)$Ic)
# cumInc.full.prop <- 
#   sapply(use.wards, function(w) {
#     trj <- trajectory(result, index = w)
#     out <- trj$Ic / rowSums(trj[, sir.comp])
#     out[is.na(out) | is.infinite(out)] <- 0
#     out
#   })

#dynamics <- trajectory(model=result)
# enden.dynamics <- dynamics %>%
#   filter(time > obvs.days) # just looking at the last 5 years. 
# 
# ShWardNs <- dynamics %>%
#   dplyr::select(time, Sm1:Sm64) %>%
#   rowwise() %>% 
#   mutate(wardNS = sum(c(Sm1:Sm64))) %>%
#   dplyr::select(time, wardNS)
# 
# ShWardNs$ward <- rep(ward.names, n.days-1)
# yearT <- which(day(seq(start.date, end.date, by = "days"))==1 & month(seq(start.date, end.date, by = "days"))==1) 
# ShWardNs <- ShWardNs %>%
#   mutate(year = ifelse(time >= 2 & time <= yearT[2]-1, 1,
#                        ifelse(time >= yearT[2] & time <= yearT[3]-1, 2,
#                               ifelse(time >= yearT[3] & time <= yearT[4]-1, 3,
#                                      ifelse(time >= yearT[4] & time <= yearT[5]-1,  4,
#                                             ifelse(time >= yearT[5] & time <= yearT[6]-1, 5,
#                                                    ifelse(time >= yearT[6] & time <= yearT[7]-1, 6,
#                                                           ifelse(time >= yearT[7] & time <= yearT[8]-1, 7,
#                                                                  ifelse(time >= yearT[8] & time <= yearT[9]-1, 8,
#                                                                         ifelse(time >= yearT[9] & time <= yearT[10]-1, 9, 10))))))))))
# 
# ShWardNs <- ShWardNs %>%
#   group_by(ward, year) %>%
#   summarise(NSyearWard = mean(wardNS))
# 
# endemicases <- yearlycases %>%
#   filter(casevec!="NA") # casevec is a numeric vector so doesn't recognise is.na()
# endemicases$year <- as.double(endemicases$year)
# 
# # our FOIs - ie percentage of susceptible individuals getting infected each year 
# # is within the credible intervals of Will's analysis. 
# 
# FOIest <- ShWardNs %>%
#   filter(year >= 5) %>%
#   left_join(endemicases, by = join_by(year, ward)) %>%
#   mutate(foi = (casevec/NSyearWard)*100)
# 
# meanFOI <- FOIest %>%
#   group_by(ward) %>%
#   summarise(meanFOI = mean(foi))
# 
# 
# 
# }, mc.cores = detectCores()/2)
# print(end <- Sys.time())
# 
# end-start
