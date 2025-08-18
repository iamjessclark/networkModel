#interventionbackup <- intervention
#intervention <- intervention[1,]
source("1. LoadItems2025.R")
# run model script
set.seed(123)

k = 1

print(start <- Sys.time())
  # run SEIR model----

#runMod <- mclapply(intervention$k, function(k) {
    # print(paste("Simulating disease transmission for scenario", k, "of", nrow(intervention)))
    # print(intervention[k, ])
 ### seed ward here
 # which ward(s) to put infected animals into
 seed.ward <-
   if(substr(intervention$seed.ward[k], 1, 4) == "rand") {
     sample(not.mkts, gsub("rand", "", intervention$seed.ward[k]))
   } else
     if(intervention$seed.ward[k] == "balanced") {
       not.mkts.rvf.risk.tab[intervention$rep[k] %% nrow(not.mkts.rvf.risk.tab) + 1, ]
     } else
       if(substr(intervention$seed.ward[k], 1, 8) == "rvf.risk") {
         sample(not.mkts, gsub("rvf.risk", "", intervention$seed.ward[k]),
                prob = wards[not.mkts, "rvf.risk"])
       } else
         intervention$seed.ward[k]

    # add infected animals to u0
    u0 <- u0.outer
    u0[seed.ward, "Ih"] <- u0[seed.ward, "Ic"] <- inf.pars$i0
    u0$Sh <- u0$Sh - u0$Ih
 
    # get vaccination plan for scenario k 
    vax.df <- vax.list[[intervention$vax[k]]]
    vax.per.ward <- 0
    if(all(vax.df$day == 0) & !(intervention$vax[k] %in% c("none"))) {
      #ward.names.ns <- ward.names[!ward.names %in% seed.ward]
      ward.names.ns <- ward.names
      vax.df <- vax.df[vax.df$ward %in% ward.names.ns, ]
      vax.df <- vax.df[1:min(intervention$n.vax[k], nrow(vax.df)), ]
      if(intervention$vax[k] == "rand") vax.df <- vax.df[sample(nrow(vax.df)), ]
      if(!(intervention$vax[k] %in% c("none", "all.lo", "all.hi"))) {
        vax.df <- vax.df[cumsum(u0[vax.df$ward, "Sh"] * vax.df$cov) < n.vax.dose, ]
      }
      vax.per.ward <- round(u0[vax.df$ward, "Sh"] * vax.df$cov)
      u0[vax.df$ward, "Rh"] <- u0[vax.df$ward, "Rh"] + vax.per.ward
      u0[vax.df$ward, "Sh"] <- u0[vax.df$ward, "Sh"] - vax.per.ward
    }
    n.ward.vax <- nrow(vax.df[vax.df$ward != "none", ])
    vax.used <- sum(vax.per.ward)
    if(is.na(vax.used)) vax.used <- 0
    print(intervention$vax[k])
    print(paste(n.ward.vax, "wards vaccinated.", vax.used, "doses used."))
    
    ## set up SEIR model ##
    
    # make sure starting values for states are integers 
    u0$Ic <- 0
    u0 <- as.data.frame(lapply(u0,as.integer))
    rownames(u0) <- ward.names
    #u0$Sh <- u0$Sh + 10
    
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
    
    # make host carrying capacity per node and set to the original starting values? 
     ldata <- 
       data.frame(
          muH = rep(intervention$muH[k], n.ward[3])
          ,muBirth = rep(intervention$muBirth[k], n.ward[3])
       )

    ldata[c(which(u0$Sh==0)),c(1,2)] <- 0
    
    # set up model 
    model <- mparse(transitions = mod.input$transitions,
                    compartments = as.character(mod.input$compartments),
                    gdata = 
                      c(recovery  = intervention$recovery[k] 
                        ,muH  = intervention$muH[k] 
                        ,muBirth = intervention$muBirth[k]
                        ,mu = intervention$mu[k]
                        ,p_hv = intervention$p_hv[k]
                        ,p_vh = intervention$p_vh[k]
                        ,laying = intervention$laying[k]
                        ,extrinInc = intervention$extrinInc[k]
                        ,development = intervention$development[k]
                        ,K = intervention$K[k]
                        ,Kh = 500000
                        ,biteRate = intervention$biteRate[k]
                        ,firstBiteDelay = intervention$firstBiteDelay[k]
                        ,pupsMort = intervention$pupsMort[k]
                        ,coupling = coupling
                        )
                    ,u0 = mod.input$u0
                    #,ldata = ldata
                    ,tspan = 1:n.days
                    ,events = input.list[[k]]
                    ,E = mod.input$Ev
                    ,N = N.mat
                    )

    # run the simulation
    print(Sys.time())
    result <- run(model)
    out <- trajectory(result)
    print(Sys.time())
    #save.image("modout2025.RData")
