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
      as.data.frame(as.matrix(
         u0$Sh*5
      ))
    
    colnames(ldata) <- "Kh"
      
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
                    #,ldata = ldata
                    ,u0 = mod.input$u0
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
    save.image("modout2025.RData")
    #save.image("networkcheck.RData")
#     
# Sh <- 
#   out %>%  
#     dplyr::select(time, node, Sh1:Sh64) %>%
#     rowwise() %>% 
#     mutate(wardNS = sum(c_across(Sh1:Sh64))) %>%
#     group_by(time) %>%
#     summarise(N = sum(wardNS)) %>%
#     mutate(state = "Sh")
# 
# Ih <- 
#   out %>%  
#   dplyr::select(time, node, Ih1:Ih64) %>%
#   rowwise() %>% 
#   mutate(wardNS = sum(c_across(Ih1:Ih64))) %>%
#   group_by(time) %>%
#   summarise(N = sum(wardNS)) %>%
#   mutate(state = "Ih")
# 
# 
# Rh <- 
#   out %>%  
#   dplyr::select(time, node, Rh1:Rh64) %>%
#   rowwise() %>% 
#   mutate(wardNS = sum(c_across(Rh1:Rh64))) %>%
#   group_by(time) %>%
#   summarise(N = sum(wardNS)) %>%
#   mutate(state = "Rh")
# 
# 
# dynamics <- 
#   bind_rows(Sh, Ih, Rh) %>%
#   pivot_wider(id_cols = time, values_from = N, names_from = state) %>%
#   mutate(propSh = (Sh/(Sh+Ih+Rh))*100, 
#          propIh = (Ih/(Sh+Ih+Rh))*100, 
#          propRh = (Rh/(Sh+Ih+Rh))*100
#          )
# 
# dynamics %>%
#   ggplot()+
#   geom_line(aes(x = time, y = propIh))
#   
# ####     u0[u0$Ic > 0, ]####
# #     # plot(result, node = which(u0$I > 0))
# #     
# #      event.tab <- table(model@events@event)
# #      event.report <-
# #        paste("There have been", event.tab[1], "deaths,", event.tab[2], "births, and",
# #              event.tab[3], "movements") # this isn't quite right anymore, becuase there is the D state
# #      print(event.report)
# #     
# #     # extract the total cumulative incidence across all wards (except markets) 
# #     # and the total number of animals
# #     use.wards <- match(not.mkts, rownames(u0))
# #     sir.comp <- mod.input$compartments[mod.input$compartments != "Ic"]
# #     i.comp <- sir.comp[substring(sir.comp, 1, 2) %in% c("Eh", "Ih")]
# #     cumInc.full <- sapply(use.wards, function(w) trajectory(result, index = w)$Ic)
# #     cumInc.full.prop <- 
# #       sapply(use.wards, function(w) {
# #         trj <- trajectory(result, index = w)
# #         out <- trj$Ic / rowSums(trj[, sir.comp])
# #         out[is.na(out) | is.infinite(out)] <- 0
# #         out
# #       })
# #   
# #     dynamics <- trajectory(model=result)
# #     enden.dynamics <- dynamics %>%
# #       filter(time > obvs.days) # just looking at the last 5 years. 
# #     
# #     Inc.full <- sapply(use.wards, function(w) rowSums(cbind(trajectory(result, index = w)[, i.comp])))
# #     Inc.full.prop <- 
# #       sapply(use.wards, function(w) {
# #         trj <- trajectory(result, index = w)
# #         out <- rowSums(cbind(trj[, i.comp])) / rowSums(trj[, sir.comp])
# #         out[is.na(out) | is.infinite(out)] <- 0
# #         out
# #       })
# #     
# #     Inc <- rowSums(Inc.full)
# #     cumInc <- rowSums(cumInc.full)
# #     N <- 
# #       rowSums(sapply(use.wards, function(w)
# #         rowSums(trajectory(result, index = w)[, sir.comp])))
# #     out <-
# #       cbind(
# #         intervention[k, ],
# #         day = 1:length(cumInc),
# #         cumInc = cumInc, 
# #         Inc = Inc, 
# #         N = N, 
# #         cumInc.ward = rowSums(cumInc.full.prop > rep.thresh),
# #         Inc.ward = rowSums(Inc.full.prop > rep.thresh),
# #         N.ward = n,
# #         cumInc.district = 
# #           sapply(1:n.days, 
# #                  function(d) length(unique(wards$district[cumInc.full[d, ] > 0]))),
# #         N.district = length(unique(wards$district)),
# #         row.names = NULL)
# #     attr(out, "cumInc.full.prop") <- cumInc.full.prop
# #     attr(out, "Inc.full.prop") <- Inc.full.prop
# #     attr(out, "n.ward.vax") <- n.ward.vax
# #     attr(out, "vax.used") <- vax.used
# #     attr(out, "event.report") <- event.report
# #     out
# #     
# #     #cumulative total
# #     out2 <- data.frame(matrix(nrow = n.ward[[3]], ncol= 2))
# #     out2$X1 <- sapply(seq(1:n.ward[[3]]), function(w) result@U[257*w,n.days])
# #     out2$X2 <- sapply(seq(1:n.ward[[3]]), function(w) ifelse(result@U[257*w,n.days]>0, 1, 0))
# #     colnames(out2) <- c("CumInc", "WardInf")
# #     out2
# #     # 
# #     out <- list(out, 
# #                 out2, 
# #                 enden.dynamics)
# #     
# #  #  }, mc.cores = detectCores()/2)
# # # print(end <- Sys.time())
# # #end-start
