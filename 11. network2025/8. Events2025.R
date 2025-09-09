# simulating movement events with parallel computing ####
# this whole bit simulates movement between differennt subnodes given some measure
# betweeness, in/out centrality etc
# loop over different policy options
# making the events (movements) under each option

start.time <- Sys.time()
input.list <- mclapply(intervention$k, function(k) {
  print(paste("Simulating events for scenario", k, "of", nrow(intervention)))
  print(intervention[k, ])
  
  # this code doesn't work if intervention$seed.ward is e.g. "balanced" or "rand"
  ward.names.ns <-
    if(exclude.seed) ward.names[!ward.names %in% intervention$seed.ward[k]] else ward.names
  
  # select movement ban
  mban.all <- mban.list[[intervention$mban[k]]]
  if(intervention$mban[k] == "rand") {
    mban.all$ward <- sample(mban.all$ward)
  }
  
  # load parameters for simulating monthly movements
  par.list <-
    lapply(all.months, function(m) {
      print(m)
      
      # predicted probability of movement
      mu.z <- mu.z.list[[m]] * movt.multiplier
      
      # predicted number moved given movement
      mu.c <- mu.c.list[[m]]
      
      # dispersion (theta) parameter for the ZTNB distribution
      theta <- theta.list[[m]]
      
      # problem:
      # because the pink slips record only onward journeys from markets (n = 111), 
      # wards with no market (n = 287) appear to produce no animals. 
      # a simple solution is to assume that they send all their animals to
      # the nearest market. this could be improved by taking account of Gemma's
      # survey data, and also the trade-off between distance and size.
      # go with the simple solution for now.
      # identify wards that sent out no animals
      zero.out <- data.frame(ward = rownames(mu.z)[rowSums(mu.z) == 0], stringsAsFactors = FALSE)
      rowSums(mu.z[zero.out$ward, ]) # zero probability of outflow
      rowSums(mu.c[zero.out$ward, ]) # but nonzero rate given any outflow
      # for those wards that send out no animals 
      # assume that this is due to missing data and assume that all
      # animals from this ward went to the nearest active primary wards (active = outflow > 0) 
      zero.out$nearest.active.primary <- 
        sapply(zero.out$ward, function(w) {
          nap.index <- 
            which.nearest(
              unlist(wards[w, c("lat", "long")]),
              wards[!wards$fullname %in% zero.out$ward & wards$market == "Primary", 
                    c("lat", "long")])
          names(nap.index)
        })
      # how many primary markets have outflow now, but no inflow?
      rownames(mu.z)[rowSums(mu.z) > 0 & !rownames(mu.z) %in% zero.out$nearest.active.primary]
      wards[rownames(mu.z)[rowSums(mu.z) > 0 & !rownames(mu.z) %in% zero.out$nearest.active.primary], "prodsys"]
      # for now, just assume that all of these wards are supplied internally.
      # now distribute the outflow from each "nearest active primary" over 
      # the upstream wards with no recorded outflow
      for(nap in unique(zero.out$nearest.active.primary)) {
        mu.z[zero.out$ward[zero.out$nearest.active.primary == nap], nap] <- 
          1/sum(zero.out$nearest.active.primary == nap)
      }; rm(nap)
      sort(rowSums(mu.z))
      # diagonals should all be zero
      all(diag(mu.z) == 0)
      
      # apply movement ban via the movement probability matrix mu.z
      if(length(mban.all$ward) > 0 & intervention$mban[k] != "none" & m %in% mban.all$mban.months) {
        keep.wards <- mban.all$ward[mban.all$ward %in% ward.names.ns]
        keep.wards <- keep.wards[1:(min(length(ward.names.ns), max(1, intervention$n.mban[k])))]
        mu.z[keep.wards, ] <- 0
        mu.z[, keep.wards] <- 0
      }
      
      # multiply matrices to give expected flow
      mu <- mu.z * mu.c
      # check all wards now have outputs
      print(table(rowSums(mu) == 0))
      # output
      list(mu.z = mu.z, mu.c = mu.c, mu = mu, theta = theta)
    })
  names(par.list) <- all.months
  
  # remove secondary market wards from the contiguous wards network
  # because these wards are assumed to have no standing population of animals
  contig.wards <- contig.wards.all[contig.wards.all$from %in% not.mkts & contig.wards.all$to %in% not.mkts, ]
  
  # simulate balanced moves
  sim.moves.list <- 
    lapply(par.list, 
           function(mpar) 
             sim.hurdle.move(mu.z = mpar$mu.z, mu.c = mpar$mu.c, alpha = mpar$theta, is.mkt = is.mkt))
  
  # create events
  # start with event type "extTrans": moves between nodes
  
  # turn the move matrix list into a table of pairwise moves
  contig.move.n <- intervention$contig.move.n[k]
  evts.list <-
    lapply(names(sim.moves.list), 
           function(mth) {
             print(mth)
             #             date1 <- as.Date(paste("2015", mth, "01", sep = "-"))
             date1 <- all.months.day1[mth]
             move <- sim.moves.list[[mth]]
             
             # make a data frame of movements, one row per movement
             move.pair <- 
               expand.grid(fr.ward = rownames(move), 
                           to.ward = colnames(move), 
                           time = NA, n.moved = NA,
                           stringsAsFactors = FALSE)
             rownames(move.pair) <- paste(move.pair$fr.ward, move.pair$to.ward, sep = ".")
             
             
             # remove rows where no animals were moved
             # plus contiguous ward pairs, if any neighbour-neighbour movement 
             dim(move.pair)
             keep.pairs <- 
               apply(which(move > 0, arr.ind = TRUE), 1, 
                     function(i) paste(rownames(move)[i[1]], colnames(move)[i[2]], sep = "."))
             if(contig.move.n > 0) keep.pairs <- unique(c(keep.pairs, rownames(contig.wards)))
             move.pair <- move.pair[keep.pairs, ]
             dim(move.pair)
             
             # separate into the four different types of movement
             move.order <- c("wm", "mm", "mw", "ww")
             move.pair$move.type <- 
               factor(tolower(paste0(c("w", "m")[(move.pair$fr.ward %in% mkts) + 1], 
                                     c("w", "m")[(move.pair$to.ward %in% mkts) + 1])), 
                      move.order)
             
             for(m in move.order) {
               move.pair$time[move.pair$move.type == m] <- match(m, move.order) + date1
               move.pair$n.moved[move.pair$move.type == m] <- 
                 diag(
                   as.matrix(
                     move[move.pair$fr.ward[move.pair$move.type == m], 
                          move.pair$to.ward[move.pair$move.type == m]]))
               move.pair
             }; rm(m)
             move.pair <- move.pair[order(move.pair$time), ]
             
             # check that the same total number of animals is being moved
             sum(move.pair$n.moved) == sum(move)
             
             # simulate movements between contiguous wards
             # by moving an extra contig.move.n animals between all contiguous wards each month
             if(contig.move.n > 0) {
               move.pair[rownames(contig.wards), "n.moved"] <- 
                 move.pair[rownames(contig.wards), "n.moved"] + 
                 rpois(nrow(contig.wards), lambda = contig.move.n)
             }
             
             print(paste(round(nrow(contig.wards) * contig.move.n), sp, 
                         "moved between neighbours in month", mth))
             
             # output data frame of movement events
             if(nrow(move.pair) > 0.5) {
               move.df <- 
                 data.frame(
                   event = "extTrans", # options: exit (death), enter (birth), intTrans (e.g. S -> I), extTrans (movement)
                   # Events scheduled at same time processed inorder: exit, enter, internal trans, external trans
                   time = move.pair$time,
                   node = match(move.pair$fr.ward, ward.names),
                   dest = match(move.pair$to.ward, ward.names),
                   n = move.pair$n.moved,
                   proportion = 0,
                   select = 2, # chooses the column (compartment) of E the event operates on
                   shift = 0,
                   stringsAsFactors = FALSE)
             } else {
               move.df <- 
                 data.frame(event = character(0), time = integer(0), node = integer(0),
                            dest = integer(0), n = integer(0), proportion = integer(0),
                            select = integer(0), shift = integer(0)) 
             }
             
             
             # what proportion of movements are of 1 animal 
             print(paste0(round(100 * mean(move.df$n == 1)), "% of market movements have batch size = 1 in month ", mth))
             # measure changes in ward populations for each node
             active.nodes <- sort(unique(c(move.df$node, move.df$dest)))
             active.nodes <- active.nodes[!active.nodes %in% which(is.mkt)]
             inout.df.list <- 
               lapply(active.nodes, function(j) {
                 # n out
                 n.out <- sum(move.df$n[move.df$node == j & move.df$event == "extTrans"])
                 # n in
                 n.in <- sum(move.df$n[move.df$dest == j & move.df$event == "extTrans"])
                 imbalance <- n.in - n.out
                 # create births or deaths to correct imbalance in in/out-flow
                 if(n.out != n.in) {
                   demog.df <- 
                     data.frame(
                       event = ifelse(imbalance > 0, "exit", "enter"), # options: exit (death), enter (birth), intTrans (e.g. S -> I), extTrans (movement)
                       time = abs(max(sign(imbalance) * move.pair$time)) + sign(imbalance),
                       node = j,
                       dest = 0,
                       n = abs(imbalance),
                       proportion = 0,
                       select = ifelse(imbalance > 0, 2, 1), # chooses the column (compartment) of E the event operates on
                       shift = 0,
                       stringsAsFactors = FALSE)
                   return(demog.df)
                 } else return(NULL)
               })
             output.df <- rbind(move.df, do.call("rbind", inout.df.list))
             # check balance between inflow, outflow, births and deaths 
             imbalances <- 
               sapply(1:n, function(i) {
                 sum(output.df$n[output.df$event == "extTrans" & output.df$dest == i]) -
                   sum(output.df$n[output.df$event == "extTrans" & output.df$node == i]) +
                   sum(output.df$n[output.df$event == "enter" & output.df$node == i]) - 
                   sum(output.df$n[output.df$event == "exit" & output.df$node == i])
               })
             stopifnot(all(imbalances == 0))
             if(nrow(output.df) > 0.5) output.df$month <- mth
             output.df
           })
  sapply(evts.list, dim)
  
  # rbind the movements into a data frame
  evts <- do.call("rbind", evts.list)
  table(evts$time)
  
  # select vaccination schedule
  vax.df <- vax.list[[intervention$vax[k]]]
  
  # add vaccination event
  if(!all(vax.df$day == 0) & any(vax.df$day <= max(evts$time)) & any(vax.df$cov > 0)) {
    print("Adding vax events")
    if(intervention$vax[k] == "rand") vax.df <- vax.df[sample(nrow(vax.df)), ]
    vax.df <- vax.df[vax.df$ward %in% ward.names.ns, ]
    vax.df <- vax.df[1:min(intervention$n.vax[k], nrow(vax.df)), ]
    vaccination <- 
      data.frame(event = "intTrans", 
                 time = vax.df$day, 
                 node = match(vax.df$ward, ward.names), 
                 dest = 0, 
                 n = 0, 
                 proportion = vax.df$cov,
                 select = 1, 
                 shift = 1,
                 month = NA,  # don't need month for vax events
                 stringsAsFactors = FALSE)
    evts <- rbind(evts, vaccination)
  }
  
  
  evts <- evts[order(evts$time), ]
  evts
  
  
  # add the seasonal pupae emergence 
  # lets roughly say based on chirps data that rainy season is 
  # march - may and october - december 
  # then lets work out those days and make that our entrance days for pups
  
  # timeHatch <-
  #   which(month(seq(start.date, end.date, by = "days"))==3|
  #           month(seq(start.date, end.date, by = "days"))==4|
  #           month(seq(start.date, end.date, by = "days"))==5|
  #           month(seq(start.date, end.date, by = "days"))==10|
  #           month(seq(start.date, end.date, by = "days"))==11|
  #           month(seq(start.date, end.date, by = "days"))==12)
  # 
  # # then create the events E matrix for the events table to work off
  #  pupemergence <-
  #    data.frame(event = "enter",
  #               time = rep(timeHatch, length(ward.names)),
  #               node = 1:length(ward.names),
  #               dest = 0,
  #               n = wards$ndvi*10000, # this is making the number of mosquitoes emerging each season dependent on NDVI
  #               proportion = 0,
  #               select = 3,
  #               shift = 0,
  #               month = NA,
  #               stringsAsFactors = FALSE) # don't need month for hatching events
  # 
  # evts <- rbind(evts, pupemergence)
  # 

  infAnimal <- 
    data.frame(event = "enter", 
               time = 1,
               node = which(wards$rvf.risk == max(wards$rvf.risk))[1], 
               dest = 0, 
               n = 1, 
               proportion = 0, 
               select = 3, 
               shift = 0,
               month = NA, 
               stringsAsFactors = FALSE)
  
  evts <- rbind(evts, infAnimal)
  
  # inspect events table
  head(evts)
  dim(evts)
  table(evts$event) # 0 = exit; 1 = enter; 2 = internal transfer; 3 = external transfer
  
  
  # output events table
  evts 
  
}, mc.cores = detectCores() / 2, mc.preschedule = TRUE)

intervention$seed.ward <- "rvf.risk01"

plans <- 
  if(all(intervention$mban == "none")) vax.plans else
    if(all(intervention$vax == "none")) mban.plans else
      unique(intervention$plan)


# check all input events objects are data frames 
stopifnot(sapply(input.list, class) == "data.frame")

