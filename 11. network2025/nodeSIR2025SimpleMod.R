# SIR model within each node ----

# function to split nodes into subnodes to create within-node heterogeneity ----
# the subnodes are arranged on a square grid, so sqrt(n.subnodes) must be an integer
split.nodes.sir<- 
  function(n.subnodes, u0, adj = NULL) { # removing seir= arguemnt
    
    # check that n.subnodes subnodes can be arranged in a square grid
    stopifnot(sqrt(n.subnodes) == round(sqrt(n.subnodes)))
    
    # define compartments
    compartments <- c("Sh", "Ih", "Rh", "Pm", "Jm", "Sm","Em","Im", "Ic")
    n.compartments <- length(compartments)-1
    hostcompartments <- c("Sh", "Ih", "Rh")
    vectorcompartments <- c("Pm", "Jm", "Sm","Em","Im") 
    
    #here I have removed the if statement saying to remove "E" if seir = false
    # and have added in eggs and pups
    
    # split the compartments into subnodes - these subnodes are labeled according to how many subnodes there are
    # so if you have a 3x3 grid in each node, then there are 9 epi models - one for each subnode, so there are 9 x s 
    # 9 x I, 9 x R etc. 
    
    sub.comp <- 
      paste0(rep(compartments[1:8], each = n.subnodes), rep(1:n.subnodes, length(compartments[1:8]))) 
    
    # how many subcompartments are there
    n.sub.comp <- length(sub.comp)
    # create adjacency matrix, if none is supplied
    # the adjacency matrix tells me which subnodes in teh matrix are next to each other
    if(is.null(adj)) {
      d <- sqrt(n.subnodes)
      # so many an empty adjacency matrix htat has ncols = n.subnodes and nrows = n.subnodes
      # then in each cell of the adjacency matrix, you're saying eg. subnode 1 is next to subnode 2 and subnode 4 so row 1, columns 2 and 4 get 1's 
      # and the subnodes not next to each other have 0s eg. row 1, col 1, 3, 5, 6, 7, 8, and 9 have 0's. 
      adj <- 
        structure(rep(0, n.subnodes^2), dim = c(n.subnodes, n.subnodes), 
                  dimnames = list(1:n.subnodes, 1:n.subnodes))
      for(i in 1:nrow(adj)) {
        for(j in 1:ncol(adj)) {
          # detect if subnodes i and j are adjacent on a d X d square matrix. if so, indicate on adjacency matrix with 1,
          # otherwise leave at 0. the spatial arrangement of subnodes is matrix(1:n.subnodes, nrow = d)
          # and adjacency is defined as being one sideways move away, so diagonal neighbours are not adjacent
          if(((i - j) == 1 && j%%d != 0) | ((i - j) == -1 && i%%d != 0) | abs(i - j) == d) adj[i, j] <- 1
        }; rm(j)
      }; rm(i)
    }
    
    stopifnot(all(dim(adj) == n.subnodes))
    # create the list of transitions between compartments, including contamination from coupled (adjacent)
    transitions.tab <-
      # using sn as an index
      sapply(1:n.subnodes, function(sn) {
        # in each row, which columns have 1s in them
        # store these as a series of numbers that get appended to the appropriate state vars below
        nbrs <- rownames(adj)[adj[sn, ] == 1]
        
        # this is saying that if there is contamination of subnode x from the surrounding
        # subnodes, then add in the contamination from the infected individuals in that/ those 
        # adjacent subnodes
        env.contamSh <- 
          if(length(nbrs) > 0) paste0(" + ((", paste0("Sh", nbrs, collapse = "+"), ")*coupling)") else NULL
        
        env.contamIh <- 
          if(length(nbrs) > 0) paste0(" + ((", paste0("Ih", nbrs, collapse = "+"), ")*coupling)") else NULL
        
        env.contamRh <- 
          if(length(nbrs) > 0) paste0(" + ((", paste0("Rh", nbrs, collapse = "+"), ")*coupling)") else NULL
        
        Nh <- paste0("(", paste0(hostcompartments, sn, collapse = "+"), ")")
        Nm <- paste0("(", paste0(vectorcompartments, sn, collapse = "+"), ")")
        
        # Hosts ####
        #### Infection processes ####
        sh.trans <-
          paste0("Sh", sn, " -> ((Sh", sn, "+Ih", sn, "+Rh", sn, ") > 0) ? biteRate*p_vh*Im", sn, "*Sh", sn , "/(", Nh , ") : 0 -> Ih", sn, " + Ic")
        
        ih.trans <-
          paste0("Ih", sn, " -> recovery*Ih", sn, " -> Rh", sn)
        
        #### Deaths ####
        sh.trans.D <-
          paste0("Sh", sn, " -> muH*Sh", sn, " -> @") 
        
        ih.trans.D <-
          paste0("Ih", sn, " -> muH*Ih", sn, " -> @")
        
        rh.trans.D <-
          paste0("Rh", sn, " -> muH*Rh", sn, " -> @")
        
        #### Host Births ####
        hbirths <-
          paste0("@ -> muBirth*(Sh", sn,"+Rh",sn,")*(1-", Nh, "/Kh) -> Sh", sn)
        
        # Mosquitos ####
        
        pups.trans <- 
          paste0("Pm", sn, " -> development*Pm", sn, "-> Jm", sn)  
        
        junior.trans <- 
          paste0("Jm", sn, "-> firstBiteDelay*Jm", sn, "-> Sm", sn)
        
        #### Infection processes ####
        sm.trans <- 
          paste0("Sm", sn, " -> (", Nh, "> 0) ? (p_hv*biteRate*Ih", sn, env.contamIh, ")/", Nh,"*Sm", sn, ": 0 -> Em", sn)

        em.trans <- 
          paste0("Em", sn, " -> extrinInc*Em", sn, "-> Im", sn)

        #### Deaths ####
        pups.trans.D <-
          paste0("Pm", sn, " -> pupsMort*Pm", sn, " -> @")
        
        jm.trans.D <-
          paste0("Jm", sn, " -> mu*Jm", sn, " -> @")
        
        sm.trans.D <-
          paste0("Sm", sn, " -> mu*Sm", sn, " -> @")
        
        em.trans.D <-
          paste0("Em", sn, " -> mu*Em", sn, " -> @")

        im.trans.D <-
          paste0("Im", sn, " -> mu*Im", sn, " -> @")
        
        #### births ####
        mbirths <-
          paste0("@ -> laying*biteRate*(1-",Nm,"/K) -> Pm", sn)
    
        c(sh.trans
          , ih.trans
          , hbirths
          , sh.trans.D
          , ih.trans.D
          , rh.trans.D
          , mbirths
          , pups.trans
          , junior.trans
          , sm.trans
          , em.trans
          , pups.trans.D
          , jm.trans.D
          , sm.trans.D
          , em.trans.D
          , im.trans.D
        )
      })
    
    transitions <- c(t(transitions.tab))
    # split the initial numbers in each compartment for each node across the subnodes.
    # numbers from each node are randomly (mulitnomially) allocated to subnodes
    u0.list.host <- 
      #lapply(hostcompartments,
      lapply(hostcompartments, 
             function(cmp) {
               out <- t(sapply(1:nrow(u0), function(i) rmultinom(1, u0[i, cmp], rep(1, n.subnodes))))
               if(n.subnodes == 1) out <- t(out)
               colnames(out) <- paste0(cmp, 1:n.subnodes)
               out
             })
    
    
    u0.list.mosq <- 
      lapply(vectorcompartments, 
             function(cmp) {
               out <- t(sapply(1:nrow(u0), function(i) rmultinom(1, u0[i, cmp], rep(1, n.subnodes))))
               if(n.subnodes == 1) out <- t(out)
               colnames(out) <- paste0(cmp, 1:n.subnodes)
               out
             })
    u0.out <- data.frame(do.call("cbind", u0.list.host))
    u0.out <- cbind(u0.out, data.frame(do.call("cbind", u0.list.mosq)))
    
    stopifnot(all(names(u0.out) %in% sub.comp) & all(sub.comp %in% names(u0.out)))
    u0.out <- u0.out[, sub.comp]
    u0.out$Ic <- u0$Ic
    
    # make the events table. this allows different types of event (e.g. movement, vaccination, etc)
    # to happen to different sets of animals. for example, movements might be taken from all subnodes 
    # and all compartments, while a movement across a border might just happen to edge subnodes.
    # event matrix births in col 1 and deaths in col 2 - everyone born suscepbtilbel so only 1s in S compartments
    # everyone can die so all 1's on col 2. 
    Ev <- matrix(rep(0, n.sub.comp * 3),
                 nrow = n.sub.comp,
                 dimnames = list(sub.comp))
    
    Ev[substr(rownames(Ev), 1, 2) == "Sh", 1] <- 1 # just susceptible cattle
    Ev[substr(rownames(Ev), 1, 2) == "Sh"|
         substr(rownames(Ev), 1, 2) == "Ih"|
         substr(rownames(Ev), 1, 2) == "Rh", 2] <- 1 # all cattle 
    Ev[substr(rownames(Ev), 1, 2) == "Ih", 3] <- 1 # just infected hosts
    
    Ev <- rbind(Ev, Ic = 0)
    
    stopifnot(all(rownames(Ev) %in% names(u0.out)) & all(names(u0.out) %in% rownames(Ev)))
    # output all elements as a list
    list(compartments = as.character(rownames(Ev)), transitions = transitions, Ev = Ev, u0 = u0.out, adj = adj)
  }
