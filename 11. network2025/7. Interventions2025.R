# movement ban starts where and when? (use day=Inf for no ban, start=1 for no movement ever)
# note that the movement ban happens before contiguous movements, so will not include them
# if contiguous movements are also being banned, set contig.move.n <- 0

# target how many wards?
n.ward <- c(few = 20, many = 80, all = n)

# how many months does the movement ban apply to? (assume it is reactive)
mban.months <- all.months[-1]

# list of the movement ban strategies ----
mban.list <- list(
  none = list(ward = "none", mban.months = "none"),
  all = list(ward = ward.names[ward.names %in% mkts.all], mban.months = mban.months),
  rand = list(ward = sample(mkts.all), mban.months = mban.months),
  between = list(ward = between.des[between.des %in% mkts.all], mban.months = mban.months),
  indegree = list(ward = indegree.des[indegree.des %in% mkts.all], mban.months = mban.months),
  outdegree =  list(ward = outdegree.des[outdegree.des %in% mkts.all], mban.months = mban.months),
  mp.comb = list(ward = mp.comb.des[mp.comb.des %in% mkts.all], mban.months = mban.months),
  mp.eigval = list(ward = mp.eigval.des[mp.eigval.des %in% mkts.all], mban.months = mban.months),
  mp.between = list(ward = mp.between.des[mp.between.des %in% mkts.all], mban.months = mban.months),
  deg.cent = list(ward = deg.cent.des[deg.cent.des %in% mkts.all], mban.months = mban.months),
  rvf = list(ward = rvf.des[rvf.des %in% mkts.all], mban.months = mban.months),
  between.rvf = list(ward = between.rvf.des[between.rvf.des %in% mkts.all], mban.months = mban.months),
  indegree.rvf = list(ward = indegree.rvf.des[indegree.rvf.des %in% mkts.all], mban.months = mban.months),
  outdegree.rvf = list(ward = outdegree.rvf.des[outdegree.rvf.des %in% mkts.all], mban.months = mban.months),
  mp.comb.rvf = list(ward = mp.comb.rvf.des[mp.comb.rvf.des %in% mkts.all], mban.months = mban.months),
  deg.cent.rvf = list(ward = deg.cent.rvf.des[deg.cent.rvf.des %in% mkts.all], mban.months = mban.months),
  n.animals = list(ward = n.animals.des[n.animals.des %in% mkts.all], mban.months = mban.months))

# list of vaccine details ----
# where, when and what proportion to vaccinate? use day=Inf or cov=0 for never

vax.list <- list(none = data.frame(ward = "none", day = Inf, cov = 0, stringsAsFactors = FALSE),
                 all.lo = data.frame(ward = ward.names, day = 0, cov = 0.1, stringsAsFactors = FALSE, row.names = NULL),
                 all.hi = data.frame(ward = ward.names, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 rand = data.frame(ward = sample(ward.names), day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 between = data.frame(ward = between.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 mp.between = data.frame(ward = mp.between.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 indegree = data.frame(ward = indegree.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 outdegree = data.frame(ward = outdegree.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 mp.comb = data.frame(ward = mp.comb.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 deg.cent = data.frame(ward = deg.cent.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 deg.cent.rvf = data.frame(ward = deg.cent.rvf.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 mp.eigval= data.frame(ward = mp.eigval.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 rvf = data.frame(ward = rvf.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 between.rvf = data.frame(ward = between.rvf.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 indegree.rvf = data.frame(ward = indegree.rvf.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 outdegree.rvf = data.frame(ward = outdegree.rvf.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 mp.comb.rvf = data.frame(ward = mp.comb.rvf.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 n.animals = data.frame(ward = n.animals.des, day = 0, cov = inf.pars$cov, stringsAsFactors = FALSE, row.names = NULL),
                 mkts.monthly = list(ward = ward.names[is.mkt], day = all.months.day1 + 3, cov = inf.pars$cov, stringsAsFactors = FALSE))

# allow vaccine coverage to vary among wards according to RVF risk
if(FALSE) {
  for(vn in names(vax.list)[!names(vax.list) %in% c("none", "all.lo")]) {
    vax.list[[vn]]$cov <- wards[vax.list[[vn]]$ward, "cov"] 
  }; rm(vn)
}

# for convenience, list likely combinations of vax and mban strategies to trial ----
vax.sets <-
  list(
    c("none"),
    c("none", "all.hi"),
    c("none", "all.hi", "all.lo", "rand", "deg.cent", "rvf", "deg.cent.rvf", "n.animals"), # used for WOHC2020
    c("none", "all.hi", "mp.between", "deg.cent", "mp.eigval", "rand"),
    c("none", "all.hi", "mp.between", "deg.cent", "mp.eigval", "rand", "n.animals"),
    c("none", "all.hi", "all.lo", "rand", "between", "indegree", "mp.comb"),
    c("none", "all.hi", "all.lo", "rand"),
    c("none", "all.hi", "between", "indegree", "mp.comb"),
    c("none", "all.hi", "rvf"),
    c("none", "all.hi", "between.rvf", "indegree.rvf", "mp.comb.rvf", "rvf"))

vax.plans <- unique(unlist(vax.sets))

mban.sets <- 
  list(
    c("none"),
    c("none", "all"),
    c("none", "all", "rand", "deg.cent", "rvf", "deg.cent.rvf", "n.animals"), # used for WOHC2020
    c("none", "all", "mp.between", "deg.cent", "mp.eigval", "rand"),
    c("none", "all", "mp.between", "deg.cent", "mp.eigval", "rand", "n.animals"),
    c("none", "all", "rand", "between", "indegree", "mp.comb"),
    c("none", "all", "rand", "rvf"),
    c("none", "all", "rand", "between.rvf", "indegree.rvf", "mp.comb.rvf", "rvf"),
    "all")

mban.plans <- unique(unlist(mban.sets))

# intervention table ----
intervention.vax <- 
  expand.grid(vax = vax.sets[[1]]
              ,mban = "none"
              ,n.vax = n.ward["few"]
              ,n.mban = 1
              ,recovery  = inf.pars$recovery #0.2 # host recovery
              ,muH  = inf.pars$muH #0.000456621 #0.005479452 gives a steady infection but 100% sero+ when not differentiating birth and death rate
              ,muBirth = inf.pars$muBirth # 0.001369863
              ,mu = inf.pars$mu # 0.03058104 # https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-023-05792-3#:~:text=The%20average%20adult%20lifespan%20for,Culex%20species%20(Table%202)
              ,p_hv = inf.pars$p_hv #0.25 # unknwown
              ,p_vh = inf.pars$p_vh #0.04 # unknown
              ,laying =inf.pars$laying #20000 # rafts of 300 eggs
              ,extrinInc = inf.pars$extrinInc #0.1094668 # unknown extrinsic incubation - some sources say as little as 1 day in high temps
              ,development = inf.pars$development #0.1  # eggs take 10 days to develop into adults
              ,K = inf.pars$K #1000000#110000
              ,Kh = inf.pars$Kh #150000#1200
              ,biteRate = inf.pars$biteRate #0.68
              ,firstBiteDelay = inf.pars$firstBiteDelay #0.14
              ,pupsMort = inf.pars$pupsMort #0.8
              ,rep = 1:nrep
              ,stringsAsFactors = FALSE)

intervention.mban <- 
  expand.grid(vax = "none"
              ,mban = mban.sets[[1]] 
              ,n.vax = 1
              ,n.mban = n.ward["few"]
              ,recovery  = inf.pars$recovery #0.2 # host recovery
              ,muH  = inf.pars$muH #0.000456621 #0.005479452 gives a steady infection but 100% sero+ when not differentiating birth and death rate
              ,muBirth = inf.pars$muBirth # 0.001369863
              ,mu = inf.pars$mu # 0.03058104 # https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-023-05792-3#:~:text=The%20average%20adult%20lifespan%20for,Culex%20species%20(Table%202)
              ,p_hv = inf.pars$p_hv #0.25 # unknwown
              ,p_vh = inf.pars$p_vh #0.04 # unknown
              ,laying =inf.pars$laying #20000 # rafts of 300 eggs
              ,extrinInc = inf.pars$extrinInc #0.1094668 # unknown extrinsic incubation - some sources say as little as 1 day in high temps
              ,development = inf.pars$development #0.1  # eggs take 10 days to develop into adults
              ,K = inf.pars$K #1000000#110000
              ,Kh = inf.pars$Kh #150000#1200
              ,biteRate = inf.pars$biteRate #0.68
              ,firstBiteDelay = inf.pars$firstBiteDelay #0.14
              ,pupsMort = inf.pars$pupsMort #0.8
              ,rep = 1:nrep
              ,stringsAsFactors = FALSE)

intervention <- rbind(intervention.vax, intervention.mban)

if(all(intervention$vax == "none")) 
  intervention$plan <- intervention$mban else
    if(all(intervention$mban == "none")) intervention$plan <- intervention$vax else
      intervention$plan <- paste0("VX: ", intervention$vax, "; MB: ", intervention$mban)

intervention$n.vax[intervention$vax %in% "none"] <- 0
intervention$n.vax[intervention$vax %in% c("all.lo", "all.hi")] <- n
intervention$n.mban[intervention$mban %in% "none"] <- 0
intervention$n.mban[intervention$mban %in% "all"] <- length(mkts.all)
intervention <- intervention[!duplicated(intervention), ]
dim(intervention)
intervention$k <- 1:nrow(intervention)
intervention$contig.move.n <- contig.move.n