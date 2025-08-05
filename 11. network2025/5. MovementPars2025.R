# import movement parameters ----

# load predicted probability of movement
mu.z.list <-
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    mu.z <- 
      as.matrix(
        read.csv(paste0("cattle.month01-12.mvt.2018-06-07-00/", sp, ".month", m, ".mvt.matrix.prob.2018-06-07-00.csv"), 
                 row.names = 1))
    # make rownames and colnames consistent                  
    colnames(mu.z) <- rownames(mu.z)
    # keep movements to and from included wards
    mu.z[ward.names, ward.names]
  })


# load predicted number moved given movement
mu.c.list <-
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    mu.c <- 
      as.matrix(
        read.csv(paste0("cattle.month01-12.mvt.2018-06-07-00/", sp, ".month", m, ".mvt.matrix.rate.2018-06-07-00.csv"), 
                 row.names = 1))
    # make rownames and colnames consistent                  
    colnames(mu.c) <- rownames(mu.c)
    # keep movements to and from included wards
    mu.c[ward.names, ward.names]
  })

# load dispersion (theta) parameter for the ZTNB distribution
theta.list <- 
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    read.csv(paste0("cattle.month01-12.mvt.2018-06-07-00/", sp, ".month", m, ".mvt.matrix.rate.theta.2018-06-07-00.csv"))$theta
  })

names(mu.z.list) <- names(mu.c.list) <- names(theta.list) <- all.months

# calculate expected livestock flow for the whole time period (all.months)
mu.yr <- Reduce("+", lapply(all.months, function(m) mu.z.list[[m]] * mu.c.list[[m]]))

# define a market as a ward with expected outflow >0 for the year
mkts.all <- unique(c(mkts, rownames(mu.yr)[rowSums(mu.yr) != 0]))

# estimate betweenness
between <- wards$betweenness
names(between) <- rownames(wards)

# estimate outdegree
outdeg <- wards$year.out.deg
names(outdeg) <- rownames(wards)


# estimate indegree
indeg <- wards$year.in.deg
names(indeg) <- rownames(wards)

sort(indeg)

# NB, in Phil Trans paper:
# betweenness = mp.between
# degree centrality = deg.cent
# eigenvector centrality = mp.eigval

# select wards for targeting interventions

# those with highest RVF risk (SPREAD, NOT EMERGENCE)
rvf.des <- ward.names[order(-wards[ward.names, "rvf.risk2"])]
mkts %in% rvf.des[1:20]

rvf.des[1:20] %in% not.mkts.rvf.risk.tab[, 1]
not.mkts.rvf.risk.tab[, 1]  %in% rvf.des[1:20]

# check they are ordered by RVF risk
wards[rvf.des[1:5], "rvf.risk2"]

# how many of the seed wards are in the RVF top 20?
not.mkts.rvf.risk.tab[not.mkts.rvf.risk.tab %in% rvf.des[1:20]]

# those with most animals (change to density, once the model includes density?)
n.animals.des <- ward.names[order(-wards[ward.names, paste0("pop.", sp)])]
mkts %in% n.animals.des[1:80]
# check they are ordered by no of animals
wards[n.animals.des[1:5], paste0("pop.", sp)]

# highest betweenness 
between.des <- names(rev(sort(between[ward.names])))
mkts %in% between.des[1:80]
# check they are ordered by betweenness
wards[between.des[1:5], "betweenness"]

# highest betweenness combined with highest RVF SPREAD risk
between.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk2) + wards$betweenness.rank)]
#  c(between.des[between.des %in% wards$fullname[wards$rvf.risk.bin == 1]],
#    between.des[between.des %in% wards$fullname[wards$rvf.risk.bin == 0]])
# check they are ordered by betweenness and rvf risk
wards[between.rvf.des[1:5], c("betweenness", "rvf.risk")]

# highest outdegree 
outdegree.des<- names(rev(sort(outdeg[ward.names])))
mkts %in% outdegree.des[1:80]
# check they are ordered by outdegree
wards[outdegree.des[1:5], "year.out.deg"]

# highest outdegree combined with highest RVF SPREAD risk
outdegree.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk2) + wards$year.out.deg.rank)]
# check they are ordered by outdegree and rvf risk
wards[outdegree.rvf.des[1:5], c("year.out.deg", "rvf.risk")]

# highest indegree
indegree.des <- names(rev(sort(indeg[ward.names])))
mkts %in% indegree.des[1:80]
# check they are ordered by indegree
wards[indegree.des[1:5], "year.in.deg"]

# highest indegree combined with highest RVF SPREAD risk
indegree.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk2) + wards$year.in.deg.rank)]
# check they are ordered by indegree and rvf risk
wards[indegree.rvf.des[1:5], c("year.in.deg", "rvf.risk")]

# highest combined multiplex rank
mp.comb.des <- rownames(wards)[order(wards$mp.comb.rank)]
mkts %in% mp.comb.des[1:80]
# check they are ordered by multiplex rank
wards[mp.comb.des[1:5], "mp.comb.rank"]

# highest combined multiplex rank combined with highest RVF SPREAD risk
mp.comb.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk2) + wards$mp.comb.rank)]
# check they are ordered by combined multiplex rank and rvf risk
wards[mp.comb.rvf.des[1:5], c("mp.comb.rank", "rvf.risk")]

# highest MP geometric mean degree (REFERRED TO AS DEGREE CENTRALITY IN PHIL TRANS PAPER)
deg.cent.des <- rownames(wards)[order(-wards$deg.cent)]
mkts %in% deg.cent.des[1:80]
mkts %in% deg.cent.des[1:20]
table(mkts.all %in% deg.cent.des[1:20])
# check they are ordered by geometric mean degree
wards[deg.cent.des[1:5], "deg.cent"]

# how many of the seed wards are in the deg.cent top 20?
not.mkts.rvf.risk.tab[not.mkts.rvf.risk.tab %in% deg.cent.des[1:20]]

# highest MP geometric mean degree combined with highest RVF SPREAD risk
deg.cent.rvf.des <-
  ward.names[order(rank(-wards$rvf.risk2) + wards$deg.cent.rank)]
# check they are ordered by indegree and rvf risk
wards[deg.cent.rvf.des[1:5], c("deg.cent", "rvf.risk")]

# how many of the seed wards are in the deg.cent.rvf top 20?
not.mkts.rvf.risk.tab[not.mkts.rvf.risk.tab %in% deg.cent.rvf.des[1:20]]

# highest MP betweenness
mp.between.des <- rownames(wards)[order(-wards$mp.betweenness)]
mkts %in% mp.between.des[1:80]
mkts %in% mp.between.des[1:20]
table(mkts.all %in% mp.between.des[1:20])

# check they are ordered by MP betweenness
wards[mp.between.des[1:5], "mp.betweenness"]


# highest MP eigenvalue centrality
mp.eigval.des <- rownames(wards)[order(-wards$mp.eigval)]
mkts %in% mp.eigval.des[1:80]
mkts %in% mp.eigval.des[1:20]
table(mkts.all %in% mp.eigval.des[1:20])
# check they are ordered by eigenvalue centrality
wards[mp.eigval.des[1:5], "mp.eigval"]


# how similar are the top wards selected by
#   mp eigenvalue centrality?
#   mp betweenness?
#   mp geomean degree?

n.compare <- 20
unique(c(mp.between.des[1:n.compare], deg.cent.des[1:n.compare], mp.eigval.des[1:n.compare]))

# between and geomean degree differ on 9 wards in top 20
setdiff(mp.between.des[1:n.compare], between.des[1:n.compare])
setdiff(between.des[1:n.compare], mp.between.des[1:n.compare])

# between and geomean degree differ on 7 wards in top 20
setdiff(mp.between.des[1:n.compare], deg.cent.des[1:n.compare])
setdiff(deg.cent.des[1:n.compare], mp.between.des[1:n.compare])

# between and eigval differ on 15 wards in top 20
setdiff(mp.between.des[1:n.compare], mp.eigval.des[1:n.compare])
setdiff(mp.eigval.des[1:n.compare], mp.between.des[1:n.compare])

# geomean and eigval differ on 10 wards in top 20
setdiff(deg.cent.des[1:n.compare], mp.eigval.des[1:n.compare])
setdiff(mp.eigval.des[1:n.compare], deg.cent.des[1:n.compare])
