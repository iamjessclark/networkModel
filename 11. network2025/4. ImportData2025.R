# import wards data - required to identify the nearest primary market ----
wards <- 
  read.csv('./networkModelData/WardSpatialData.csv', 
           stringsAsFactors = FALSE)

# drop tanga
wards <- droplevels(wards[wards$region %in% c("Arusha", "Manyara", "Kilimanjaro"), ])
all.ward.names <- wards$fullname

  # wards to include in the analysis: all wards enough animals ----
  ### check sensitivity to choice of minimum number (500, 1000, etc)
  min.ward.pop <- 1000
  wards[wards[, paste0("pop.", sp)] < min.ward.pop, paste0("pop.", sp)] <- min.ward.pop
  ward.names <- wards$fullname
  rownames(wards) <- ward.names
  
  all(all.ward.names == ward.names)

# import network measures for wards ----
netmeasures <-
  read.csv("./networkModelData/AllWards.csv",
           stringsAsFactors = FALSE, quote = "", 
           row.names = 1)
head(netmeasures)

# import list of contiguous wards ----
#, to allow spread via local diffusion as well as 
# via the movement network
contig.wards.all <- 
  read.csv(".//networkModelData/contiguous.wards.csv", 
           stringsAsFactors = FALSE, quote = "")
rownames(contig.wards.all) <- paste(contig.wards.all$from, contig.wards.all$to, sep = ".")


# epidemic dates ----    
# the date the epidemic starts
start.date <- as.Date("2015-01-01")

# runs through which months?
# this is only running for a year so no concern of leap years 
n.years <-  1
n.months <- n.years*12
end.date <- start.date %m+% months(n.months)
obvs.date <- start.date+1 #end.date-years(5) 
n.days <- -(as.numeric(difftime(start.date, end.date, units = "days")))
obvs.days <- n.days-as.numeric(difftime(end.date, obvs.date, units = "days"))# this has a negative on it as it is a negative number that comes out of diff time
all.months <- formatC(1:n.months, width = 2, flag = "0")

# which day of the year (1:365) is the first of each month?

all.months.day1 <-
  sapply(all.months, function(mm) {
    print(mm)
    yyyy <- 2014 + ceiling(as.numeric(mm)/12)
    mm <- ((as.numeric(mm) - 1) %% 12) + 1
    as.Date(paste(yyyy, mm, "01", sep = "-")) - start.date + 1
  })

# make risk score for RVF spread ----
#(as opposed to emergence, which is defined as proportion of
# ward with NDVI between 0.15 and 0.4 and exists in wards DF as wards$rvf.risk).
# there is a negative correlation 
# Will de G: "I think we did find a strong and negative relationship with continuous NDVI
# for seroprevalence risk" - he's looking into this to confirm.
wards$rvf.risk2 <- 1 - wards$ndvi

# which wards contain secondary markets?
# these wards are considered as all market, so they have no population of animals resident 
# so there shouldn't be much time for transmission
is.mkt <- ward.names %in% c("arusha/monduli/meserani", "arusha/arusha/bwawani", 
                            "kilimanjaro/hai/machamekusini")
mkts <- ward.names[is.mkt]
not.mkts <- ward.names[!is.mkt]

# network measures ----

# make geometric mean of mp.indegree and mp.outdegree
netmeasures$deg.cent <- 
  exp(apply(log(netmeasures[, c("mp.indegree", "mp.outdegree")]), 1, mean))


# create ranks for selected measures

# sets of closely correlated measures:
#   year.out.deg, mp.outdegree, mp.degree, cy.outdegree, cy.degree
#   year.in.degree, mp.indegree, cy.indegree
#   betweenness, cy.betweenness

meas.to.rank <-
  c("mp.eigval", "mp.indegree",	"mp.outdegree", "mp.betweenness", "betweenness", 
    "year.in.deg", "year.out.deg", "deg.cent")
for(x in meas.to.rank) netmeasures[, paste0(x, ".rank")] <- rank(-netmeasures[, x])
rm(x)

# create sum of ranks for selected measures 
netmeasures$mp.comb.rank <- rank(rowSums(netmeasures[, paste0(meas.to.rank[1:4], ".rank")]))

# restrict to included wards
netmeasures <- netmeasures[ward.names, ]
setdiff(rownames(netmeasures), wards$fullname)
setdiff(wards$fullname, rownames(netmeasures))
#plot(netmeasures[, c("betweenness.rank", "year.in.deg.rank", "mp.comb.rank")])

# add net measures to the wards data frame
wards <- cbind(wards, netmeasures[rownames(wards), ])


# choose top 5 rvf risk wards (FOR EMERGENCE, NOT SPREAD) ----
# by combining out.degree and rvf.risk
#not.mkts.rvf.risk <- 
#  not.mkts[order(-rank(wards[not.mkts, "year.out.deg"]) * rank(wards[not.mkts, "rvf.risk"]))][1:20]
#not.mkts.rvf.risk <- 
#  not.mkts[order(-rank(wards[not.mkts, "rvf.risk"]))][1:20]
not.mkts.rvf.risk <- sample(not.mkts, 5)
not.mkts.rvf.risk.tab <- matrix(not.mkts.rvf.risk, ncol = 1)
