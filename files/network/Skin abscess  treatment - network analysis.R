library(pcnetmeta)
setwd("../Network")

data.all <- read.table(textConnection('
    Study,        Year,	PMID,	    s.id, r,	n,	 t.id, trt,
    Duam*,       	2017,	28657870,	1,   48,	263,	4, TMP-SMX,
    Duam*,     	  2017,	28657870,	1,   45,	266,	3, Clindamycin,
    Duam*,       	2017,	28657870,	1,   80,	257,	1, Placebo,
    Talan,    	  2016,	26962903,	2,  123,	630,	4, TMP-SMX,
    Talan,    	  2016,	26962903,	2,  163,	617,	1, Placebo,
    Schmitz,  	  2010,	20346539,	3,   15,	 88,	4, TMP-SMX,
    Schmitz,  	  2010,	20346539,	3,   27,	102,	1, Placebo,
    Rajendran,	  2007,	17846141,	4,   13,	 82,	2, Ceph.1st,
    Rajendran,	  2007,	17846141,	4,    8,	 84,	1, Placebo,
    Duong,        2010,	19409657,	5,    3,	 73,	4, TMP-SMX,
    Duong,        2010,	19409657,	5,    4,	 76,	1, Placebo,
    Chen*,        2011, 1339275,  6,    3,   97,  3, Clindamycin,
    Chen*,        2011, 1339275,  6,    6,   97,  2, Ceph.1st,
    '), sep = ",", strip.white = TRUE, header=TRUE)
newdata <- data.all[which(data.all$trt!='ZZZ'),]
sum(newdata$n)
max(newdata$r/newdata$n)
median(newdata$r/newdata$n)
min(newdata$r/newdata$n)
data.adults <- read.table(textConnection('
    Study,        Year,	PMID,	    s.id, r,	n,	 t.id, trt,
    Duam*,       	2017,	28657870,	1,   48,	263,	4, TMP-SMX,
    Duam*,     	  2017,	28657870,	1,   45,	266,	3, Clindamycin,
    Duam*,       	2017,	28657870,	1,   80,	257,	1, Placebo,
    Talan,    	  2016,	26962903,	2,  123,	630,	4, TMP-SMX,
    Talan,    	  2016,	26962903,	2,  163,	617,	1, Placebo,
    Schmitz,  	  2010,	20346539,	3,   15,	 88,	4, TMP-SMX,
    Schmitz,  	  2010,	20346539,	3,   27,	102,	1, Placebo,
    Rajendran,	  2007,	17846141,	4,   13,	 82,	2, Ceph.1st,
    Rajendran,	  2007,	17846141,	4,    8,	 84,	1, Placebo,
    '), sep = ",", strip.white = TRUE, header=TRUE)
newdata <- NULL
newdata <- data.adults[which(data.adults$trt == 'Placebo'),]
sum(newdata$n)
sum(newdata$r)
max(newdata$r/newdata$n)
median <- median(newdata$r/newdata$n)
placebo.rate <- newdata$r/newdata$n
min(newdata$r/newdata$n)
std <- sd(newdata$r/newdata$n)
error <- qnorm(0.975)*std/sqrt(sum(newdata$n))
ci.left <- median-error
ci.right <- median+error
data.pediatrics <- read.table(textConnection('
    Study,        Year,	PMID,	    s.id, r,	n,	 t.id, trt,
    Duong,        2010,	19409657,	5,    3,	 73,	4, TMP-SMX,
    Duong,        2010,	19409657,	5,    4,	 76,	1, Placebo,
    Chen*,        2011, 1339275,  6,    3,   97,  3, Clindamycin,
    Chen*,        2011, 1339275,  6,    6,   97,  2, Ceph.1st,
    '), sep = ",", strip.white = TRUE, header=TRUE)
newdata <- data.pediatrics[which(data.pediatrics$trt == 'Placebo'),]
sum(newdata$n)
max(newdata$r/newdata$n)
median(newdata$r/newdata$n)
min(newdata$r/newdata$n)
#Old studies
#Llera,    	  1985,	3880635, 	7,    1,	 27,	2, Ceph.1st,
#Llera,    	  1985,	3880635, 	7,    1,	 23,	1, Placebo,
#Macfie,   	  1977,	322789,	  8,    0,	 57,	3, Clindamycin,
#Macfie,   	  1977,	322789,	  8,    3,	 41,	1, Placebo,

# arm-based network meta-analysis proposed by Zhang et al (2014).
set.seed(1234)
nma.all.out    <- nma.ab.bin(Study, t.id, r, n, data = data.all, param = c("RR"),
                      trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"),
                      model = "het_cor", higher.better = FALSE, prior.type = "invwishart", n.adapt = 400, n.iter = 100, n.chains = 1, digits = 3)
set.seed(1234)
nma.adults.out <- nma.ab.bin(Study, t.id, r, n, data = data.adults, param = c("RR"),
                      trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"),
                      model = "het_cor", higher.better = FALSE, prior.type = "invwishart", n.adapt = 400, n.iter = 100, n.chains = 1, digits = 3)
set.seed(1234)
nma.pediatrics.out <- nma.ab.bin(Study, t.id, r, n, data = data.pediatrics, param = c("RR"),
                      trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"),
                      model = "het_cor", higher.better = FALSE, prior.type = "invwishart", n.adapt = 400, n.iter = 100, n.chains = 1, digits = 3)

#All
#nma.all.out$AbsoluteRisk$Mean_SD
nma.all.out$AbsoluteRisk$Median_CI
absolute.plot(nma.all.out, alphabetic = FALSE, width = 5, network.name = "Abscess")
nma.all.out$RelativeRisk$Mean_SD
nma.all.out$RelativeRisk$Median_CI
contrast.plot(nma.all.out, width = 5, network.name = "Abscess")  
nma.networkplot(s.id, t.id, data = data.all, title = "Network map: skin abscesses treatment with antibiotics",
                trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"))
nma.all.out <- nma.ab.bin(Study, t.id, r, n, data = data.adults, param= "rank.prob",
                trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"),
                model = "het_cor", n.adapt = 400, n.iter = 100, n.chains = 1)
nma.all.out$TrtRankProb
rank.prob(nma.all.out)

#Adults
nma.adults.out$AbsoluteRisk$Mean_SD
nma.adults.out$AbsoluteRisk$Median_CI
absolute.plot(nma.adults.out, alphabetic = FALSE, width = 5, network.name = "Abscess(Adults)")
nma.adults.out$RelativeRisk$Median_CI
contrast.plot(nma.adults.out, width = 5, network.name = "Abscess(Adults)")  
nma.networkplot(s.id, t.id, data = data.adults, title = "Network map: skin abscesses in adults treatment with antibiotics",
                trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"))
nma.adults.out <- nma.ab.bin(Study, t.id, r, n, data = data.adults, param= "rank.prob",
                      trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"),
                      model = "het_cor", n.adapt = 400, n.iter = 100, n.chains = 1)
nma.adults.out$TrtRankProb
rank.prob(nma.adults.out)

#Pediatrics
nma.pediatrics.out$AbsoluteRisk$Median_CI
absolute.plot(nma.pediatrics.out, alphabetic = FALSE, width = 5, network.name = "Abscess(Pediatrics)")
nma.pediatrics.out$RelativeRisk$Median_CI
contrast.plot(nma.pediatrics.out, width = 5, network.name = "Abscess(Pediatrics)")  
nma.networkplot(s.id, t.id, data = data.pediatrics, title = "Network map: skin abscesses in children treatment with antibiotics",
                trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"))
# Prettify name for Rank chart
nma.pediatrics.out <- nma.ab.bin(Study, t.id, r, n, data = data.pediatrics, param= "rank.prob",
                      trtname = c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"),
                      model = "het_cor", n.adapt = 400, n.iter = 100, n.chains = 1)
nma.pediatrics.out$TrtRankProb
rank.prob(nma.pediatrics.out)

# WORKS
library (metafor)
### calculate log odds for each study arm
dat <- escalc(measure="PLO", xi=r, ni=n, add=1/2, to="all", weights=freq, data=data.adults)
### create network graph (using 'plyr' and 'igraph' packages if installed)
if (require(plyr) && require(igraph)) {
  pairs <- do.call(rbind, sapply(split(dat$t.id, dat$s.id), function(x) t(combn(x,2))))
  pairs <- ddply(data.frame(pairs), .(X1, X2), count)
  g <- graph.edgelist(as.matrix(pairs[,1:2]), directed=FALSE)
  plot(g, edge.curved=FALSE, edge.width=pairs$freq, vertex.label.dist=.7,
       vertex.label=c("Placebo", "Ceph.1st", "Clindamycin", "TMP-SMX"))
}

#WORKS
library(gemtc)
data <- NULL
data$study      <- data.adults$Study
data$treatment  <- as.factor(data.adults$t.id)
data$responders <- data.adults$r
data$sampleSize <- data.adults$n
data <- data.frame(data)
treatments <- read.table(textConnection('
  id description
  1 "Placebo"
  2 "Ceph.1st"
  3 "Clindamycin"
  4 "TMP-SMX"
  '), header=TRUE)
network.adults <- mtc.network(data.ab = data, description="Skin abscess in adults", treatments=treatments)
plot(network.adults)
print(network.adults)
# Analysis of heterogeneity (ANOHE)
hetero <- mtc.anohe(network.adults)
print(hetero)
summary(hetero)
# Hetero from node splitting NOWORK
result.ns <- mtc.nodesplit(network.adults, thin=50)
summary.ns <- summary(result.ns)
print(summary.ns)
plot(summary.ns)
# Summarize the network (generate some interesting network properties)
summary(network.adults)
model <- mtc.model(network.adults)
result <- mtc.run(model, sampler = NA, n.adapt = 5000, n.iter = 20000, thin = 1)
#plot(results)
forest(result, use.description=TRUE, "")
#library(grid)
#grid.text("MAIN TITLE", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
gridtbl <- relative.effect.table(result, covariate=NA)
print(tbl)
