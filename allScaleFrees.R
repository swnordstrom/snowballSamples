library(dplyr)
library(igraph)

load("10scalefree_n740.RData")
source("trialGenerators.R")

set.seed(171872898)

for (x in 1:10) {
  
  assign(paste0("samp",x), budgetedTrial(graphs[[x]], 150, 50))
  
  results = calcErrorDivergence(get(paste0("samp", x)), graphs[[x]], 50, 150)
  assign(paste0("indivs",x), cbind(gr = rep(x, nrow(results$indivs)), results$indivs))
  assign(paste0("errDiv",x), cbind(gr = rep(x, nrow(results$errDiv)), results$errDiv))
  
  assign(paste0("biasVar", x), 
         aggregate(get(paste0("indivs",x))$kSamp, 
                   by = list(trial = get(paste0("indivs",x))$trial), mean) %>%
           merge(unique(indivs1[,c("trial", "maxL")]), by = "trial") %>%
           rename(mean.deg = x) %>%
           group_by(maxL) %>%
           summarise(mean = mean(mean.deg), sd = sd(mean.deg)) )
  assign(paste0("biasVar",x),
         cbind(gr = rep(x, nrow(get(paste0("biasVar",x)))), 
               get(paste0("biasVar",x)),
               errMean = get(paste0("biasVar",x))$mean - mean(degree(graphs[[x]]))))
  
}

errDiv = bind_rows(mget(grep("errDiv", ls(), value = TRUE)))

allIndivs = bind_rows(mget(grep("indivs", ls(), value = TRUE)))
allIndivs = allIndivs[with(allIndivs, order(gr,node)),]

biasVar = bind_rows(mget(grep("biasVar", ls(), value = TRUE)))
biasVar = biasVar[with(biasVar, order(gr,maxL)),]

# save(allIndivs, errDiv, biasVar, file = "data_original10graphs.RData")
load('data_original10graphs.RData')

# plot mean in error estimates for all ten networks

plot(errMean ~ maxL, biasVar,
     xlab = "Steps per snowball",
     ylab = "Mean smaple degree - true mean")
abline(h = 0)
for (j in 1:10) with(biasVar[biasVar$gr %in% j,], lines(maxL, errMean, col = 'gray66'))

# plot standard deviation in mean-degree estimate

plot(sd ~ maxL, biasVar,
     xlab = "Steps per snowball", 
     ylab = "Std. Dev. of mean sample degree")

## Plot divergence + loess fit for all ten graphs

colz = c("red", "orange", "green", "blue")
plot(div ~ kTrue, errDiv,
     pch = 19, cex = 0.15, col = colz[maxL],
     log = 'x')
abline(h = 0)

divLo = loess(div ~ kTrue + maxL, errDiv, span = 1)
prdcDivLo = predict(divLo, expand.grid(maxL = 1:4, kTrue = 1:300))
for (i in 1:4) points(prdctDivLo[i,], type = 'l', col = colz[i], lwd = 2)
legend('topleft', cex = 0.5, col = colz, 
       pch = 19, bty = "n", 
       legend = paste(1:4, 'step'))

## Summary plot of individual error for all four maxL settings

allIndivs = allIndivs[allIndivs$kTrue > 1,]

aierr = allIndivs %>%
  group_by(maxL, kTrue) %>%
  summarise(median = median(err), q25 = quantile(err, .25), q75 = quantile(err, .75))

par(mfrow = c(2,2))

for (i in 1:4) {
  
  plot(2, 1, xlim = c(2, max(aierr$kTrue)), ylim = c(min(aierr$q25), 1), type = 'n')
  with(aierr[aierr$maxL %in% i,], 
       polygon(x = c(kTrue, rev(kTrue)), y = c(q75, rev(q25)), 
               col = 'gray88', border = 'gray88'))
  points(median ~ kTrue, aierr[aierr$maxL %in% i,],
         type = 'l', lwd = 2)
    
}

#### Begin heuristics

## halo degree heuristic

trialsPer = 50
thresh = (3:7)/10

for (x in 1:10) {

  assign(paste0("samp",x), haloThreshold(graphs[[x]], thresh, trialsPer, 150))
  
  assign(paste0("sampDf",x), data.frame(gr = x, th = rep(thresh, each = trialsPer),
                                        mean = sapply(get(paste0("samp",x)),
                                                      function(x) mean(x$kSamp)),
                                        sd = sapply(get(paste0("samp",x)),
                                                    function(x) sd(x$kSamp))))
  
}

haloThreshSamps = bind_rows(mget(grep("sampDf", ls(), value = TRUE)))
trueMean = data.frame(gr = 1:10, trueMean = sapply(graphs, function(x) mean(degree(x))))
haloThreshSamps = haloThreshSamps %>%
  merge(trueMean, by = "gr") %>%
  mutate(errMean = mean - trueMean)

# save(haloThreshSamps, file = "haloThreshSampleData.RData")

plot(errMean ~ jitter(th), haloThreshSamps,
     xlab = "Mean halo degree threshold (jittered)",
     ylab = "Mean smaple degree - true mean")
haloThreshloess = loess(errMean ~ th, haloThreshSamps)
points(thresh, unique(haloThreshloess$fitted), type = 'l', col = 'red')

