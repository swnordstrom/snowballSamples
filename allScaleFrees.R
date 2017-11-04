load("10scalefree_n740.RData")
source("trialGenerators.R")

set.seed(171872898)

for (x in 1:10) {
  
#  assign(paste0("samp",x), budgetedTrial(graphs[[x]], 150, 50))
  
#  results = calcErrorDivergence(get(paste0("samp", x)), graphs[[x]], 50, 150)
#  assign(paste0("indivs",x), cbind(gr = rep(x, nrow(results$indivs)), results$indivs))
#  assign(paste0("errDiv",x), cbind(gr = rep(x, nrow(results$errDiv)), results$errDiv))
  
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

biasVar = bind_rows(mget(grep("biasVar", ls(), value = TRUE)))
biasVar = biasVar[with(biasVar, order(gr,maxL)),]

plot(errMean ~ maxL, biasVar)
abline(h = 0)
for (j in 1:10) with(biasVar[biasVar$gr %in% j,], lines(maxL, errMean, col = 'gray66'))

plot(sd ~ sqrt(maxL), biasVar)
