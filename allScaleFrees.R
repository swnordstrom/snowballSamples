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

biasVar = bind_rows(mget(grep("biasVar", ls(), value = TRUE)))
biasVar = biasVar[with(biasVar, order(gr,maxL)),]

# save(errDiv, biasVar, file = "data_original10graphs.RData")

plot(errMean ~ maxL, biasVar)
abline(h = 0)
for (j in 1:10) with(biasVar[biasVar$gr %in% j,], lines(maxL, errMean, col = 'gray66'))

plot(sd ~ sqrt(maxL), biasVar)

# par(mfrow = c(1,1))

plot(div ~ kTrue, qerr[qerr$maxL %in% 1,], type = 'l', col = 'red', log = 'x', 
     ylim = c(-.05, .15), xlab = 'Degree', ylab = 'KLD')
points(div ~ kTrue, qerr[qerr$maxL %in% 2,], type = 'l', col = 'orange')
points(div ~ kTrue, qerr[qerr$maxL %in% 4,], type = 'l', col = 'blue')
points(div ~ kTrue, qerr[qerr$maxL %in% 3,], type = 'l', col = 'green')
abline(h = 0)
legend('topleft', col = c('red', 'orange', 'green', 'blue'),
       pch = '*', cex = 0.8, legend = paste(1:4, 'steps'))

# maybe show divergence using opaque polygons
# annoying problem: degrees are different for each graph (ugh!)
# plot(div ~ kTrue, errDiv[errDiv$maxL %in% 1 & errDiv$gr %in% 1,],
#      log = 'x', ylim = c(range(errDiv$div)[1], range(errDiv$div)[2]), 
#      type = 'l', col = 'white')
# colz = c('red', 'orange', 'green', 'blue')
# for (i in 1:4) {
#   for (j in 2:10) {
#     points(div ~ kTrue, errDiv[errDiv$maxL %in% i & errDiv$gr %in% j,], 
#            type = 'h', col = colz[i])
#   }
# }


# with(errDiv[errDiv$ plot(div ~ )