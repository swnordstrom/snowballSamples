library(igraph)

makeP = function(d, cumulative = FALSE) {
  
  n = as.numeric(names(d))
  missing = setdiff(min(n):max(n), n)
  if (any(missing)) {
    d = c(d, rep(0, length(missing)))
    names(d)[names(d) == ""] = missing
    d = d[order(as.numeric(names(d)))]
  } 
  
  if (cumulative) return( cumsum(d) / max(cumsum(d)))
  else return(d)
  
}

evenStubs = function(x) return(crossprod(makeP(table(x)), 1:length(makeP(table(x))))[1] %% 2)

source("config_graph.R")

set.seed(143220311)

nNodes = 900
nGraphs = 10
graphs = list()
degs = list()

for (i in 1:nGraphs) {
  
  deg = round(1 / (1 - runif(nNodes)))
  while(any(deg > nNodes) | evenStubs(deg)) deg = round(1 / (1 - runif(nNodes)))
  
  print(i)
  graphs[[i]] = configGraph(makeP(table(deg)))
  degs[[i]] = deg
}

save(degs, graphs, file = "10scalefree_n740.RData")

# cauchy-sequence graphs
# set.seed(112911)
# 
# nNodes = 150
# nGraphs = 30
# graphs = list()
# 
# for (i in 1:nGraphs) {
# 
#   deg = ceiling(abs(rcauchy(nNodes, scale = 1)))
#   while(any(deg > nNodes) | evenStubs(deg)) deg = ceiling(abs(rcauchy(nNodes, scale = 1)))
#     
#   print(i)
#   graphs[[i]] = configGraph(makeP(table(deg)))
#   
# }
#
# save(graphs, file = "30cauchy_n150.RData")
