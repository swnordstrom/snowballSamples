library(igraph)
source('trialGenerators.R')

export.dir = "/users/scottnordstrom/documents/Homework/rotation 1/presentation/"

set.seed(9464841)

rand.er = erdos.renyi.game(100, .04, 'gnp')
rand.er = delete.vertices(rand.er, V(rand.er)[!degree(rand.er)])

png(paste0(export.dir, 'fullGraph.png'))
plot.igraph(rand.er, vertex.size = 2, vertex.label = NA)
dev.off()

grel = as_edgelist(rand.er)
deg.gr = degree(rand.er)
deg.df = data.frame(name = as.numeric(V(rand.er)), deg = deg.gr)

### A function for initializing sample

initSample = function(gr, grel) {
  
  e.l = matrix(nrow = 0, ncol = 2)
  # randomly select seed node
  sn = sample(V(gr), 1) 
  # begin first snowball by taking step out from seed
  e.l = stepOut(grel, e.l, sn)
  # remove the seed from the halo
  e.l$halo = e.l$halo[e.l$halo != sn]
  
  return(e.l)
  
}

### A function to clean up sample (crop to appropriate size)

cropSample = function(e.l, targ) {
  
  samp = e.l
    
  toRm = sample(samp$halo, length(unique(array(samp$g))) - targ)
  e.l$g = samp$g[!(samp$g[,1] %in% toRm | samp$g[,2] %in% toRm),]
  
  samp.g = graph.edgelist(samp$g, directed = FALSE)
  samp.g = delete.vertices(samp.g, !degree(samp.g))
  
  return(samp.g)
  
}

###
### Generate uniform vertex sample
###

e.l.u = initSample(rand.er, grel)

while (length(unique(array(e.l.u$g))) < 40) {
  incompl = merge(deg.df, data.frame(table(e.l.u$g)), 
                  by.x = 'name', by.y = 'Var1', all = TRUE)
  # randomly choose a seed node from the set of incompletely-sampled nodes
  sn = as.numeric(sample(incompl$name[incompl$deg != incompl$Freq |
                                        is.na(incompl$Freq)], 1))
  e.l.u = stepOut(grel, e.l.u$g, sn)
  e.l.u$halo = e.l.u$halo[e.l.u$halo != sn]
}  

# removing nodes to get us to the target size
samp.unif = cropSample(e.l.u, 40)

png(paste0(export.dir, 'uniformSample.png'))
plot.igraph(samp.unif, vertex.size = 2, vertex.label = NA)
dev.off()

###
### generate snowball sample
###

e.l.s = initSample(rand.er, grel)

while (length(unique(array(e.l.s$g))) < 40) {
  if (length(e.l.s$halo)) {
    e.l.s = stepOut(grel, e.l.s$g, e.l.s$halo)
  } else {
    incompl = merge(deg.df, data.frame(table(e.l.s$g)), 
                    by.x = 'name', by.y = 'Var1', all = TRUE)
    # randomly choose a seed node from the set of incompletely-sampled nodes
    sn = as.numeric(sample(incompl$name[incompl$deg != incompl$Freq |
                                          is.na(incompl$Freq)], 1))
    e.l.s = stepOut(grel, e.l.s$g, sn)
    e.l.s$halo = e.l.s$halo[e.l.s$halo != sn]
  }
}  

# removing nodes to get us to the target size
samp.snow = cropSample(e.l.s, 40)

png(paste0(export.dir, 'snowball.png'))
plot.igraph(samp.snow, vertex.size = 2, vertex.label = NA)
dev.off()

###
### two step-sample
###

e.l.2step = initSample(rand.er, grel)
steps = 1

while (length(unique(array(e.l.2step$g))) < 40) {
  if (steps < 2 & length(e.l.2step$halo)) {
    e.l.2step = stepOut(grel, e.l.2step$g, e.l.2step$halo)
    steps = steps + 1
  } else {
    incompl = merge(deg.df, data.frame(table(e.l.2step$g)), 
                    by.x = 'name', by.y = 'Var1', all = TRUE)
    # randomly choose a seed node from the set of incompletely-sampled nodes
    sn = as.numeric(sample(incompl$name[incompl$deg != incompl$Freq |
                                          is.na(incompl$Freq)], 1))
    e.l.2step = stepOut(grel, e.l.2step$g, sn)
    e.l.2step$halo = e.l.2step$halo[e.l$halo != sn]
    steps = 1
  }
}  
    
samp.2step = cropSample(e.l.2step, 40)

png(paste0(export.dir, 'twoStep.png'))
plot.igraph(samp.2step, vertex.size = 2, vertex.label = NA)
dev.off()
