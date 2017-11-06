library(igraph)
library(dplyr)

stepOut = function(true.el, cur.el, halo) {
  ### Function for taking a step "out" within a snowball
  e.l.s = rbind(cur.el, true.el[true.el[,1] %in% halo | true.el[,2] %in% halo,])
  e.l.s = unique(e.l.s)
  halo = setdiff(unique(as.vector(e.l.s)), unique(as.vector(cur.el)))
  
  return(list(g = e.l.s, halo = halo))
}


budgetedTrial = function(gr, targ, trialsPer) {

  li = list()
  grel = as_edgelist(gr)
  deg.gr = degree(gr)
  deg.df = data.frame(name = as.numeric(names(deg.gr)), deg = deg.gr)
  
  for (i in 1:4) {
    
    for (trial in 1:trialsPer) {
      
      e.l = matrix(nrow = 0, ncol = 2)
      
      # randomly select seed node
      sn = sample(V(gr)$name, 1) 
      # begin first snowball by taking step out from seed
      e.l = stepOut(grel, e.l, sn)
      # remove the seed from the halo
      e.l$halo = e.l$halo[e.l$halo != sn]
      steps = 1
      balls = 1
      
      while (length(unique(array(e.l$g))) < targ) {
        if (steps < i& length(e.l$halo)) {
          e.l = stepOut(grel, e.l$g, e.l$halo)
          steps = steps + 1
        } else {
          incompl = merge(deg.df, data.frame(table(e.l$g)), 
                          by.x = 'name', by.y = 'Var1', all = TRUE)
          # randomly choose a seed node from the set of incompletely-sampled nodes
          sn = as.numeric(sample(incompl$name[incompl$deg != incompl$Freq |
                                                is.na(incompl$Freq)], 1))
          e.l = stepOut(grel, e.l$g, sn)
          e.l$halo = e.l$halo[e.l$halo != sn]
          steps = 1
          balls = balls + 1
        }
      }  
      
      # removing nodes to get us to the target size
      toRm = sample(e.l$halo, length(unique(array(e.l$g))) - targ)
      e.l$g = e.l$g[!(e.l$g[,1] %in% toRm | e.l$g[,2] %in% toRm),]
      
      # degree distribution of the sampled graph
      deg = table(e.l$g)
      
      li[[50*(i-1) + trial]] = list(stats = data.frame(maxL = i, balls = balls),
                                    kSamp = deg, 
                                    trueDeg = deg.gr[names(deg)])
      
    }
    
  }
  
  return(li)

}

calcDiv = function(trials, true) {
  true = table(true)
  p = trials / true
  p = p / sum(p)
  q = true / sum(true)
  ret = p * log(p / q)
  ret[is.nan(ret)] = 0
  return(ret)
}

calcErrorDivergence = function(dataset, gr, trialsPer, targ) {
  
  deg.gr = degree(gr)
  
  df = unlist(sapply(dataset, function(x) return(x$kSamp)))
  df = data.frame(node = as.numeric(names(df)), kSamp = df)
  df = data.frame(df, kTrue = unlist(sapply(dataset, function(x) return(x$trueDeg))))
  
  lens = sapply(dataset, function(x) length(x$kSamp))
  lens = data.frame(lens = lens, trial = 1:length(lens), maxL = rep(1:4, each = trialsPer))
  lens$cf = cumsum(lens$lens)
  
  # (ugly) add max steps
  df$maxL = NA
  df$maxL[1:lens$cf[max(which(lens$maxL == 1))]] = 1
  for (i in 2:4) {
    df$maxL[(lens$cf[min(which(lens$maxL == i))-1]+1):lens$cf[max(which(lens$maxL == i))]] = i
  }
  
  # (ugly) add trial number
  df$trial = NA
  df$trial[1:lens$cf[1]] = 1
  for (i in 2:length(dataset)) {
    df$trial[(lens$cf[i-1]+1):lens$cf[i]] = i
  }
  
  df$kt = ifelse(df$kTrue > targ-1, targ-1, df$kTrue)
  
  df$err = with(df, (kt - kSamp)/(kt - 1))
  rm(lens)
  
  qerr = df %>%
    group_by(kTrue, maxL) %>%
    summarise(trials = length(err), median = median(err, na.rm = TRUE), 
              q25 = quantile(err, .25, na.rm = TRUE), q75 = quantile(err, .75, na.rm = TRUE))
  
  qerr = qerr %>%
    group_by(maxL) %>%
    mutate(div = calcDiv(trials, deg.gr))
  
  return(list(indivs = df, errDiv = qerr))
  
}


haloThreshold = function(gr, thresh, trialsPer, targ) {
  
  # initialize
  lip = list()
  grel = as_edgelist(gr)
  deg.gr = degree(gr)
  deg.df = data.frame(name = as.numeric(names(deg.gr)), deg = deg.gr)
  
  for (i in thresh) {
    
    for (trial in 1:trialsPer) {
      
      e.l = matrix(nrow = 0, ncol = 2)
      
      # set initial seed and begin snowball
      sn = sample(V(gr)$name, 1)
      e.l = stepOut(grel, e.l, sn)
      e.l$halo = e.l$halo[e.l$halo != sn]
      
      while (length(unique(array(e.l$g))) < targ) {
        # if the halo exists, and its mean degree is less than the threshold, 
        # then step out
        if (length(e.l$halo) & mean(deg.gr[as.character(e.l$halo)]) < i) {
          e.l = stepOut(grel, e.l$g, e.l$halo)
        # else, begin a new seed by sampling from the set of partially/fully unsampled nodes
        } else {
          incompl = merge(deg.df, data.frame(table(e.l$g)), by.x = 'name', by.y = 'Var1', all = TRUE)
          sn = as.numeric(sample(incompl$name[incompl$deg != incompl$Freq | is.na(incompl$Freq)], 1))
          e.l = stepOut(grel, e.l$g, sn)
          e.l$halo = e.l$halo[e.l$halo != sn]
        }
      }
      
      # if target exceeded, remove nodes until sample is appropriate size
      toRm = sample(e.l$halo, length(unique(array(e.l$g))) - targ)
      e.l$g = e.l$g[!(e.l$g[,1] %in% toRm | e.l$g[,2] %in% toRm),]
      
      deg = table(e.l$g)
      
      thIdx = which(thresh == i)
      lip[[(thIdx-1)*trialsPer + trial]] = list(kSamp = as.numeric(deg), 
                                                trueDeg = as.numeric(deg.gr[names(deg)]))
      
    }
  }
  
  return(lip)
  
}


haloEdgeSampling = function(gr, thresh, trialsPer, targ) {

  ### This function uses the proportion of realized within-halo edges to decide
  ### whether to step out or begin another snowball
  ### A within-halo edge is one with both incident nodes in the halo
  ### the proportion of realized edges is the number of edges divided by the
  ### size of the halo (in nodes) choose two
    
  # initialize
  liv = list()
  grel = as_edgelist(gr)
  deg.gr = degree(gr)
  deg.df = data.frame(name = as.numeric(names(deg.gr)), deg = deg.gr)
  
  
  for (i in thresh) {
    
    for (trial in 1:trialsPer) {
      
      e.l = matrix(nrow = 0, ncol = 2)
      
      # randomly select seed node
      sn = sample(V(gr)$name, 1)
      # begin first snowball by taking step out from seed
      e.l = stepOut(grel, e.l, sn)
      # remove the seed from the halo
      e.l$halo = e.l$halo[e.l$halo != sn]
      
      while (length(unique(array(e.l$g))) < targ) {
        
        # calculate number of within-halo edges
        winedge = sum((apply(grel, 1, function(x) sum(x %in% e.l$halo) == 2)))
        # output proportion of possible within-halo edges realized
        p_edge = winedge / choose(length(e.l$halo), 2)
        
        if (!is.nan(p_edge) & p_edge < i) {
          e.l = stepOut(grel, e.l$g, e.l$halo)
        } else {
          incompl = merge(deg.df, data.frame(table(e.l$g)),
                          by.x = 'name', by.y = 'Var1', all = TRUE)
          # randomly choose a seed node from the set of incompletely-sampled nodes
          sn = as.numeric(sample(incompl$name[incompl$deg != incompl$Freq |
                                                is.na(incompl$Freq)], 1))
          e.l = stepOut(grel, e.l$g, sn)
          e.l$halo = e.l$halo[e.l$halo != sn]
        }
      }
      
      # removing nodes to get us to the target size
      toRm = sample(e.l$halo, length(unique(array(e.l$g))) - targ)
      e.l$g = e.l$g[!(e.l$g[,1] %in% toRm | e.l$g[,2] %in% toRm),]
      
      # degree distribution of the sampled graph
      deg = table(e.l$g)
      
      thIdx = which(thresh == i)
      liv[[trialsPer*(thIdx-1) + trial]] = list(kSamp = deg, trueDeg = deg.gr[names(deg)])
      
    }
    
  }
  
  return(liv)
  
}