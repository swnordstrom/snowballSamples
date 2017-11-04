configGraph = function(d, allVerts = FALSE) {

  nStub = sum(d * 1:length(d))
	if (nStub %% 2) stop("Need an even number of stubs")	
	
	el = matrix(nrow = 0, ncol = 2)
	csd = c(0, cumsum(d))
	
	for (i in 1:(length(d)-1)) for (j in (i+1):length(d)) {
	  if (i == j) edges = runif(choose(d[i], 2)) < i^2 / nStub
	  else	      edges = runif(d[i] * d[j]) < i*j / nStub
	  localNodes = expand.grid(csd[i] + 1:d[i], csd[j] + 1:d[j])
	  el = rbind(el, localNodes[localNodes[,1] < localNodes[,2],][edges,])
	}
	
	g.out = graph_from_edgelist(as.matrix(el), FALSE)
	
	if (!allVerts) {
	  deg = degree(g.out)
	  g.out = delete.vertices(g.out, !deg)
	  g.out = set_vertex_attr(g.out, "name", value = which(deg>0))
	}
	
  return(g.out)

}