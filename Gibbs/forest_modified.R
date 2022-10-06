rgumbel<- function(n){
  -log(-log(runif(n)))
}

findMST<- function(logS){
  G<- graph_from_adjacency_matrix(-logS, mode = "undirected",weighted = TRUE,diag=FALSE)
  G_mst<- mst(graph = G)
  A_T<- as.matrix(get.adjacency(G_mst))
  A_T
}

drawT_approx<- function(logS){
  n = nrow(logS)-1 #adjusts for ease of next parts
  gumbelMat<- matrix(0,n+1,n+1)
  gumbelMat[lower.tri(gumbelMat,diag = FALSE)]<- rgumbel((n+1)/2*(n))
  gumbelMat<- gumbelMat+ t(gumbelMat)
  A_T<- findMST(logS+gumbelMat)
  A_T
}





gumbelMax<- function(logA){
  which.max( logA+ rgumbel(length(logA)))
}


drawT_exact<- function(logS){
  
  A_T<- matrix(0,n+1,n+1)
  
  InTree<- list()
  Next <- list()
  
  for (i in 1:(n+1)){
    InTree[[i]]<- FALSE
  }
  
  r<- n+1
  InTree[[r]]<- TRUE
  Next[[r]]<- n+1
  
  for (i in (n+1):1){
    u = i
    while(!InTree[[u]]){
      Next[[u]]<- gumbelMax(logS[u,])
      u <- Next[[u]]
    }
    u = i
    while(!InTree[[u]]){
      InTree[[u]]= TRUE
      u <- Next[[u]]
    }
  }
  
  for (u in 1:n){
    A_T[u, Next[[u]]]=1
    A_T[Next[[u]],u]=1
  } 
  A_T
}

