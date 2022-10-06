########################
# GIBBS SAMPLING FOR LATENT COORDINATE SPANNING TREE MODEL
########################

path = 'C:/Users/josmi/OneDrive/Documents/Graduate School/Research/7979 Fall 2022'
#replace with your own path to run

source(sprintf('%s/stuff/Bayesian_forest_clustering/forest_modified.R',path))

########################
# Utility functions
########################

normalize_Z = function(Z){
  if( sum( colSums(Z**2) )  != as.numeric(ncol(Z)) ){ #logical is not evaluating as intended
    lengths = diag(t(Z)%*%Z) |> sqrt()
    for(i in 1:ncol(Z)){ #standardize columns of Z
      Z[,i] = Z[,i]/lengths[i]
    }
    warning('Columns of Z will be normalized.')
  }
  return(Z)
} #not currently implemented

compute_Laplace = function(Z,A){ #maybe 1 sec for n=959
  n = nrow(Z)
  
  Adj = matrix(0, nrow=n, ncol=n)
  W = Z%*%t(Z)
  
  vcosine = function(v1,v2){
    crossprod(v1,v2)/( norm(v1,'2')*norm(v2,'2'))
  }
  
  #d=1 adaptation
  
  P = vcosine(t(Z),t(Z))
    Adj[lower.tri(I(n))] = dbinom(A[lower.tri(I(n))], 1, P[lower.tri(I(n))])*W[lower.tri(I(n))]
      # apply(A[lower.tri(I(n))], 1, dbinom, P[lower.tri(I(n))])

  Adj = Adj + t(Adj)
  
  D = diag(rowSums(Adj)-diag(Adj))
  return(D-Adj)
}

########################
# Load a graph and generate initial points for the sampler
########################

{

require(igraph)
require(matrixNormal)
g1 = read_graph(sprintf('%s/Data/talairach_subject1.graphml', path), format='graphml')

A = get.adjacency(g1, type='both', attr = NULL)
n = nrow(A) # number of nodes in graph
set.seed(1)
t = mst(graph=g1)
U = rmatnorm(s=1, M=matrix(0,n,2), U=I(n), V=I(2))
d = 1 #dimension of latent coordinate vectors
# Z = matrix(runif(n*d), nrow=n, ncol=d) # elements too small in size, multiply to 0
# Z = matrix(rbeta(n*d,shape1=14,shape2=1),ncol=d,nrow=n) #push Z values towards 1

z_shape = 2
z_scale = 2 #tuning parameters?
Z = matrix(rgamma(n*d, shape=z_shape, scale=z_scale),ncol=d,nrow=n)

# p = matrix(runif(n*n), ncol=n, nrow=n) #similar problem as z
# p = matrix(rbeta(n*n, shape1=9,shape2=1), ncol=n, nrow=n) #push probabilities towards 1

}

########################
# Gibbs sampler
########################

gen_U = function(Z,A){ #should this depend on A? I think it's fine, since we just need full conditionals
  # normalize_Z(Z)
  
  L = compute_Laplace(Z,A)
  n = nrow(Z)
  # need to make sure Laplace is computed beforehand
  # will also want to try to speed compute_Laplace() up dramatically
  return( rmatnorm(M=matrix(0, nrow=n, ncol=2), U=L + J(n)/(n^2), V=I(2)))
}

gen_ST = function(Z,verbose=F){
  if(verbose) print('Starting Tree Generation:')
  logS = log(Z%*%t(Z))
  drawT_approx(logS=logS)
}

gen_Z = function(t,U,A,verbose=F){ #currently only works for d=1 and gamma prior on Z
  if(verbose) print('Starting Z generation:')
  n = nrow(A)
  U_scale = function(U_i){
    row_diff = function(U_j, U_i){
      norm(U_i-U_j,'2')^2
    }
    
    sum_comp = apply(U, 1, row_diff, U_i) |> sum()
    return((1/z_scale + .5*sum_comp)^(-1))
  }
  scale = apply(U, 1, U_scale) #this is the biggest time bottleneck
  
  # result = rep(NA,n)
  
  result = apply(matrix(1,nrow=n), 1, rgamma, shape=z_shape+1, scale=scale)
  
  # for(i in 1:n){
  #   # result[i] = rgamma(n=1, shape=z_shape+1, scale=U_scale(U[i,]))
  #   if(verbose & i%%50==0) print(i)
  # }
  result = matrix(result,nrow=n)
  result
}

sample_lcst = function(A,t,U,Z,verbose=F, max_iter=500){
  i = 1 #iteration counter
  converged=0
  # t_tilde = #how to store several igraph objects?
  # U_tilde = array(NA, dim=c(dim(U)[1], dim(U)[2], max_iter+1))
  Z_tilde = array(Z, dim=c(dim(Z)[1], dim(Z)[2], max_iter+1))
  
  while(!converged){
    if(verbose) print(sprintf("Gibbs sampler iteration %d",i))
    U = gen_U(Z_tilde[,,i],A)
    t = gen_ST(Z_tilde[,,i],verbose=verbose)
    Z_tilde[,,i+1] = gen_Z(t,U,A,verbose=verbose)
    
  if(i >= max_iter){
    converged=1
  }
    i=i+1
  }
  Z_tilde
}



########################
# Density functions
########################

lkl = function(A, t, Z, w=1, p){ #currently adapted only for 1 tree
  if(nrow(A) != nrow(p) | ncol(A) != ncol(p) ){
    stop("A and p must be matrices of the same dimensions.")
  }
  if(sum(w)!=1){
    w = w/sum(w)
    warning('Rescaling weight vector to sum to 1 across entries')
  }
  
  edges_t = get.adjacency(t, type='both') |> as.matrix()
  edge_ind = which(edges_t==1, arr.ind = T)
  A_edge = A[edge_ind] #all ones, trivially
  p_edge = p[edge_ind]
  

  # return( sum( w*prod( dbinom(x=A,size=1,prob=p)) ) )
  return( sum( w*prod(p_edge) ) )
} #not currently implemented

pi_tree = function(t, Z, A){
  # ZZt = Z%*%t(Z)
  
 # Z = normalize_Z(Z)
 
 L = compute_Laplace(Z,A,verbose=T)
 
 n = nrow(Z)
 
 edges_t = (Z%*%t(Z))*( get.adjacency(t, type='both') |> as.matrix())
 weight_t = edges_t[edges_t!=0] |> as.vector()
 eps = 0.1887
 # prod(weight_t+eps)
 prior = prod(weight_t+eps)/determinant(L+J(n)/(n^2), logarithm = F)$modulus
 
 return( list(prior=prior, Laplace=L) )
 # determinant is extremely large
  
} #not currently implemented

########################
# Ancillary Code
########################

# Example graph

# e = matrix(0, nrow=100, ncol=2)
# for(i in 1:nrow(e)){
#   e[i,] = c(ceiling(i/10),ifelse(i%%10>0, i%%10, 10))
# }
# g = graph(t(e), directed=F)
# 

