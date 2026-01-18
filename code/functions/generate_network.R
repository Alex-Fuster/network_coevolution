#### FUNCTION TO GENERATE NETWORKS FROM DEGREE DISTRIBUTIONS ######


generate_distrib <- function(connectance, n1, n2, distribution_type, nestedness, modularity) {
  # Random Network (Poisson Distribution)
  if (distribution_type == "random") {
    ones <- ceiling(connectance * n1 * n2)
    p1 <- (1:n1) / (n1 + 1)
    p2 <- (1:n2) / (n2 + 1)
    lambda1 <- connectance * n2
    lambda2 <- connectance * n1
    ED1 <- 1 + qpois(p = p1, lambda = lambda1 - 1)
    ED2 <- 1 + qpois(p = p2, lambda = lambda2 - 1)
    if (sum(ED1) < sum(ED2)) {
      ED1[1] <- ED1[1] + 1
    }
    if (sum(ED2) < sum(ED1)) {
      ED2[1] <- ED2[1] + 1
    }
  }
  # Scale-Free Network (Pareto Distribution)
  else if (distribution_type == "scale-free") {
    # Implement the generation of a scale-free network (Pareto distribution)
    
    graph <- igraph::barabasi.game(n = n1 + n2, m = 1)
    # Extract degrees from the graph and split them into two vectors ED1 and ED2
    ED1 <- igraph::degree(graph)[1:n1]
    ED2 <- igraph::degree(graph)[(n1 + 1):(n1 + n2)]
  }
  
  # Nested Network
  else if (distribution_type == "nested") {
    nested_matrix <- matrix(0, nrow = n1, ncol = n2)
    for (i in 1:n1) {
      nested_matrix[i, sample(1:n2, nestedness)] <- 1
    }
    ED1 <- rowSums(nested_matrix)
    ED2 <- colSums(nested_matrix)
  }
  else {
    stop("Invalid distribution_type. Supported types: random, scale-free, nested")
  }
  
  
  # Introduce modularity
  if (modularity > 0) {
    # Adjust connections within modules
    ED1[1:floor(n1/2)] <- ED1[1:floor(n1/2)] + modularity
    ED2[1:floor(n2/2)] <- ED2[1:floor(n2/2)] + modularity
  }
  
  # Return the result as a list
  return(list(ED1, ED2))
}


# Example usage:
#set.seed(123)
#result <- generate_distrib(connectance = 0.2, n1 = 15, n2 = 25, distribution_type = "scale-free", nestedness = 3, modularity = 8)

#mat <- generate.BEDD(result[[1]], result[[2]])



randomize.BEDD<-function(bipmat){ #bipmat is an existing bipartite matrix, the function randomizes it with the aim of keeping the expected degrees
  s<-sum(bipmat)
  m1<-apply(bipmat,1,sum)
  m2<-apply(bipmat,2,sum)
  if(max(m1)*max(m2) > s) {
    print("Warning, some probabilities are higher than 1")
  }
  apply(m1%o%m2/s,c(1,2),function(x) rbinom(1,size = 1,prob = min(x,1)))
}

generate.BEDD<-function(ED1,ED2){ #ED1 and ED2 are the vectors of expected degrees
  s1<-sum(ED1)
  s2<-sum(ED2)
  if(s1 != s2) {
    print("Warning, the two sequences have different sums")
  }
  s<-max(s1,s2)
  if(max(ED1)*max(ED2) > s) {
    print("Warning, some probabilities are higher than 1")
  }
  apply(ED1%o%ED2/s,c(1,2),function(x) rbinom(1,size = 1,prob = min(x,1)))
}

generate.EDs.poisson<-function(connectance, n1, n2){ #this generates two "deterministic" sequences of degrees following a Poisson distribution, using the quantiles
  ones<-ceiling(connectance*n1*n2)
  p1<-(1:n1)/(n1+1)
  p2<-(1:n2)/(n2+1)
  lambda1<-connectance*n2
  lambda2<-connectance*n1
  ED1<- qpois(p = p1, lambda = lambda1)
  ED2<- qpois(p = p2, lambda = lambda2)
  if(sum(ED1) < sum(ED2)) {
    ED1[1]<-ED1[1] + 1
  }
  if(sum(ED2) < sum(ED1)) {
    ED2[1]<-ED2[1] + 1
  }
  list(ED1,ED2)
}

generate.EDs.truncatedpoisson<-function(connectance, n1, n2){ #as above, but takes care that none of the nodes has a degree of zero
  ones<-ceiling(connectance*n1*n2)
  p1<-(1:n1)/(n1+1)
  p2<-(1:n2)/(n2+1)
  lambda1<-connectance*n2
  lambda2<-connectance*n1
  ED1<- 1+qpois(p = p1, lambda = lambda1-1)
  ED2<- 1+qpois(p = p2, lambda = lambda2-1)
  if(sum(ED1) < sum(ED2)) {
    ED1[1]<-ED1[1] + 1
  }
  if(sum(ED2) < sum(ED1)) {
    ED2[1]<-ED2[1] + 1
  }
  list(ED1,ED2)
}
