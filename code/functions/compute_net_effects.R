##############################################
#   FUNCTIONS TO COMPUTE NET EFFECTS FROM BINARY MATRIX
##############################################


# Function to compute net effects that species receive and give in a network


compute_total_effects <- function(int_matrix) {
  
  f_mat <- create_F(int_matrix)
  
  # generating a matrix that all rows sum to one. In this matrix the element describes the proportional effect of a species in the column on the species in the row
  
  P = f_mat
  
  
  for (i in 1:nrow(f_mat)) {
    P[i,] <- P[i,] / sum(P[i,])
  }
  
  
  #weighting each row for the total contribution of interactions for a given process (e.g., the amount of selection imposed by mutualistic partners on the mean trait of a given species
  
  m=0.6 # m is the strength and needs to be between zero and one
  
  Q=m*P 
  
  #Note: you can change m for a diagonal matrix if you want m to be species specific
  
  #computing the matrix of total effects (direct and indirect)
  I <- diag(nrow(f_mat)) #identity matrix
  T <- solve(I - Q)  #compute the inverse of I-Q
  
  return(T)
  
}




#Function to calculate normalized Shannon Index for a row, which we call net flow 1 (i.e. total effects that each species receives):
  
# log2 or log?

compute_shannon <- function(vec) {
  eps <- 1e-10  # Small constant to avoid taking the logarithm of zero
  vec_proportions <- vec / sum(vec)
  vec_proportions[vec_proportions == 0] <- eps  # Replace zeros with eps
  shannon_index <- -sum(vec_proportions * log(vec_proportions))
  normalized_shannon <- shannon_index / log(length(vec_proportions))
  return(normalized_shannon)
}



compute_shannon_nonorm <- function(vec) {
  eps <- 1e-10  # Small constant to avoid log(0)
  vec_proportions <- vec / (sum(vec) + eps)  # Avoid division by zero
  vec_proportions[vec_proportions < eps] <- eps  # Avoid log(0)
  shannon_index <- -sum(vec_proportions * log(vec_proportions))
  normalized_shannon <- shannon_index / log(length(vec_proportions) + eps)
  return(shannon_index)
}










