##############################################
#   FUNCTIONS TO TRIM EMPIRICAL NETWORKS AND COMPUTE METRICS
##############################################


# ------------------------------------------------------------
# Generate a random bipartite binary matrix with no isolated rows
# ------------------------------------------------------------
# Creates a rows × cols binary incidence matrix (e.g. victims × exploiters).
# Each row is forced to contain at least one interaction to avoid isolated
# victim species at initialization. Columns are not explicitly constrained.
#
# Input:
#   rows  - number of victim species
#   cols  - number of exploiter species
#
# Output:
#   Binary matrix with dimensions rows × cols

generate_binary_matrix <- function(rows, cols) {
  mat <- matrix(0, nrow = rows, ncol = cols)
  
  # Generating random binary matrix
  for (i in 1:rows) {
    while (all(mat[i,] == 0)) {
      mat[i,] <- sample(0:1, cols, replace = TRUE)
    }
  }
  
  return(mat)
}


# ------------------------------------------------------------
# Construct the full square adjacency matrix ("F matrix")
# ------------------------------------------------------------
# Converts a bipartite incidence matrix into a square adjacency matrix
# of size (n_v + n_e) × (n_v + n_e), with interactions in the off-diagonal
# blocks and zeros on the diagonal blocks:
#
#     [ 0   B ]
#     [ B'  0 ]
#
# This representation is useful for computing indirect effects and
# graph-based metrics on bipartite networks.
#
# Input:
#   binmat - binary incidence matrix (victims × exploiters)
#
# Output:
#   Square adjacency matrix

create_F <- function(binmat) {
  mat_F<-rbind (
    cbind( matrix(0, nrow=nrow(binmat), ncol=nrow(binmat)),binmat),
    cbind( t(binmat), matrix(0, nrow=ncol(binmat),ncol=ncol(binmat)))) 
  return(mat_F)
}





# ------------------------------------------------------------
# Remove isolated species from a bipartite matrix
# ------------------------------------------------------------
# Identifies and removes rows (victims) and/or columns (exploiters)
# that have zero interactions. This ensures that downstream network
# metrics are computed only on species participating in the network.
#
# Input:
#   network_matrix - bipartite incidence matrix
#
# Output:
#   Pruned incidence matrix with no all-zero rows or columns

remove_isolated_species <- function(network_matrix) {
  
  # Convert dataframe to matrix
  network_matrix <- as.matrix(network_matrix)
  
  # Identify rows and columns with all zeros
  zero_rows <- which(rowSums(network_matrix) == 0)
  zero_cols <- which(colSums(network_matrix) == 0)
  
  if(length(zero_cols) == 0 & length(zero_rows) > 0) {
    
    # Remove rows and columns with all zeros
    pruned_matrix <- network_matrix[-zero_rows,]
    
  }else if (length(zero_cols) > 0 & length(zero_rows) == 0){
    
    pruned_matrix <- network_matrix[, -zero_cols]
    
  }else if (length(zero_cols) > 0 & length(zero_rows) > 0){
    
    # Remove rows and columns with all zeros
    pruned_matrix <- network_matrix[-zero_rows, -zero_cols]
    
  }else if (length(zero_cols) == 0 & length(zero_rows) == 0){
    
    pruned_matrix = network_matrix
  }
  
  
  
  return(pruned_matrix)
}



# ------------------------------------------------------------
# Estimate power-law exponent of a degree distribution
# ------------------------------------------------------------
# Fits a discrete power-law model to strictly positive degree values
# using the poweRlaw package and returns the scaling exponent alpha.
#
# Input:
#   vec_degrees - numeric vector of node degrees
#
# Output:
#   Estimated power-law exponent alpha


compute_alpha_pl <- function(vec_degrees) {
  
  fit <- displ$new(vec_degrees[vec_degrees > 0])
  alpha <- estimate_pars(fit)$pars
  
  return(alpha)
  
}


# ------------------------------------------------------------
# Compute Shannon entropy of a degree distribution
# ------------------------------------------------------------
# Treats degree frequencies as a probability distribution and computes
# Shannon entropy (natural logarithm). Higher values indicate more even
# degree distributions.
#
# Input:
#   degrees - numeric vector of degrees
#
# Output:
#   Shannon entropy of the degree distribution


compute_entropyd <- function(degrees) {
  probabilities <- table(degrees) / length(degrees)
  entropy <- -sum(probabilities * log(probabilities, base = exp(1)))
  return(entropy)
}



# ------------------------------------------------------------
# Compute SVD-based entropy of a bipartite matrix
# ------------------------------------------------------------
# Computes an entropy measure based on the singular value spectrum
# of the bipartite incidence matrix. Singular values are normalized
# and used to quantify structural complexity.
#
# Input:
#   binary_matrix - bipartite incidence matrix
#
# Output:
#   Normalized SVD entropy


compute_SVD_entropy <- function(binary_matrix) {
  svd_result <- svd(binary_matrix)  # Compute SVD
  
  # Extract singular values
  singular_values <- svd_result$d
  
  # Filter out zero singular values
  non_zero_singular_values <- singular_values[singular_values > 0]
  
  # Normalize singular values to have values between 0 and 1
  normalized_singular_values <- non_zero_singular_values / max(non_zero_singular_values)
  
  # Compute k (rank of the matrix)
  k <- length(non_zero_singular_values)
  
  # Compute SVD entropy
  entropy <- - (1 / log(k)) * sum(normalized_singular_values * log(normalized_singular_values))
  
  return(entropy)
}


# ------------------------------------------------------------
# Compute species-level centrality metrics in bipartite networks
# ------------------------------------------------------------
# Converts a bipartite incidence matrix into an undirected igraph object
# and computes degree and betweenness centrality for all species, as well
# as separately for victims and exploiters.
#
# Input:
#   binmatrix - bipartite incidence matrix
#   n_v       - number of victim species (rows)
#   n_e       - number of exploiter species (columns)
#
# Output:
#   List containing degree and betweenness vectors for all species,
#   victims only, and exploiters only


compute_centrality <- function(binmatrix, n_v, n_e) {
  
  nspp = n_v + n_e
  
  graph <- graph_from_incidence_matrix(binmatrix, directed = FALSE)
  
  betw <- igraph::betweenness(graph)
  betw_v <- igraph::betweenness(graph)[1:n_v]
  betw_e <- igraph::betweenness(graph)[(n_v+1):nspp]
  
  degree <- igraph::degree(graph)
  degree_v <- degree[1:n_v]
  degree_e <- degree[(n_v+1):nspp]
  
  list <- list("betw" = betw, "betw_v" = betw_v, "betw_e" = betw_e, "degree" = degree, "degree_v" = degree_v, "degree_e" = degree_e)
  
  return(list)
  
}




# ------------------------------------------------------------
# Plot relationships between network structure and net effects
# ------------------------------------------------------------
# Produces a multi-panel figure linking structural network metrics
# (modularity, connectance, richness, degree heterogeneity, complexity,
# entropy) to a user-specified response variable (e.g. mean or entropy
# of biotic effects).
#
# Input:
#   variable      - response variable to map to color or y-axis
#   lab_variable  - label for the color scale
#   dataset       - dataframe containing network metrics
#
# Output:
#   ggarrange object combining multiple structural plots


plot_str_neteffects <- function(variable, lab_variable, dataset) {
  
  
  # Plot modularity against nestedness with T_mean and T_nf1_mean
  p_mod_nest_meanT <-   ggplot(dataset, aes(x = modularity, y = nestedness)) +  
    geom_point(aes(color = variable, size = nspp), alpha = 0.8) +  
    labs(x = expression(italic(Q)), y = "NODF") +
    scale_color_viridis_c(name = lab_variable) +
    my_theme
  
  
  # Plot connectance against NSPP with T_mean and T_nf1_mean
  plot_conn_meanT <- ggplot(dataset, aes(x = connectance, y = nspp)) +
    geom_point(aes(color = variable, size = T_nf1_mean), alpha = 0.8, size = 3) +  # Add T_mean and T_nf1_mean to color and shape aesthetics
    labs(x = "C", y = "S") +
    scale_color_viridis_c(name = lab_variable) +
    my_theme
  
  
  # Plot nexpl and nvict against nvictim with T_mean and T_nf1_mean
  plot_nexpl_nvictim_T_mean <- ggplot(dataset, aes(x = log(nexpl), y = log(nvictim))) +
    geom_point(aes(color = variable), alpha = 0.8, size = 3) +  # Add T_mean and T_nf1_mean to color and shape aesthetics
    labs(x = "log(Se)", y = "log(Sv)") +
    scale_color_viridis_c(name = lab_variable) +
    my_theme
  
  # Plot distribution of alpha values with T_mean and T_nf1_mean as vertical lines
  plot_alpha_tmean <- ggplot(dataset, aes(x = alpha, y = variable)) +
    geom_point( alpha = 0.5, size = 3) +
    labs(x = expression(alpha)) +
    my_theme
  
  # Plot distribution of complexity index with T_mean and T_nf1_mean as vertical lines
  plot_complexity_tmean <- ggplot(dataset, aes(x = complexity_index, y = variable)) +
    geom_point( alpha = 0.5, size = 3) +
    labs(x = expression(sqrt(C * S))) +
    my_theme
  
  
  
  # Plot distribution of entropy_d with T_mean and T_nf1_mean as vertical lines
  plot_entropy_tmean <- ggplot(dataset, aes(x = entropy, y = variable)) +
    geom_point( alpha = 0.5, size = 3) +
    labs(x = "entropy") +
    my_theme
  
  
  arranged_plots_str <- ggarrange(
    
    p_mod_nest_meanT,
    plot_conn_meanT,
    plot_nexpl_nvictim_T_mean,
    plot_alpha_tmean,
    plot_complexity_tmean,
    plot_entropy_tmean,
    
    
    ncol = 2,
    nrow = 3,
    
    common.legend = TRUE,
    
    labels = c("A", "B", "C", "D", "E", "F", "G")
    
    
  )
  
  return(arranged_plots_str)
  
}
