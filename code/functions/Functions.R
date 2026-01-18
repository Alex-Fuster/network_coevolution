##############################################
#   FUNCTIONS FOR THE SIMULATION - EVOLUTION OF TRAITS IN FIXED VICTIM-EXPLOITER NETWORK STRUCTURE
##############################################

# ============================================================
# Biotic selection values
# ============================================================

# ------------------------------------------------------------
# Update biotic-selection strength from percentage differences
# ------------------------------------------------------------
# Takes an original selection strength and a vector of % changes and returns a
# new vector of selection strengths:
#   new_bs = original_bs + (perc_difference/100) * original_bs
#
# Inputs:
#   original_bs          - scalar baseline biotic selection value
#   vec_perc_difference  - vector of percentage deviations (e.g., -10, 0, +20)
#
# Output:
#   Numeric vector of updated biotic selection values


compute_bs_from_perc.difference <- function(original_bs, vec_perc_difference) {
  
  
  vec_new_bs <- numeric(length(vec_perc_difference))
  
  for (i in 1:length(vec_perc_difference)) {
    
    vec_new_bs[i] = vec_perc_difference[i]/100*original_bs+original_bs
    
  }
  
  
  return(vec_new_bs)
  
}




# ============================================================
# Network generation (fixed size / connectance; optional trait matching)
# ============================================================

# ------------------------------------------------------------
# Generate a bipartite network with a target connectance (probabilistic sampling)
# ------------------------------------------------------------
# Uses a trait-matching probability matrix to sample a fixed number of realized
# interactions (n_ones = round(n_host * n_exploiter * conn)). Returns the binary
# incidence matrix L and the realized connectance C.
#
# Inputs:
#   n_exploiter - number of exploiters (columns)
#   n_host      - number of victims/hosts (rows)
#   conn        - target connectance (0-1)
#
# Output:
#   List with:
#     [[1]] L  - binary incidence matrix (hosts × exploiters)
#     [[2]] C  - realized connectance


generate_network_given_conn <- function(n_exploiter, n_host, conn) {
  
  
  prob_matrix<-generate_prob_matrix(n.host = n_host, n.exploiter = n_exploiter)
  
  index<-which(prob_matrix>0, arr.ind = TRUE)
  
  n_ones <- round(n_exploiter*n_host*conn)
  
  realized.int<-sample(1:(n_exploiter*n_host), size = n_ones, prob = prob_matrix[index])
  
  L <- matrix(0, nr = n_host, nc = n_exploiter)
  
  L[index[realized.int, ]] <-1
  
  
  C<-sum(L)/(n_exploiter*n_host)
  
  result<-list(L, C)
  
  return(result)
  
}


# ------------------------------------------------------------
# Generate a trait-matching probability matrix (used for network sampling)
# ------------------------------------------------------------
# Creates random trait values for all species, computes pairwise trait differences,
# and converts host-exploiter differences into interaction probabilities using a
# Gaussian matching function:
#   p_ij = exp( - ( (x_i - x_j) / scale )^2 )
#
# Inputs:
#   n.host      - number of hosts (victims)
#   n.exploiter - number of exploiters
#   mean.trait  - mean of the trait distribution used to initialize traits
#   sd.trait    - SD of the trait distribution used to initialize traits
#   scale       - matching width (smaller = sharper preference for similarity)
#
# Output:
#   Matrix of size n.host × n.exploiter with interaction probabilities



generate_prob_matrix<- function(n.host, n.exploiter, mean.trait=0, sd.trait=.1, scale=.1){
  
  nsp = n.host + n.exploiter
  
  trait<-rnorm(nsp, mean=mean.trait, sd = sd.trait)
  
  dif.traits <- outer(trait, trait, FUN= "-") 
  
  # Create trait matching matrix
  
  matching_matrix<-dif.traits[(n.exploiter+1):nsp, 1:n.exploiter]
  
  
  # compute probabilities based on trait matching
  prob_matching<-exp(-(matching_matrix/scale)^2)
  
  return (prob_matching)
  
}



# ------------------------------------------------------------
# Generate a network at each simulation step from current traits (connectance by quantile)
# ------------------------------------------------------------
# Builds a probability matrix from the current trait vector and realizes links by
# taking the top fraction 'conn' of probabilities (thresholded by quantile).
# Ensures that every exploiter has at least one interaction by forcing one link
# to its best-matching host if needed.
#
# Inputs:
#   n_exploiter - number of exploiters
#   n_host      - number of hosts
#   conn        - target connectance (fraction of potential links realized)
#   trait_sim   - numeric trait vector of length n_host + n_exploiter
#
# Output:
#   Binary incidence matrix L (hosts × exploiters)



generate_network_given_conn_each_sim <- function(n_exploiter, n_host, conn, trait_sim) {

 
  
  prob_matrix<-generate_prob_matrix_each_sim(n.host = n_host, n.exploiter = n_exploiter,
                                             trait_vec = trait_sim)
  
  
  #n_ones <- round(n_exploiter*n_host*conn)
  
  q_prob_int<- quantile(prob_matrix, prob = 1 - conn)
  
  index<-which(prob_matrix > q_prob_int, arr.ind = TRUE)
  
  L <- matrix(0, nr = n_host, nc = n_exploiter)
  
  L[index] <-1
  
  
  expl_no_int<-colSums(L)
  
  ## Check what columns and what rows have no ones
  
  
  
  vec_e_no_int<-which(expl_no_int == 0)
  
  for (i in vec_e_no_int) {
    
    
    v_toeat <- which.max(prob_matrix[,i])
    
    L[v_toeat,i] <- 1
    
  }
  
  return(L)
  
  
}


# ------------------------------------------------------------
# Trait-driven probability matrix for a given simulation step
# ------------------------------------------------------------
# Same trait-matching rule as above, but using a provided trait vector.
#
# Inputs:
#   n.host      - number of hosts
#   n.exploiter - number of exploiters
#   trait_vec   - numeric trait vector (length n.host + n.exploiter)
#   scale       - matching width
#
# Output:
#   Matrix of size n.host × n.exploiter of matching probabilities


generate_prob_matrix_each_sim <- function(n.host, n.exploiter, trait_vec, scale = 0.1) {
  
  nsp = n.host + n.exploiter
  
  dif.traits <- outer(trait_vec, trait_vec, FUN= "-")
  
  # Create trait matching matrix
  
  matching_matrix<-dif.traits[(n.exploiter+1):nsp, 1:n.exploiter]
  
  # compute probabilities based on trait matching
  prob_matching<-exp(-(matching_matrix/scale)^2)
  
  return (prob_matching)
  
}



# ------------------------------------------------------------
# Generate a trait-driven network ensuring all species have ???1 interaction
# ------------------------------------------------------------
# Same as generate_network_given_conn_each_sim(), but additionally enforces that
# every host also has at least one interaction (i.e., no isolated hosts).
#
# Inputs / Output: same as generate_network_given_conn_each_sim()

generate_network_given_conn_each_sim_allsp.int <- function(n_exploiter, n_host, conn, trait_sim) {
  
  
  prob_matrix<-generate_prob_matrix_each_sim(n.host = n_host, n.exploiter = n_exploiter,
                                             trait_vec = trait_sim)
  
  
  #n_ones <- round(n_exploiter*n_host*conn)
  
  q_prob_int<- quantile(prob_matrix, prob = 1 - conn)
  
  index<-which(prob_matrix > q_prob_int, arr.ind = TRUE)
  
  L <- matrix(0, nr = n_host, nc = n_exploiter)
  
  L[index] <-1
  
  
  expl_no_int<-colSums(L)
  

  
  vec_e_no_int<-which(expl_no_int == 0)
  
  for (i in vec_e_no_int) {
    
    
    v_toeat <- which.max(prob_matrix[,i])
    
    L[v_toeat,i] <- 1
    
  }
  
  
  host_no_int<-rowSums(L)
  


  
  vec_h_no_int<-which(host_no_int == 0)
  
  for (i in vec_h_no_int) {
    
    
    e_predating <- which.max(prob_matrix[i,])
    
    L[i,e_predating] <- 1
    
  }
  
  return(L)
  
  
}




# ============================================================
# Matrices used for selection calculations
# ============================================================

# ------------------------------------------------------------
# Create the full square bipartite adjacency matrix ("F matrix")
# ------------------------------------------------------------
# Embeds a host × exploiter incidence matrix B into a square matrix:
#     [ 0   B ]
#     [ B'  0 ]
# used to mask trait differences only along realized interactions.
#
# Input:
#   binmat - binary incidence matrix (hosts × exploiters)
#
# Output:
#   Square adjacency matrix (nsp × nsp)



create_F <- function(binmat) {
  mat_F<-rbind (
    cbind( matrix(0, nrow=nrow(binmat), ncol=nrow(binmat)),binmat),
    cbind( t(binmat), matrix(0, nrow=ncol(binmat),ncol=ncol(binmat)))) 
  return(mat_F)
}




# ------------------------------------------------------------
# Create the "a matrix": trait differences masked by interactions
# ------------------------------------------------------------
# Computes pairwise trait differences between all species, then masks these
# differences by the adjacency matrix F (so only interacting pairs contribute):
#   a_ij = F_ij * (x_i - x_j)
#
# Inputs:
#   trait_vec - numeric vector of species traits (length nsp)
#   Fmat      - square adjacency matrix (nsp × nsp)
#
# Output:
#   Square matrix a.matrix (nsp × nsp) with trait differences on edges and 0 elsewhere

create_a_matrix <- function(trait_vec, Fmat) {
  dif.traits <- outer(trait_vec, trait_vec, FUN= "-") # dif.traits is the outer product of the arrays trait[sim,] and trait[sim,] (https://www.youtube.com/watch?v=FCmH4MqbFGs&t=492s) (https://math.stackexchange.com/questions/973559/outer-product-of-two-matrices)
  a.matrix <- Fmat * dif.traits # 
  
  return(a.matrix)
}




# ============================================================
# Selection pressures (interaction-mediated selection gradients)
# ============================================================

# ------------------------------------------------------------
# Selective pressure exerted by one exploiter on a host (victim)
# ------------------------------------------------------------
# Implements a thresholded matching rule: if the trait difference between an
# interacting exploiter and host lies within ±threshold, selection pushes the
# host away from the exploiter's trait (escape-by-mismatch). If the mismatch is
# already beyond the threshold, the exploiter imposes no selection (Sh = 0).
#
# Inputs:
#   one_trait_match - scalar trait difference (masked edge value in a.matrix)
#   threshold       - matching threshold (epsilon)
#
# Output:
#   Scalar selection contribution Sh for that exploiter???host pair


Calculate_selpres_exploiters_onhost <- function(one_trait_match, threshold) {
  
  s <- round(one_trait_match, digits=8)
  
  if (s > 0 && s < threshold) {
    #Sh[exploiter] in the origninal code. With only Sh, because I go through
    #the vector of exploiters, shoul add into the Sh vector?
    Sh <- (-1*(one_trait_match) + threshold)
    
  } else if (s < 0 && s > -threshold) {
    
    Sh <- (-1*(one_trait_match) - threshold)
    
  } else if (s == 0) {
    
    Sh <- (-1*(one_trait_match) + sample(c(-threshold, threshold),1))
    
  } else if (s >= threshold || s <= -threshold)  {
    Sh <- 0
  }
  
  return(Sh)
}





# ------------------------------------------------------------
# Total interaction selection on hosts (victims): Rh
# ------------------------------------------------------------
# For each host, extracts trait differences only to its interacting exploiters.
# Each exploiter contributes Sh based on the thresholded rule above, and host-level
# selection Rh is the (scaled) sum of those contributions (optionally normalized by
# the number of partners within the threshold).
#
# Inputs:
#   a_matrix    - square a.matrix with masked trait differences (nsp × nsp)
#   big_mat     - square adjacency (F) indicating which pairs interact (nsp × nsp)
#   KSI.b       - scalar scaling coefficient for biotic selection strength
#   n_host      - number of hosts
#   e_threshold - interaction threshold (epsilon)
#
# Output:
#   Numeric vector Rh of length n_host


calc_Rh <- function(a_matrix, big_mat, KSI.b, n_host, e_threshold){
  # # Create a matrix with each host' differences in trait values
  # (a.matrix.host) only for those interactions that happen with exploiters
  # (bigmat_F)
  
  Rh <- numeric(n_host)
  
  for (h in 1:n_host) {
    
    
    a.matrix.host.n.host <- a_matrix[h,big_mat[h,] > 0] # PROBLEM HERE
    
    
    if (length(a.matrix.host.n.host) == 0) {
      
      Rh[h] <- 0
      
    } else {
      
      n_partners <- sum(abs(round(a.matrix.host.n.host, digits=8)) <= e_threshold ) # K in f(5)
      
      Sh <- c()
      
      for (exploiter in 1:length(a.matrix.host.n.host)) {
        
        Sh[exploiter]<- Calculate_selpres_exploiters_onhost(
          one_trait_match =  a.matrix.host.n.host[exploiter], threshold = e_threshold)
        
      }
      
      
      Rh[h] <- KSI.b * sum(Sh, na.rm=TRUE) / n_partners # biot se;ection in f(5)
      # Why should I get NAs?
      Rh[h] <- replace(Rh[h],is.nan(Rh[h]),0)
    }
    
    Rh[Rh == -Inf] <- 0
    Rh[Rh == Inf] <- 0  
    
  }
  
  return(Rh)
  
}	


# ------------------------------------------------------------
# Corrected variant of Rh (no normalization by n_partners)
# ------------------------------------------------------------
# Same as calc_Rh(), but does not divide by the number of partners within
# the threshold. This treats Rh as a cumulative (not per-partner averaged)
# selection pressure.


calc_Rh_corrected <- function(a_matrix, big_mat, KSI.b, n_host, e_threshold){
  # # Create a matrix with each host' differences in trait values
  # (a.matrix.host) only for those interactions that happen with exploiters
  # (bigmat_F)
  
  Rh <- numeric(n_host)
  
  for (h in 1:n_host) {
    
    
    a.matrix.host.n.host <- a_matrix[h,big_mat[h,] > 0] # PROBLEM HERE
    
    
    if (length(a.matrix.host.n.host) == 0) {
      
      Rh[h] <- 0
      
    } else {
      
      n_partners <- sum(abs(round(a.matrix.host.n.host, digits=8)) <= e_threshold ) # K in f(5)
      
      Sh <- c()
      
      for (exploiter in 1:length(a.matrix.host.n.host)) {
        
        Sh[exploiter]<- Calculate_selpres_exploiters_onhost(
          one_trait_match =  a.matrix.host.n.host[exploiter], threshold = e_threshold)
        
      }
      
      
      Rh[h] <- KSI.b * sum(Sh, na.rm=TRUE) # biot se;ection in f(5)
      # Why should I get NAs?
      Rh[h] <- replace(Rh[h],is.nan(Rh[h]),0)
    }
    
    Rh[Rh == -Inf] <- 0
    Rh[Rh == Inf] <- 0  
    
  }
  
  return(Rh)
  
}	






# ------------------------------------------------------------
# Selective pressure exerted by each host on one exploiter: Sp (pairwise)
# ------------------------------------------------------------
# Computes the contribution of each interacting host to exploiter selection,
# weighted by a Gaussian function of mismatch. This implements the idea that
# better-matched hosts contribute more strongly to selection on exploiters.
#
# Inputs:
#   a.matrix.exploiter_nhost - vector of exploiter-host trait differences for interacting hosts
#   exp_sum                  - sum(exp(-B * diff^2)) across interacting hosts (normalization)
#   B                        - strength of Gaussian weighting by mismatch
#
# Output:
#   Vector Sp of per-host contributions to exploiter selection


Calculate_selpres_victims_onexploiter <- function(a.matrix.exploiter_nhost, 
                                                  exp_sum, B) {
  
  Sp <- numeric(length(a.matrix.exploiter_nhost)) 
  
  for (host in 1:length(a.matrix.exploiter_nhost))
  {
    Sp[host] <- (exp(-B*(a.matrix.exploiter_nhost[host]^2))/exp_sum)*
      (-1*a.matrix.exploiter_nhost[host]) #(?) why multiplying by -1?
    
  }
  return(Sp)
}



# ------------------------------------------------------------
# Total interaction selection on exploiters: Rp (original version)
# ------------------------------------------------------------
# For each exploiter, extracts trait differences only to its interacting hosts and
# computes a weighted sum of host contributions (Sp). Rp is then scaled by KSI.b
# and (in this version) divided by the exploiter's degree (number of partners).
#
# Inputs:
#   a_matrix     - a.matrix (nsp × nsp)
#   big_mat      - adjacency F (nsp × nsp)
#   KSI.b        - scalar biotic selection strength
#   B            - Gaussian weighting parameter (mismatch sensitivity)
#   n_host       - number of hosts
#   n_exploiter  - number of exploiters
#   nsp          - total number of species (n_host + n_exploiter)
#
# Output:
#   Numeric vector Rp of length n_exploiter

calc_Rp <- function(a_matrix, big_mat, KSI.b, B, n_host,n_exploiter,nsp) {
  
  Rp <- numeric(n_exploiter)
  
  for (p in as.numeric(n_host+1):nsp) {
    
    a.matrix.exploiter.n.host <- a_matrix[p,big_mat[p,] > 0]
    
    exp.sum.p <- sum(exp(-B*(a.matrix.exploiter.n.host^2))) # denominator of formula (3)
    
    Sp<-Calculate_selpres_victims_onexploiter(a.matrix.exploiter_nhost=
                                                a.matrix.exploiter.n.host, 
                                              exp_sum = exp.sum.p,
                                              B=b)
    
    
    Rp[p-n_host] <- KSI.b * (sum(Sp, na.rm=TRUE)/sum(big_mat[p,])) #biotic 
    #selection part in formula (5)
    
  }
  
  return(Rp)
  
}


# ------------------------------------------------------------
# Corrected variant of Rp (no normalization by degree)
# ------------------------------------------------------------
# Same as calc_Rp(), but does not divide by the number of interacting hosts,
# treating Rp as a cumulative selection pressure.

calc_Rp_corrected <- function(a_matrix, big_mat, KSI.b, B, n_host,n_exploiter,nsp) {
  
  Rp <- numeric(n_exploiter)
  
  for (p in as.numeric(n_host+1):nsp) {
    
    a.matrix.exploiter.n.host <- a_matrix[p,big_mat[p,] > 0]
    
    exp.sum.p <- sum(exp(-B*(a.matrix.exploiter.n.host^2))) # denominator of formula (3)
    
    Sp<-Calculate_selpres_victims_onexploiter(a.matrix.exploiter_nhost=
                                                a.matrix.exploiter.n.host, 
                                              exp_sum = exp.sum.p,
                                              B=b)
    
    
    Rp[p-n_host] <- KSI.b * sum(Sp, na.rm=TRUE) #biotic 
    #selection part in formula (5)
    
  }
  
  return(Rp)
  
}





# ============================================================
# Summary measures of mean adaptation at the end of a simulation
# ============================================================

# ------------------------------------------------------------
# Mean final "biotic adaptation" of hosts (victims)
# ------------------------------------------------------------
# Computes trait differences at the final timestep, masks by the interaction
# matrix bigmat.F, then averages (rowMeans) across victims only.
# Zeros (non-interactions) are treated as NA so that only interacting pairs
# contribute to the mean.
#
# Note: this function assumes the existence of objects nsims and bigmat.F
# in the global environment (as in the original simulation script).


compute_mean_adapt_biotsel_final_v <- function(matrix_trait, n_sims) {
  
  dif.traits.final <- outer(matrix_trait[nsims,], matrix_trait[nsims,], FUN= "-") 
  a.matrix_final <- bigmat.F*dif.traits.final
  
  a.matrix_final_nozeros<-a.matrix_final
  
  is.na(a.matrix_final_nozeros) <- a.matrix_final_nozeros==0
  
  v_Adapt.biot.sel_final<-mean(rowMeans(a.matrix_final_nozeros, na.rm=TRUE))
  
  return(v_Adapt.biot.sel_final)
}

compute_mean_adapt_biotsel_final_e <- function(matrix_trait, n_sims) {
  
  dif.traits.final <- outer(matrix_trait[nsims,], matrix_trait[nsims,], FUN= "-") 
  a.matrix_final <- bigmat.F*dif.traits.final
  
  a.matrix_final_nozeros<-a.matrix_final
  
  is.na(a.matrix_final_nozeros) <- a.matrix_final_nozeros==0
  
  v_Adapt.biot.sel_final<-mean(colMeans(a.matrix_final_nozeros, na.rm=TRUE))
  
  return(v_Adapt.biot.sel_final)
}







####################################################################

# GENERATE NETWORKS #

####################################################################


  


  
################### Simulate network built from trait matching function 
################### with certain connectance and modularity
  

######## Generate network by trait matching function  with a certain connectance





######## Simulate network with certain connectance & modularity


network_conn_mod <- function(n_host, n_exploiter, connectance, lowerMod, upperMod) {
  
  mod <- 0
  count <- 0
  
  #list_item <- list()
  
  
  # while (length(list_item[1]) < 10 ) {
  
  
  
  while (!(mod < upperMod && mod > lowerMod)) {
    
    
    ntw <- generate_net(n.host = n_host, n.exploiter = n_exploiter, C = connectance, mean.trait = 0, sd.trait = 0.1, scale_dlogis = 0.06)
    
    mod <- bipartite::networklevel(ntw, index = "modularity")
    
    count <- count + 1
    
  }
  
  list_item <- list(ntw, mod, count)
  return(list_item)
  
}

#}




######### Generate n networks, store them in list and save connectance, mod, and nest metrics


generate_n_networks <- function(n, nexploiter, nhost, connectance) {
  
  
  list_results <- list()
  
  for (i in 1:n) {
    
    print(sprintf("loop n  %s", i))
    
    
    net <- generate_network_given_conn_each_sim_allsp.int (n_exploiter = nexploiter, n_host = nhost, conn = connectance, trait_sim = rnorm(nhost + nexploiter, mean=0, sd=0.1))
    
    list_netw[[i]] <- net
    
    vec_mod_netw[i]<-networklevel(net, index = "modularity")
    
    vec_nest_netw[i]<-networklevel(net, index = "nestedness")
    
    vec_conn_netw[i]<-networklevel(net, index = "connectance")
    
  }
  
  list_results <- list(list_netw, vec_mod_netw, vec_nest_netw, vec_conn_netw)
  
  return(list_results)
}



# Function to add random changes in teta of the community
random_change_teta <- function(value_teta, SD) {
  new_teta <- sample(rnorm(1000, mean = value_teta, sd = SD), 1)
  return(new_teta)
}


# Add at least 1 interaction to disconnected species


assign_interactions <- function(binary_matrix) {
  
  
  
  # Check for disconnected species
  disconnected_rows <- which(rowSums(binary_matrix) == 0)
  disconnected_cols <- which(colSums(binary_matrix) == 0)
  
  
  if(length(disconnected_cols) > 0){
    
    for (col in 1:length(disconnected_cols)) {
      
      binary_matrix[sample(1:nrow(binary_matrix), 1), disconnected_cols[col]] <- 1
      
    }
    
    
  }
  
  if(length(disconnected_rows) > 0){
    
    
    for (row in 1:length(disconnected_rows)) {
      
      binary_matrix[disconnected_rows[row], sample(1:ncol(binary_matrix), 1)] <- 1
      
    }
    
    
  }
  

  return(binary_matrix)
}






back_transform_sqrt <- function(x) {
  return(x^2)  # Square back the transformed value
}

  







