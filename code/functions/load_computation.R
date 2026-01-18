# ------------------------------------------------------------
# Fitness (load) functions for exploiters and victims
# ------------------------------------------------------------
# These functions compute "fitness load" terms (lower is better) by decomposing
# selection into abiotic mismatch, biotic mismatch with interaction partners,
# and (for exploiters) an additional network term capturing heterogeneity among
# partner traits. Functions return both the total load and its components to
# allow diagnostics and plotting.



# ============================================================
# EXPLOITERS
# ============================================================
# Total exploiter load = abiotic + biotic + network
# - Abiotic load: squared distance between exploiter trait (x_s) and abiotic optimum (x_0)
# - Biotic load: squared distance between x_s and the mean trait of its partners (x_optbiotic)
# - Network load: increases when partners are trait-diverse (pairwise trait dispersion)
#
# Inputs:
#   x_s        - focal exploiter trait value
#   x_0        - abiotic optimum trait value
#   x_partners - numeric vector of partner (victim) trait values for this exploiter
#

fitness_load_exploiters <- function(x_s, x_0, x_partners) {
  # Abiotic load
  abiotic_term <- (x_s - x_0)^2
  
  # Biotic load
  x_optbiotic <- mean(x_partners)  # Mean trait value of partners
  biotic_term <- (x_s - x_optbiotic)^2
  
  # Network load
  network_term <- 0.5 * sum(outer(x_partners, x_partners, FUN = "-")^2) / (length(x_partners)^2)
  
  # Total load
  total_load <- abiotic_term + biotic_term + network_term
  list("load" = total_load, 
       "abiotic_L" = abiotic_term, 
       "biotic_L" = biotic_term, 
       "network_L" = network_term)
}


# ============================================================
# VICTIMS (threshold / escape-by-mismatch)
# ============================================================
# Victim load combines:
# - Abiotic load: squared distance to abiotic optimum (x_0)
# - Biotic load: imposed only by exploiters whose traits are within an interaction
#   "matching window" of size epsilon around the victim trait (x_s). If no exploiters
#   fall within this window, the victim experiences no biotic load (escape).
#
# Interpretation:
#   Victims are penalized when exploiters are close enough in trait space to interact;
#   they can reduce biotic load by moving outside exploiter matching ranges.


# Abiotic load for victims: squared distance from abiotic optimum
# Inputs:
#   x_s - focal victim trait value
#   x_0 - abiotic optimum
# Output:

#   Scalar abiotic load

compute_abiotic_load_victims <- function(x_s, x_0) {
  (x_s - x_0)^2  # Penalizes distance from abiotic optimum
}


# Biotic load for victims with a matching threshold epsilon
# Steps:
#   1) Keep only exploiters within epsilon of the victim trait (potential interactors)
#   2) If none exist, return 0 (no interaction pressure)
#   3) Otherwise, compute a penalty based on distance to the nearest edge of each exploiter's
#      interaction window (x_i ± epsilon), and average across in-range exploiters
#
# Inputs:
#   x_s        - focal victim trait value
#   x_partners - numeric vector of exploiter trait values
#   epsilon    - matching threshold (interaction window half-width)
#
# Output:

#   Scalar biotic load

compute_biotic_load_victims <- function(x_s, x_partners, epsilon) {
  # Filter exploiters by trait proximity
  x_in_range <- x_partners[abs(x_partners - x_s) <= epsilon]
  
  if (length(x_in_range) == 0) {
    return(0)  # No exploiters within range, no biotic load
  }
  
  # Compute minimum distance to either (x_i - epsilon) or (x_i + epsilon)
  mean(sapply(x_in_range, function(x_i) min((x_s - (x_i - epsilon))^2, (x_s - (x_i + epsilon))^2)))
}

# Total victim load = abiotic + biotic (no explicit network term here)
# Inputs:
#   x_s        - focal victim trait value
#   x_0        - abiotic optimum
#   x_partners - numeric vector of exploiter trait values
#   epsilon    - matching threshold defining exploiter interaction range
#
# Output:
#   List with total load plus abiotic and biotic components


fitness_load_victims <- function(x_s, x_0, x_partners, epsilon) {
  # Compute abiotic load
  abiotic_term <- compute_abiotic_load_victims(x_s, x_0)
  
  # Compute biotic load
  biotic_term <- compute_biotic_load_victims(x_s, x_partners, epsilon)
  
  # Total load
  total_load <- abiotic_term + biotic_term
  
  list(
    "load" = total_load,
    "abiotic_L" = abiotic_term,
    "biotic_L" = biotic_term
  )
}


