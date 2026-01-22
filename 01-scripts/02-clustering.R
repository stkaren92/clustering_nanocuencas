library(tidyverse)
library(factoextra)
library(dad)

# Utils functions ----

#' Calculate Jensen-Shannon Divergence between two probability distributions
#'
#' @param p First probability distribution
#' @param q Second probability distribution  
#' @param base Logarithm base (default: 2, giving result in bits)
#' @return Jensen-Shannon divergence between p and q
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' q <- c(0.2, 0.3, 0.5)
#' js_divergence(p, q)
js_divergence <- function(p, q, base = 2) {
  # Input validation and normalization
  if (length(p) != length(q)) {
    stop("Probability distributions must have the same length")
  }
  
  if (abs(sum(p) - 1) > 1e-10 || abs(sum(q) - 1) > 1e-10) {
    warning("Input distributions don't sum to 1. Normalizing...")
    p <- p / sum(p)
    q <- q / sum(q)
  }
  
  # Remove zero probabilities to avoid log(0)
  epsilon <- 1e-10
  p_safe <- p + epsilon
  q_safe <- q + epsilon
  p_safe <- p_safe / sum(p_safe)
  q_safe <- q_safe / sum(q_safe)
  
  # Calculate the average distribution
  m <- 0.5 * (p_safe + q_safe)
  
  # Calculate KL divergences
  kl_pm <- sum(p_safe * log(p_safe / m, base = base))
  kl_qm <- sum(q_safe * log(q_safe / m, base = base))
  
  # Jensen-Shannon divergence
  jsd <- 0.5 * kl_pm + 0.5 * kl_qm
  
  return(jsd)
}

#' Calculate the square root of Jensen-Shannon Divergence
#'
#' @param p First probability distribution
#' @param q Second probability distribution
#' @param base Logarithm base (default: 2)
#' @return Square root of Jensen-Shannon divergence
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' q <- c(0.2, 0.3, 0.5)
#' js_distance(p, q)
js_distance <- function(p, q, base = 2) {
  return(sqrt(js_divergence(p, q, base)))
}

#' Compute a distance matrix
#'
#' @description
#' Calculates a pairwise distance matrix between observations using variables
#' represented as a probability distribution or as numeric parameters (mean, sd).
#' For each pair of rows, the function:
#' * Computes a Jensen–Shannon Divergence (`js_distance`) for variables classified
#'   as `"distribution"`
#' * Computes a Hellinger distance (`hellingerpar`) for variables classified
#'   as `"parameters"`
#'
#' The total distance between two observations is the sum of all
#' variable-specific distances.
#'
#' @param X A data frame containing the variables used to compute
#'   distances. Each row corresponds to one observation and each column is a
#'   variable represented either as a probability distribution (vector that sums
#'   to 1) or as numeric parameters (vector with (mean, sd)).
#' @param variable_type A character vector describing the type of each column
#'   in `X`. Must have the same length as `ncol(X)`. Valid values:
#'   * `"distribution"` – column contains a probability distribution;
#'     Jensen–Shannon Divergence is applied.
#'   * `"parameters"` – column contains numeric parameters (mean, sd);
#'     Hellinger distance on parameterized distributions is applied.
#'
#' @return
#' A numeric matrix with dimensions `nrow(X) × nrow(X)` where cell `(i, j)`
#' represents the total distance between observations `i` and `j`.
#'
#' @examples
#' X <- data.frame(
#'   a = list(c(0.1, 0.9), c(0.15, 0.85), c(5.1, 2.2)),
#'   b = list(c(0.3, 0.7), c(0.35, 0.65), c(4.2, 1.8))
#' )
#' variable_type <- c("distribution", "distribution", "parameters")
#' S <- get_distance_matrix(X, variable_type)
#'
get_distance_matrix <- function(X, variable_type){
  
  # Number of observations (rows) and variables (columns)
  n <- nrow(X)
  n_var <- ncol(X)
  # Initialize an empty n × n matrix to store distances
  S <- matrix(0, nrow = n, ncol = n)
  
  # Loop over all row pairs (i, j)
  for (i in 1:n){
    for(j in 1:n){
      # Extract row vectors for observation i and j
      a <- X[i,]
      b <- X[j,]
      dist <- 0 # initialize distance accumulator
      
      # Loop over variables
      for (z in 1:n_var) {
        # Extract the z-th variable
        p <- unlist(a[z], use.names=FALSE)
        q <- unlist(b[z], use.names=FALSE)
        
        # -------- Variable type: categorical distribution --------
        if(variable_type[z] == "distribution") {
          # Jensen–Shannon distance between probability vectors
          dist_z <- js_distance(p, q)
          dist <- dist + dist_z # Add to total distance
        
          # -------- Variable type: parameters (mean & sd) --------  
        } else if (variable_type[z] == "parameters"){
          # Hellinger distance between parameterized distributions
          # Note: second parameter is variance, so sd^2
          dist_z <- hellingerpar(p[1], p[2]^2,
                                 q[1], q[2]^2)
          # Add only if the value is not missing
          if(!is.na(dist_z)){
            dist <- dist + dist_z
          }
        }
      }
      # Store total distance in the matrix
      S[i,j] <- dist
    }
    cat(sprintf("\rProgreso: %.1f%%", round((i / n) * 100, 1))) # Print progress
    flush.console() 
  }
  # Return the completed pairwise distance matrix
  return(S)
}


# Clustering ----
INPUT_DIR <- '02-processed_data'
OUTPUT_DIR <- '03-clustering_output'
current_date <- now() %>% date()

data <- read_csv(fs::path_join(c(INPUT_DIR, "2026-01-12_dataset.csv")))

# Normalize count variables to get a probability distribution
data <- data %>%
  rowwise() %>%
  mutate(across(starts_with('glg'),
                ~ . / sum(c_across(starts_with('glg')))),
         across(starts_with('suelo1'),
                ~ . / sum(c_across(starts_with('suelo1')))),
         across(starts_with('suelo2'),
                ~ . / sum(c_across(starts_with('suelo2')))),
         across(starts_with('usv'),
                ~ . / sum(c_across(starts_with('usv'))))
  ) %>%
  ungroup()

# Group columns in vector variables
data <- data %>%
  rowwise() %>%
  mutate(tipo_rocas = list(c_across(starts_with('glg'))),
         tipo_suelo_gen = list(c_across(starts_with('suelo1'))),
         tipo_suelo_gen2 = list(c_across(starts_with('suelo2'))),
         uso_suelo = list(c_across(starts_with('usv'))),
         pendiente = list(c_across(all_of(c("slope_mean","slope_sd"))))
         ) %>%
  ungroup()

# Get distance matrix 
X <- data %>% 
  dplyr::select(c("tipo_rocas",
                  "tipo_suelo_gen",
                  "tipo_suelo_gen2",
                  "uso_suelo",
                  "pendiente"))

variable_type <- c("distribution",
                   "distribution",
                   "distribution",
                   "distribution",
                   "parameters")

S <- get_distance_matrix(X, variable_type)
rownames(S) <- data %>% pull(clave_sht)
colnames(S) <- data %>% pull(clave_sht)

# Select optimal k
fviz_nbclust(S, hcut, method = "wss")
fviz_nbclust(S, hcut, method = "silhouette",
             hc_method = "complete")

# Clustering
k <- 4
hclust <- hclust(as.dist(S), method = 'complete')
cut <- cutree(hclust, k = k)
data$cluster <- cut
data <- data %>% 
  relocate(cluster)

# Save dandogram
jpeg(fs::path_join(c(OUTPUT_DIR, paste(current_date, "dendrogram.jpg", sep = "_"))),
     width = 1800, height = 600, quality = 90)
plot(hclust)
rect.hclust(hclust , k = k, border = 2:6)
dev.off() # Close the device and save the file

# Description by cluster 
data_sum <- data %>% 
  group_by(cluster) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  pivot_longer(
    cols = -c("cluster"), 
    names_to = "category",           
    values_to = "value"           
  ) %>% 
  separate_wider_delim(category, 
                       delim = "_", 
                       names = c('var', 'category'), too_many = 'merge')

ggplot(data_sum %>% 
         filter(var != "slope"), 
       aes(x=category, y=value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(var ~ cluster, scales="free_y",
             space = "free")
ggsave(fs::path_join(c(OUTPUT_DIR, paste(current_date, "cluster_description.jpg", sep = "_"))),
       width = 30,
       height = 50,
       units = "cm")

# Save cluster shapefile
nanocuecas_data <- sf::read_sf('00-raw_data/nanocuencas/SHT_BAconsensuadov2_uw.shp')

nanocuecas_data <- nanocuecas_data %>% 
  left_join(data %>% 
              select(cluster,clave_sht), 
            by = c("CLAVE_SHT"="clave_sht"))

ggplot(nanocuecas_data, 
       aes(fill = as.factor(cluster))) +
  geom_sf()

sf::write_sf(nanocuecas_data,
             fs::path_join(c(OUTPUT_DIR,
                             paste(current_date, 
                                   "nanocuecas_cluster.shp", 
                                   sep = "_"))))