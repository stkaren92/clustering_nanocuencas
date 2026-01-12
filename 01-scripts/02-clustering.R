library(tidyverse)

# Utils functions ----

# Hellinger Distance
hellinger_discrete <- function(p, q) {
  # p <- p / sum(p)
  # q <- q / sum(q)
  sqrt( sum( (sqrt(p) - sqrt(q))^2 ) / 2 )
}

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

#' Calculate the square root of Jensen-Shannon Divergence (a proper metric)
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

# Clustering ----
OUTPUT_DIR <- '02-processed_data'
current_date <- now() %>% date()

data <- read_csv(fs::path_join(c(OUTPUT_DIR, "2026-01-12_dataset.csv")))

# Categories in each variable
tipo_rocas_var <- starts_with('glg')
# permeabilidad_var <- names(df)[22:31]
tipo_suelo_gen_var <- starts_with('suelo1')
tipo_suelo_gen2_var <- starts_with('suelo2')
# tipo_suelo_esp_var <- names(df)[63:88]
uso_suelo_var <- starts_with('usv')

# Normalize by variable
data <- data %>%
  rowwise() %>%
  mutate(across(starts_with('glg'),
                ~ . / sum(c_across(starts_with('glg')))),
         # across(all_of(permeabilidad_var),
         #         ~ . / sum(c_across(all_of(permeabilidad_var)))),
         across(starts_with('suelo1'),
                ~ . / sum(c_across(starts_with('suelo1')))),
         across(starts_with('suelo2'),
                ~ . / sum(c_across(starts_with('suelo2')))),
         # across(all_of(tipo_suelo_esp_var),
         #        ~ . / sum(c_across(all_of(tipo_suelo_esp_var)))),
         across(starts_with('usv'),
                ~ . / sum(c_across(starts_with('usv'))))
  ) %>%
  ungroup()

# Group categories in vector variables
data <- data %>%
  rowwise() %>%
  mutate(tipo_rocas = list(c_across(starts_with('glg'))),
         # permeabilidad = list(c_across(all_of(permeabilidad_var))),
         tipo_suelo_gen = list(c_across(starts_with('suelo1'))),
         tipo_suelo_gen2 = list(c_across(starts_with('suelo2'))),
         # tipo_suelo_esp = list(c_across(all_of(tipo_suelo_esp_var))),
         uso_suelo = list(c_across(starts_with('usv')))) %>%
  ungroup()

# Get distance matrix 
X <- data %>% 
  dplyr::select(c("tipo_rocas",
                  # "permeabilidad",
                  "tipo_suelo_gen",
                  "tipo_suelo_gen2",
                  # "tipo_suelo_esp",
                  "uso_suelo"))

n <- nrow(X)
n_c <- ncol(X)
S <- matrix(0, nrow = n, ncol = n)

for (i in 1:n){
  for(j in 1:n){
    a <- X[i,]
    b <- X[j,]
    dist <- 0
    for (z in 1:n_c) {
      dist_z <- js_distance(unlist(a[z], use.names=FALSE), 
                            unlist(b[z], use.names=FALSE))
      dist <- dist + dist_z
    }
    S[i,j] <- dist
  }
  cat(sprintf("\rProgreso: %.1f%%", round((i / n) * 100, 1))) 
  flush.console() 
}

rownames(S) <- data %>% pull(clave_sht)
colnames(S) <- data %>% pull(clave_sht)

# Clustering 
k <- 4
hclust <- hclust(as.dist(S), method = 'complete')
cut <- cutree(hclust, k = k)
data$cluster <- cut
data <- data %>% 
  relocate(cluster)

jpeg(fs::path_join(c(OUTPUT_DIR, paste(current_date, "dendrogram.jpg", sep = "_"))), 
     width = 1800, height = 600, quality = 90)
plot(hclust)
rect.hclust(hclust , k = k, border = 2:6)
dev.off() # Close the device and save the file

# Description by cluster 
# Mean by cluster
data_sum <- data %>% 
  group_by(cluster) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Reshape
data_sum <- data_sum %>%
  pivot_longer(
    cols = -c("cluster"), 
    names_to = "category",           
    values_to = "value"           
  )

data_sum <- data_sum %>% 
  separate_wider_delim(category, 
                       delim = "_", 
                       names = c('var', 'category'), too_many = 'merge')


ggplot(data_sum, 
       aes(x=category, y=value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(var ~ cluster, scales="free",
             space = "free")
ggsave(fs::path_join(c(OUTPUT_DIR, paste(current_date, "cluster_description.jpg", sep = "_"))),
       width = 30,
       height = 50,
       units = "cm")
