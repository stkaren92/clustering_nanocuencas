library(tidyverse)
library(readxl)

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

# Read file
df <- read_excel("data/matriz_v0.xlsx", "matriz v0")

#### Preprocess ####
# Categories in each variable
tipo_rocas_var <- names(df)[4:21]
permeabilidad_var <- names(df)[22:31]
tipo_suelo_gen_var <- names(df)[32:46]
tipo_suelo_gen2_var <- names(df)[47:62]
tipo_suelo_esp_var <- names(df)[63:88]
uso_suelo_var <- names(df)[89:131]

# Normalize by variable
df <- df %>%
  rowwise() %>%
  mutate(across(all_of(tipo_rocas_var),
                ~ . / sum(c_across(all_of(tipo_rocas_var)))),
         across(all_of(permeabilidad_var),
                ~ . / sum(c_across(all_of(permeabilidad_var)))),
         across(all_of(tipo_suelo_gen_var),
                ~ . / sum(c_across(all_of(tipo_suelo_gen_var)))),
         across(all_of(tipo_suelo_gen2_var),
                ~ . / sum(c_across(all_of(tipo_suelo_gen2_var)))),
         across(all_of(tipo_suelo_esp_var),
                ~ . / sum(c_across(all_of(tipo_suelo_esp_var)))),
         across(all_of(uso_suelo_var),
                ~ . / sum(c_across(all_of(uso_suelo_var))))
  ) %>%
  ungroup()

# Group categories in vector variables
df <- df %>%
  rowwise() %>%
  mutate(tipo_rocas = list(c_across(all_of(tipo_rocas_var))),
         permeabilidad = list(c_across(all_of(permeabilidad_var))),
         tipo_suelo_gen = list(c_across(all_of(tipo_suelo_gen_var))),
         tipo_suelo_gen2 = list(c_across(all_of(tipo_suelo_gen2_var))),
         tipo_suelo_esp = list(c_across(all_of(tipo_suelo_esp_var))),
         uso_suelo = list(c_across(all_of(uso_suelo_var)))) %>%
  ungroup()

#### Get distance matrix ####
X <- df %>% 
  dplyr::select(c("tipo_rocas",
                  "permeabilidad",
                  "tipo_suelo_gen",
                  "tipo_suelo_gen2",
                  "tipo_suelo_esp",
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
rownames(S) <- df$CLAVE_SHT
colnames(S) <- df$CLAVE_SHT

#### Clustering ####
k <- 4
hclust <- hclust(as.dist(S), method = 'complete')
cut <- cutree(hclust, k = k)
df$cluster <- cut
df <- df %>% 
  relocate(cluster)

jpeg("dendrogram.jpg", width = 1800, height = 600, quality = 90)
plot(hclust)
rect.hclust(hclust , k = k, border = 2:6)
dev.off() # Close the device and save the file

#### Description by cluster ####
# Mean by cluster
df_sum <- df %>% 
  group_by(cluster) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  select(-c("sup_ha_SHT",
            "id_SHT"))

# Reshape
df_sum <- df_sum %>%
  pivot_longer(
    cols = -c("cluster"), 
    names_to = "category",           
    values_to = "value"           
  )
df_sum <- df_sum %>% 
  mutate(var = case_when(category %in% tipo_rocas_var ~ "tipo_rocas", 
                         category %in% permeabilidad_var ~ "permeabilidad",
                         category %in% tipo_suelo_gen_var ~ "tipo_suelo_gen",
                         category %in% tipo_suelo_gen2_var ~ "tipo_suelo_gen2",
                         category %in% tipo_suelo_esp_var ~ "tipo_suelo_esp",
                         category %in% uso_suelo_var ~ "uso_suelo"))


ggplot(df_sum, 
       aes(x=category, y=value)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(var ~ cluster, scales="free",
             space = "free")
ggsave("cluster_description.jpg",
       width = 30,
       height = 50,
       units = "cm")

# ggplot(df_sum %>% 
#          filter(category %in% tipo_suelo_esp_var), 
#        aes(x=category, y=value, 
#            fill = as.factor(cluster))) + 
#   geom_bar(position="dodge", stat = "identity") +
#   coord_flip()
# 
# ggplot(df_sum %>% 
#          filter(category %in% tipo_suelo_esp_var), 
#        aes(x=category, y=value)) + 
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   facet_wrap(~cluster)

# which(S == max(S), arr.ind = TRUE)
# x <- X[159,]
# y <- X[27,]
# 
# p <- unlist(x[3], use.names=FALSE)
# q <- unlist(y[3], use.names=FALSE)
# 
# hellinger_discrete(p,q)
# js_distance(p,q)
