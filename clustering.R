library(tidyverse)
library(readxl)
library(statip)

# Hellinger Distance
hellinger_discrete <- function(p, q) {
  # p <- p / sum(p)
  # q <- q / sum(q)
  sqrt( sum( (sqrt(p) - sqrt(q))^2 ) / 2 )
}

# Read file
df <- read_excel("data/matriz_v0.xlsx", "matriz v0")

# Categories by variable
tipo_rocas_var <- names(df)[4:21]
permeabilidad_var <- names(df)[22:31]
tipo_suelo_gen_var <- names(df)[32:46]
tipo_suelo_gen2_var <- names(df)[47:62]
tipo_suelo_esp_var <- names(df)[63:88]
uso_suelo_var <- names(df)[89:131]


#### Normalize ####
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

#### Create vectors ####
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
      dist_z <- hellinger_discrete(unlist(a[z], use.names=FALSE), 
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
k <- 6
hclust <- hclust(as.dist(S), method = 'complete')
cut <- cutree(hclust, k = k)
df$cluster <- cut
df <- df %>% 
  relocate(cluster)

plot(hclust)
rect.hclust(hclust_avg , k = k, border = 2:6)
abline(h = 3, col = 'red')


df_sum <- df %>% 
  group_by(cluster) %>% 
  

# which(S == min(S), arr.ind = TRUE)
# x <- X[1,]
# y <- X[6,]
# 
# p <- unlist(x[3], use.names=FALSE)
# q <- unlist(y[3], use.names=FALSE)
# 
# hellinger(p,q)
# hellinger_discrete(p,q)
