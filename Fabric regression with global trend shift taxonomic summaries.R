# Taxonomic summaries of significant Fabric node effects

library(geiger)
library(dplyr)

setwd("C:/Users/hy929145/Documents/Cactus height")

beta_df <- read.csv("Cactus_Fabric_reg_GT_beta_effects.csv")
beta_df$X <- NULL
scalar_df <- read.csv("Cactus_Fabric_reg_GT_scalar_effects.csv")
scalar_df$X <- NULL
head(beta_df)
head(scalar_df)
tree <- readRDS("Cactus_Fabric_reg_GT_time_tree.rds")

# Descendant species and genera for each shift
beta_df$descendant_species <- sapply(beta_df$ape_node, function(x) paste(tips(tree, x), collapse = "; "))
beta_df$descendant_genera <- sapply(beta_df$ape_node, function(x) paste(unique(sub("_.*", "", tips(tree, x))), collapse = "; "))
beta_df$n_descendant_species <- sapply(beta_df$ape_node, function(x) length(tips(tree, x)))
beta_df$n_descendant_genera <- sapply(beta_df$ape_node, function(x) length(unique(sub("_.*", "", tips(tree, x)))))

scalar_df$descendant_species <- sapply(scalar_df$ape_node, function(x) paste(tips(tree, x), collapse = "; "))
scalar_df$descendant_genera <- sapply(scalar_df$ape_node, function(x) paste(unique(sub("_.*", "", tips(tree, x))), collapse = "; "))
scalar_df$n_descendant_species <- sapply(scalar_df$ape_node, function(x) length(tips(tree, x)))
scalar_df$n_descendant_genera <- sapply(scalar_df$ape_node, function(x) length(unique(sub("_.*", "", tips(tree, x)))))

# Assign tribes and subfamilies
taxonomy <- read.csv("Cactus_fabric_taxonomy.csv")
taxonomy$X <- NULL

genus_tax <- taxonomy %>%
  dplyr::select(genus, tribe, subfamily) %>%
  distinct()

scalar_df <- scalar_df %>%
  rowwise() %>%
  mutate(tribe = paste(unique(genus_tax$tribe[match(strsplit(descendant_genera, "; ")[[1]], genus_tax$genus)]), collapse = "; "),
         subfamily = paste(unique(genus_tax$subfamily[match(strsplit(descendant_genera, "; ")[[1]], genus_tax$genus)]), collapse = "; ")) %>%
  ungroup()

beta_df <- beta_df %>%
  rowwise() %>%
  mutate(tribe = paste(unique(na.omit(genus_tax$tribe[match(strsplit(descendant_genera, "; ")[[1]], genus_tax$genus)])), collapse = "; "),
         subfamily = paste(unique(na.omit(genus_tax$subfamily[match(strsplit(descendant_genera, "; ")[[1]], genus_tax$genus)])), collapse = "; ")) %>%
  ungroup()

write.csv(scalar_df, "Cactus_Fabric_reg_GT_scalar_effects_descendants.csv")
write.csv(beta_df, "Cactus_Fabric_reg_GT_beta_effects_descendants.csv")

table(scalar_df$tribe)
table(scalar_df$subfamily)
table(beta_df$tribe)
table(beta_df$subfamily)
