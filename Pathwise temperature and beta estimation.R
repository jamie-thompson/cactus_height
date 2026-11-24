# Cactus pathwise temperature and beta estimation

library(ape)
library(castor)
library(dplyr)
library(phytools)
library(tibble)

setwd("C:/Users/hy929145/Documents/Cactus height")

# Pathwise temperature change
geo_tree <- readNexus("Cacti_BT_geo_THH_tree.nex.trees")
anc_temp <- read.csv("Full_paleotemp_cacti.csv")

root_node <- Ntip(geo_tree) + 1
paths <- lapply(seq_len(Ntip(geo_tree)), function(i) {
  nodepath(geo_tree, from = root_node, to = i)
})
names(paths) <- geo_tree$tip.label

temp_medians <- anc_temp %>%
  group_by(ape_node) %>%
  summarise(median_temperature = median(temperature, na.rm = TRUE), .groups = "drop")

node_climate <- temp_medians

pathwise_climate_change <- lapply(names(paths), function(sp) {
  node_path <- paths[[sp]]
  climate <- node_climate %>%
    filter(ape_node %in% node_path) %>%
    mutate(path_order = match(ape_node, node_path)) %>%
    arrange(path_order)
  temp_change <- diff(climate$median_temperature)
  temp_abs <- sum(abs(temp_change), na.rm = TRUE)
  temp_pos <- sum(temp_change[temp_change > 0], na.rm = TRUE)
  temp_neg <- sum(temp_change[temp_change < 0], na.rm = TRUE)
  temp_net <- sum(temp_change, na.rm = TRUE)
  tibble(species = sp,
    temp_abs = temp_abs,
    log_temp_abs = log1p(temp_abs),
    temp_pos = temp_pos,
    log_temp_pos = log1p(temp_pos),
    temp_neg = temp_neg,
    log_temp_neg = log1p(abs(temp_neg)),
    temp_net = temp_net,
    log_temp_net = ifelse(temp_net != 0, log(abs(temp_net)), NA_real_))
}) %>% bind_rows()

head(pathwise_climate_change)
write.csv(pathwise_climate_change, "Cacti_pathwise_climate_changes_July_log1p.csv", row.names = FALSE)

# Pathwise beta effects
fabric_dir <- file.path(getwd(), "cactus_height_Fabric_reg_July_gt")
tree <- read.nexus("Cactus_height_gf_BT_tree_July_correction_final.nex.prue.trees")
pp_df <- read.csv(file.path(fabric_dir, "cactus_height_Fabric_reg_July_gt_master.csv"), header = TRUE)
cactus_height <- read.table("cactus_height_Fabric_reg_July/Data-0000000.txt", sep = "\t", header = FALSE)
cactus_height = cactus_height[c("V1", "V2")]
colnames(cactus_height) <- c("species", "max_size")

cactus_height <- cactus_height %>%
  rename(plant_height = max_size) %>%
  filter(species %in% tree$tip.label)

pp_df$ape_node <- NA
for (i in seq_len(nrow(pp_df))) {
  taxa_list <- as.character(unlist(pp_df[i, grep("^Taxa\\.", colnames(pp_df))]))
  taxa_list <- taxa_list[!is.na(taxa_list)]
  taxa_list <- trimws(taxa_list)
  taxa_list <- intersect(taxa_list, tree$tip.label)
  mrca_node <- tryCatch({
    get_mrca_of_set(tree, taxa_list)
  }, error = function(e) {
    warning(paste("Error in row", i, "with taxa:", paste(taxa_list, collapse = ", ")))
    NA
  })
  pp_df$ape_node[i] <- mrca_node
}

beta_signed_lookup <- setNames(pp_df$Mean..Beta...BL..NZ, pp_df$ape_node)

get_lineage_nodes <- function(tree, tip_node) {
  current_node <- tip_node
  lineage_nodes <- integer(0)
  while (current_node != Ntip(tree) + 1) {
    edge_row <- which(tree$edge[, 2] == current_node)
    lineage_nodes <- c(current_node, lineage_nodes)
    current_node <- tree$edge[edge_row, 1]
  }
  lineage_nodes
}

beta_history_df <- lapply(seq_along(tree$tip.label), function(i) {
  lineage_nodes <- get_lineage_nodes(tree, i)
  effects <- beta_signed_lookup[as.character(lineage_nodes)]
  effects <- as.numeric(effects)
  effects <- effects[!is.na(effects)]
  dir_abs <- sum(abs(effects))
  dir_net <- sum(effects)
  dir_pos <- sum(effects[effects > 0])
  dir_neg <- sum(effects[effects < 0])
  tibble(species = tree$tip.label[i],
    dir_abs = dir_abs,
    dir_net = dir_net,
    dir_pos = dir_pos,
    dir_neg = dir_neg,
    log_dir_abs = log1p(dir_abs),
    log_dir_pos = log1p(dir_pos),
    log_dir_neg = log1p(abs(dir_neg)),
    n_dir_pos = sum(effects > 0),
    n_dir_neg = sum(effects < 0),
    n_dir_total = length(effects))
}) %>% bind_rows()

head(beta_history_df)
write.csv(beta_history_df, "Cacti_pathwise_beta_July_log1p.csv", row.names = FALSE)