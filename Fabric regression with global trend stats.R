library(ggplot2)
library(tidyverse)
library(ape)
library(castor)
library(dplyr)
library(ggnewscale)
library(ggtree)
library(ggtreeExtra)
library(reshape2)
library(viridis)

setwd("C:/Users/hy929145/Documents/Cactus height")

# Non-Fabric regressions
read_bt_log <- function(f) {
  lines <- readLines(f)
  start <- grep("^Iteration\\s+Lh", lines)
  if (length(start) == 0) stop(paste("Header not found:", f))
  df <- read.delim(text = lines[start:length(lines)], header = TRUE, sep = "\t", na.strings = "N\\A", check.names = FALSE)
  colnames(df) <- make.names(colnames(df))
  df$run <- basename(f)
  return(df)
}

log_non_fabric <- map_dfr(1:6, ~ read_bt_log(paste0("Cactus_height_ind_cont_regression_no_barrel_", .x, ".Log.txt")))

ggplot(log_non_fabric, aes(x = Iteration, y = Lh, colour = run)) +
  geom_line(alpha = 0.7) +
  theme_bw()

ggplot(log_non_fabric, aes(x = Iteration, y = Var)) +
  geom_line(alpha = 0.5) +
  theme_bw()

ggplot(log_non_fabric, aes(x = Iteration, y = Alpha)) +
  geom_line(alpha = 0.7) +
  theme_bw()

median(log_non_fabric$R.2)

beta_cols <- paste0("Beta.", 1:8)
beta_data_non_fabric <- log_non_fabric[, beta_cols]
colnames(beta_data_non_fabric) <- c("arborescent", "columnar", "cushion_forming", "Epiphytic", "geophytic", "globose_caespitose", "globose_solitary", "shrubby")
beta_long_non_fabric <- melt(beta_data_non_fabric)
colnames(beta_long_non_fabric) <- c("Parameter", "Value")

crossing_summary_non_fabric <- aggregate(Value ~ Parameter, data = beta_long_non_fabric, function(x) {
  p_pos <- mean(x > 0, na.rm = TRUE)
  p_neg <- mean(x < 0, na.rm = TRUE)
  c(p_pos = p_pos, p_neg = p_neg, cross = min(p_pos, p_neg))
})

crossing_summary_non_fabric <- data.frame(
  Parameter = crossing_summary_non_fabric$Parameter,
  p_pos = crossing_summary_non_fabric$Value[, "p_pos"],
  p_neg = crossing_summary_non_fabric$Value[, "p_neg"],
  cross = crossing_summary_non_fabric$Value[, "cross"]
)

low_crossing_non_fabric <- crossing_summary_non_fabric[crossing_summary_non_fabric$cross < 0.05, ]
low_crossing_non_fabric

ggplot(beta_long_non_fabric, aes(x = Parameter, y = Value)) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal() +
  labs(x = NULL, y = "Value")

beta_mat_non_fabric <- as.data.frame(beta_data_non_fabric)

pairwise_beta_tests_non_fabric <- combn(colnames(beta_mat_non_fabric), 2, function(pair) {
  a <- pair[1]
  b <- pair[2]
  diff_ab <- beta_mat_non_fabric[[a]] - beta_mat_non_fabric[[b]]
  p_gt_0 <- mean(diff_ab > 0, na.rm = TRUE)
  p_lt_0 <- mean(diff_ab < 0, na.rm = TRUE)
  p_two_sided <- 2 * min(p_gt_0, p_lt_0)
  data.frame(beta_1 = a,
    beta_2 = b,
    mean_diff = mean(diff_ab, na.rm = TRUE),
    median_diff = median(diff_ab, na.rm = TRUE),
    lower_95 = quantile(diff_ab, 0.025, na.rm = TRUE),
    upper_95 = quantile(diff_ab, 0.975, na.rm = TRUE),
    p_gt_0 = p_gt_0,
    p_lt_0 = p_lt_0,
    p_two_sided = p_two_sided,
    significant = p_two_sided < 0.05)
  }, simplify = FALSE) %>%
  bind_rows()

pairwise_beta_tests_non_fabric
write.csv(pairwise_beta_tests_non_fabric, "Cactus_height_non_Fabric_betas_differences.csv")

# Fabric regressions
log_dir <- "C:/Users/hy929145/Documents/Cactus height/cactus_height_Fabric_reg_July_gt"

files <- list.files(log_dir, pattern = "\\.Log\\.txt$", full.names = TRUE)

read_bt_log <- function(f) {
  lines <- readLines(f)
  start <- grep("^Iteration\\s+Lh", lines)
  if (length(start) == 0) stop(paste("Header not found:", f))
  header <- strsplit(lines[start], "\t")[[1]]
  df <- read.delim(text = paste(lines[(start + 1):length(lines)], collapse = "\n"), header = FALSE, sep = "\t", na.strings = "N\\A", check.names = FALSE, fill = TRUE)
  if (ncol(df) > length(header)) df <- df[, seq_along(header)]
  names(df) <- make.names(header, unique = TRUE)
  df$run <- basename(f)
  df
}

logs <- purrr::map_dfr(files, read_bt_log)

ggplot(logs, aes(x = Iteration, y = Lh, colour = run)) +
  geom_line(alpha = 0.7) +
  theme_bw()

logs2 <- logs %>%
  arrange(run, Iteration) %>%
  group_by(run) %>%
  mutate(iter_rel = Iteration - min(Iteration)) %>%
  ungroup() %>%
  mutate(run_id = as.numeric(factor(run))) %>%
  arrange(run_id, iter_rel) %>%
  group_by(run_id) %>%
  mutate(offset = max(iter_rel)) %>%
  ungroup()

median(logs$R.2)

offsets <- logs2 %>%
  group_by(run_id) %>%
  summarise(max_iter = max(iter_rel)) %>%
  mutate(cum_offset = cumsum(lag(max_iter, default = 0)))

logs2 <- logs2 %>%
  left_join(offsets, by = "run_id") %>%
  mutate(iter_cum = iter_rel + cum_offset)

ggplot(logs2, aes(x = iter_cum, y = Lh, colour = run)) +
  geom_line(alpha = 0.7) +
  theme_bw()

ggplot(logs, aes(x = Global.Trend, fill = run)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  theme_bw()

ggplot(logs, aes(x = Global.Trend)) +
  geom_density(alpha = 0.4, fill = "green") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  xlim(-0.05, 0.05) +
  labs(x = "Global trend", y = "Posterior density") +
  theme_bw() +
  theme(panel.grid = element_blank())

mean(logs$Global.Trend)
median(logs$Global.Trend)
quantile(logs$Global.Trend, c(0.025, 0.5, 0.975))
t.test(logs$Global.Trend, mu = 0)
100 * sum(logs$Global.Trend > 0) / length(logs$Global.Trend)

ggplot(logs, aes(x = Var, fill = run)) +
  geom_density(alpha = 0.4) +
  theme_bw()

ggplot(logs2, aes(x = iter_cum, y = Var, colour = run)) +
  geom_line(alpha = 0.7) +
  theme_bw()

mean(logs$Var)
median(logs$Var)
quantile(logs$Var, c(0.025, 0.5, 0.975))

ggplot(logs2, aes(x = iter_cum, y = Alpha, colour = run)) +
  geom_line(alpha = 0.7) +
  theme_bw()

beta_cols <- paste0("Beta.", 1:8)
beta_data <- logs[, beta_cols]
colnames(beta_data) <- c("arborescent", "columnar", "cushion_forming", "Epiphytic", "geophytic", "globose_caespitose", "globose_solitary", "shrubby")

beta_long <- melt(beta_data)
colnames(beta_long) <- c("Parameter", "Value")

crossing_summary <- aggregate(Value ~ Parameter, data = beta_long, function(x) {
  p_pos <- mean(x > 0, na.rm = TRUE)
  p_neg <- mean(x < 0, na.rm = TRUE)
  c(p_pos = p_pos, p_neg = p_neg, cross = min(p_pos, p_neg))
})

crossing_summary <- data.frame(Parameter = crossing_summary$Parameter,
  p_pos = crossing_summary$Value[, "p_pos"],
  p_neg = crossing_summary$Value[, "p_neg"],
  cross = crossing_summary$Value[, "cross"])

low_crossing <- subset(crossing_summary, cross < 0.05)
low_crossing

ggplot(beta_long, aes(x = Parameter, y = Value)) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal() +
  labs(x = NULL, y = "Value")

beta_mat <- as.data.frame(beta_data)

pairwise_beta_tests <- combn(colnames(beta_mat), 2, function(pair) {
  a <- pair[1]
  b <- pair[2]
  diff_ab <- beta_mat[[a]] - beta_mat[[b]]
  p_gt_0 <- mean(diff_ab > 0, na.rm = TRUE)
  p_lt_0 <- mean(diff_ab < 0, na.rm = TRUE)
  p_two_sided <- 2 * min(p_gt_0, p_lt_0)
  data.frame(beta_1 = a,
    beta_2 = b,
    mean_diff = mean(diff_ab, na.rm = TRUE),
    median_diff = median(diff_ab, na.rm = TRUE),
    lower_95 = quantile(diff_ab, 0.025, na.rm = TRUE),
    upper_95 = quantile(diff_ab, 0.975, na.rm = TRUE),
    p_gt_0 = p_gt_0,
    p_lt_0 = p_lt_0,
    p_two_sided = p_two_sided,
    significant = p_two_sided < 0.05)
}, simplify = FALSE) %>%
  bind_rows()

pairwise_beta_tests

write.csv(pairwise_beta_tests, "Cactus_height_Fabric_betas_differences.csv")

# Phylogenetic node effects
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

pp_df$evol_high <- ifelse(pp_df$Mean.Scalar > 1 &
                            pp_df$Sig.Scalar == 1 &
                            pp_df$Z.Scalar >= 2 &
                            pp_df$No.Sig.Nodes == 6, 1, 0)

pp_df$evol_low <- ifelse(pp_df$Mean.Scalar < 1 & pp_df$Sig.Scalar == 1 &
                            pp_df$Z.Scalar >= 2 & pp_df$No.Sig.Nodes == 6, 1, 0)

pp_df$direct_increase <- ifelse(pp_df$Mean..Beta...BL. > 0 &
                                  pp_df$Sig..Beta...BL. == 1 &
                                  pp_df$Z..Beta...BL. >= 2 &
                                  pp_df$No.Sig.Betas == 6, 1, 0)

pp_df$direct_decrease <- ifelse(pp_df$Mean..Beta...BL. < 0 & pp_df$Sig..Beta...BL. == 1 &
                                  pp_df$Z..Beta...BL. >= 2 & pp_df$No.Sig.Betas == 6, 1, 0)

pp_df$tradeoff_positive <- ifelse(pp_df$Beta.Scalar.Trade.off == "Retain Both - positively correlated" &
                                    pp_df$No.Sig.Nodes == 6, 1, 0)

pp_df$tradeoff_beta_only <- ifelse(pp_df$Beta.Scalar.Trade.off == "Retain Beta" &
                                     pp_df$No.Sig.Nodes == 6, 1, 0)

pp_df$tradeoff_indeterminate <- ifelse(pp_df$Beta.Scalar.Trade.off == "Retain Both - indeterminate" &
                                         pp_df$No.Sig.Nodes == 6, 1, 0)

pp_nodes <- pp_df[, c("ape_node", "Height", "evol_high", "evol_low", "direct_increase", "direct_decrease", "tradeoff_positive", "tradeoff_beta_only", "tradeoff_indeterminate")]

effect_summary <- pp_nodes %>%
  summarise(evol_high = sum(evol_high, na.rm = TRUE),
    evol_low = sum(evol_low, na.rm = TRUE),
    direct_increase = sum(direct_increase, na.rm = TRUE),
    direct_decrease = sum(direct_decrease, na.rm = TRUE),
    tradeoff_positive = sum(tradeoff_positive, na.rm = TRUE),
    tradeoff_beta_only = sum(tradeoff_beta_only, na.rm = TRUE),
    tradeoff_indeterminate = sum(tradeoff_indeterminate, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "effect", values_to = "n_nodes")

effect_summary

evol_high_nodes <- pp_nodes$ape_node[pp_nodes$evol_high == 1]
evol_low_nodes <- pp_nodes$ape_node[pp_nodes$evol_low == 1]
direct_increase_nodes <- pp_nodes$ape_node[pp_nodes$direct_increase == 1]
direct_decrease_nodes <- pp_nodes$ape_node[pp_nodes$direct_decrease == 1]
tradeoff_positive_nodes <- pp_nodes$ape_node[pp_nodes$tradeoff_positive == 1]
tradeoff_indeterminate_nodes <- pp_nodes$ape_node[pp_nodes$tradeoff_indeterminate == 1]

present_cats <- c(if (length(evol_high_nodes) > 0) "High Evolvability (Positive)",
  if (length(evol_low_nodes) > 0) "Low Evolvability (Negative)",
  if (length(direct_increase_nodes) > 0) "Positive Directionality",
  if (length(direct_decrease_nodes) > 0) "Negative Directionality",
  if (length(tradeoff_positive_nodes) > 0) "Trade-off (Positive Correlation)",
  if (length(tradeoff_indeterminate_nodes) > 0) "Trade-off (Indeterminate)")

shape_vals <- c("High Evolvability (Positive)" = 16, "Low Evolvability (Negative)" = 16,
                "Positive Directionality" = 17, "Negative Directionality" = 17,
                "Trade-off (Positive Correlation)" = 18, "Trade-off (Indeterminate)" = 18)

color_vals <- c("High Evolvability (Positive)" = "red", "Low Evolvability (Negative)" = "blue",
                "Positive Directionality" = "red", "Negative Directionality" = "blue",
                "Trade-off (Positive Correlation)" = "green", "Trade-off (Indeterminate)" = "orange")

q <- ggtree(tree, layout = "rectangular") +
  geom_point2(aes(subset = node %in% evol_high_nodes, color = "High Evolvability (Positive)",
                  shape = "High Evolvability (Positive)"), size = 2) +
  geom_point2(aes(subset = node %in% evol_low_nodes, color = "Low Evolvability (Negative)",
                  shape = "Low Evolvability (Negative)"), size = 2) +
  geom_point2(aes(subset = node %in% direct_increase_nodes, color = "Positive Directionality",
                  shape = "Positive Directionality"), size = 2) +
  geom_point2(aes(subset = node %in% direct_decrease_nodes, color = "Negative Directionality",
                  shape = "Negative Directionality"), size = 2) +
  geom_point2(aes(subset = node %in% tradeoff_positive_nodes, color = "Trade-off (Positive Correlation)",
                  shape = "Trade-off (Positive Correlation)"), size = 2) +
  geom_point2(aes(subset = node %in% tradeoff_indeterminate_nodes, color = "Trade-off (Indeterminate)",
                  shape = "Trade-off (Indeterminate)"), size = 2) +
  geom_fruit(
    data = cactus_height,
    geom = geom_bar,
    mapping = aes(y = species, x = plant_height, fill = plant_height),
    orientation = "y",
    stat = "identity",
    width = 0.6) +
  scale_color_manual(values = color_vals[present_cats],
    breaks = present_cats,
    name = "Node Categories" ) +
  scale_shape_manual(values = shape_vals[present_cats],
    breaks = present_cats,
    name = "Node Categories") +
  scale_fill_gradient(low = "blue",
    high = "red",
    name = "Plant height (log cm)") +
  guides(color = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      keywidth = 1,
      keyheight = 0.5,
      override.aes = list(
        shape = unname(shape_vals[present_cats]),
        size = rep(2, length(present_cats)),
        linetype = rep(0, length(present_cats)))),
    shape = "none",
    fill = guide_colorbar(order = 2,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 7,
      barheight = 0.5)) +
  theme_tree2() +
  theme(legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank())
q

# Effect sizes
scalar_effects <- pp_df %>%
  filter(Sig.Scalar == 1,
    Z.Scalar >= 2,
    No.Sig.Nodes == 6) %>%
  mutate(scalar_direction = case_when(
      Mean.Scalar > 1 ~ "higher_evolvability",
      Mean.Scalar < 1 ~ "lower_evolvability",
      TRUE ~ "neutral"),
    scalar_magnitude = Mean.Non.1.Scalar,
    scalar_distance_from_1 = abs(Mean.Non.1.Scalar - 1)) %>%
  dplyr::select(ape_node, ID, scalar_direction, scalar_magnitude,
                scalar_distance_from_1, Mean.Scalar, Mean.Non.1.Scalar, Median.Non.1,
                SD.Non.1.Scalar, Z.Scalar, No.Sig.Nodes)

beta_effects <- pp_df %>%
  filter(Sig..Beta...BL. == 1,
    Z..Beta...BL. >= 2,
    No.Sig.Betas == 6) %>%
  mutate(beta_direction = case_when(
      Mean..Beta...BL. > 0 ~ "positive_direction",
      Mean..Beta...BL. < 0 ~ "negative_direction",
      TRUE ~ "neutral"),
    beta_magnitude_signed = Mean..Beta...BL..NZ,
    beta_magnitude_abs = abs(Mean..Beta...BL..NZ)) %>%
  dplyr::select(ape_node, ID, beta_direction, beta_magnitude_signed,
                beta_magnitude_abs, Mean..Beta...BL., Mean..Beta...BL..NZ, SD..Beta...BL..NZ,
                Z..Beta...BL., No.Sig.Betas)

tradeoff_effects <- pp_df %>%
  filter(Beta.Scalar.Trade.off != "") %>%
  dplyr::select(ape_node, ID, Beta.Scalar.Trade.off, Mean.Non.1.Scalar,
                Mean..Beta...BL..NZ, Mean.Scalar, Mean..Beta...BL., Z.Scalar,
                Z..Beta...BL., No.Sig.Nodes, No.Sig.Betas) %>%
  mutate(scalar_distance_from_1 = abs(Mean.Non.1.Scalar - 1),
    beta_magnitude_abs = abs(Mean..Beta...BL..NZ))

ggplot(scalar_effects, aes(x = scalar_magnitude)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(x = "Scalar magnitude |Mean.Non.1.Scalar - 1|",
       y = "Count")

ggplot(beta_effects, aes(x = beta_magnitude_signed)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(x = "Directional magnitude |Mean..Beta...BL..NZ|",
       y = "Count")

# Branch rates and phylogeny
res1 <- readRDS("res1.rds")
res2 <- readRDS("res2.rds")
res3 <- readRDS("res3.rds")
res4 <- readRDS("res4.rds")
res5 <- readRDS("res5.rds")
res6 <- readRDS("res6.rds")

median_tree <- read.tree("Cactus_height_fabric_reg_gt_VR_scaled_median.tre")

plot(median_tree, show.tip.label = FALSE)
axisPhylo()

time_tree <- tree
runs <- list(res1, res2, res3, res4, res5, res6)

scalar_mat <- sapply(runs, function(r) r$ratesummary$medianscalar)

log_scalar_median <- apply(log(scalar_mat), 1, median)

edge_rate_df <- data.frame(node = as.integer(time_tree$edge[, 2]),
  log_scalar = log_scalar_median[-1])

p_bt <- ggtree(time_tree, layout = "rectangular") %<+% edge_rate_df +
  geom_tree(aes(color = log_scalar), linewidth = 0.6) +
  scale_color_viridis_c(option = "plasma") +
  theme_tree2() +
  theme(axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank())
p_bt

p <- ggtree(time_tree, layout = "rectangular") %<+% edge_rate_df +
  geom_tree(aes(color = log_scalar), linewidth = 0.6) +
  scale_color_viridis_c(option = "plasma", name = "log(rate)")
p <- p + ggnewscale::new_scale_color()
p <- p +
  geom_point2(aes(subset = node %in% direct_increase_nodes, color = "Positive Directionality",
                  shape = "Positive Directionality"), size = 2) +
  geom_point2(aes(subset = node %in% direct_decrease_nodes, color = "Negative Directionality",
                  shape = "Negative Directionality"), size = 2) +
  geom_point2(aes(subset = node %in% tradeoff_positive_nodes, color = "Trade-off (Positive Correlation)",
                  shape = "Trade-off (Positive Correlation)"), size = 2) +
  geom_point2(aes(subset = node %in% tradeoff_indeterminate_nodes, color = "Trade-off (Indeterminate)",
                  shape = "Trade-off (Indeterminate)"), size = 2)
p <- p +scale_color_manual(values = color_vals[present_cats],
    breaks = present_cats,
    name = "Node Categories") +
  scale_shape_manual(values = shape_vals[present_cats],
    breaks = present_cats,
    name = "Node Categories")
p <- p + geom_fruit(data = cactus_height,
    geom = geom_bar,
    mapping = aes(y = species, x = plant_height),
    orientation = "y",
    stat = "identity",
    width = 0.6,
    fill = "black")
p <- p + theme_tree2() +
  theme(legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank())
p

# Resolve beta-scalar trade-offs
pp_nodes_2 <- pp_df %>%
  mutate(evol_high = ifelse(Mean.Scalar > 1 & Sig.Scalar == 1 & Z.Scalar >= 2 & No.Sig.Nodes == 6, 1, 0),
    evol_low = ifelse(Mean.Scalar < 1 & Sig.Scalar == 1 & Z.Scalar >= 2 & No.Sig.Nodes == 6, 1, 0),
    direct_increase = ifelse(Mean..Beta...BL. > 0 & Sig..Beta...BL. == 1 & Z..Beta...BL. >= 2 &
                               No.Sig.Betas == 6, 1, 0),
    direct_decrease = ifelse(Mean..Beta...BL. < 0 & Sig..Beta...BL. == 1 & Z..Beta...BL. >= 2 &
                               No.Sig.Betas == 6, 1, 0),
    tradeoff_total = ifelse(!is.na(Beta.Scalar.Trade.off) & Beta.Scalar.Trade.off != "" &
                              No.Sig.Nodes == 6 & No.Sig.Betas == 6, 1, 0),
    beta_occ = `Observed.P..Beta...BL.`,
    scalar_occ = `Observed.P.Scalar`,
    N = `Sample.Size`,
    se_diff = sqrt((beta_occ * (1 - beta_occ) + scalar_occ * (1 - scalar_occ)) / N),
    se_diff = ifelse(is.na(se_diff) | se_diff == 0, NA, se_diff),
    z_diff = (beta_occ - scalar_occ) / se_diff,
    tradeoff_beta_only = ifelse(tradeoff_total == 1 & !is.na(z_diff) & z_diff >= 2, 1, 0),
    tradeoff_scalar_only = ifelse(tradeoff_total == 1 & !is.na(z_diff) & z_diff <= -2, 1, 0),
    tradeoff_indeterminate = ifelse(tradeoff_total == 1 & (is.na(z_diff) | abs(z_diff) < 2), 1, 0)) %>%
  dplyr::select(ape_node, Height, evol_high, evol_low, direct_increase, direct_decrease,
                tradeoff_total, tradeoff_beta_only, tradeoff_scalar_only, tradeoff_indeterminate,
                beta_occ, scalar_occ, N, se_diff, z_diff)

effect_summary_2 <- pp_nodes_2 %>%
  summarise(evol_high = sum(evol_high, na.rm = TRUE),
    evol_low = sum(evol_low, na.rm = TRUE),
    direct_increase = sum(direct_increase, na.rm = TRUE),
    direct_decrease = sum(direct_decrease, na.rm = TRUE),
    tradeoff_total = sum(tradeoff_total, na.rm = TRUE),
    tradeoff_beta_only = sum(tradeoff_beta_only, na.rm = TRUE),
    tradeoff_scalar_only = sum(tradeoff_scalar_only, na.rm = TRUE),
    tradeoff_indeterminate = sum(tradeoff_indeterminate, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "effect", values_to = "n_nodes")

effect_summary_2

evol_high_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$evol_high == 1]
evol_low_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$evol_low == 1]

direct_increase_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$direct_increase == 1]
direct_decrease_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$direct_decrease == 1]

tradeoff_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$tradeoff_total == 1]
tradeoff_beta_only_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$tradeoff_beta_only == 1]
tradeoff_scalar_only_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$tradeoff_scalar_only == 1]
tradeoff_indeterminate_nodes_2 <- pp_nodes_2$ape_node[pp_nodes_2$tradeoff_indeterminate == 1]

tradeoff_beta_lookup <- pp_df[pp_df$ape_node %in% tradeoff_beta_only_nodes_2, c("ape_node", "Mean..Beta...BL.")]

beta_only_to_increase <- tradeoff_beta_lookup$ape_node[tradeoff_beta_lookup$Mean..Beta...BL. > 0]
beta_only_to_decrease <- tradeoff_beta_lookup$ape_node[tradeoff_beta_lookup$Mean..Beta...BL. < 0]

tradeoff_scalar_lookup <- pp_df[pp_df$ape_node %in% tradeoff_scalar_only_nodes_2, c("ape_node", "Mean.Scalar")]

scalar_only_to_high <- tradeoff_scalar_lookup$ape_node[tradeoff_scalar_lookup$Mean.Scalar > 1]
scalar_only_to_low  <- tradeoff_scalar_lookup$ape_node[tradeoff_scalar_lookup$Mean.Scalar < 1]

evol_high_nodes_3 <- sort(unique(c(evol_high_nodes_2, scalar_only_to_high)))
evol_low_nodes_3  <- sort(unique(c(evol_low_nodes_2,  scalar_only_to_low)))

direct_increase_nodes_3 <- sort(unique(c(direct_increase_nodes_2, beta_only_to_increase)))
direct_decrease_nodes_3 <- sort(unique(c(direct_decrease_nodes_2, beta_only_to_decrease)))

effect_summary_3 <- data.frame(effect = c("evol_high", "evol_low", "direct_increase", "direct_decrease"),
  n_nodes = c(length(evol_high_nodes_3),
    length(evol_low_nodes_3),
    length(direct_increase_nodes_3),
    length(direct_decrease_nodes_3)))

effect_summary_3

print(list(tradeoff_beta_only_nodes_2 = tradeoff_beta_only_nodes_2,
  beta_only_to_increase = beta_only_to_increase,
  beta_only_to_decrease = beta_only_to_decrease,
  tradeoff_scalar_only_nodes_2 = tradeoff_scalar_only_nodes_2,
  scalar_only_to_high = scalar_only_to_high,
  scalar_only_to_low = scalar_only_to_low))

scalar_effects
beta_effects

write.csv(scalar_effects, "Cactus_Fabric_reg_GT_scalar_effects.csv")
write.csv(beta_effects, "Cactus_Fabric_reg_GT_beta_effects.csv")

saveRDS(evol_high_nodes_3, "Cactus_Fabric_reg_GT_evol_high_nodes_3.rds")
saveRDS(evol_low_nodes_3, "Cactus_Fabric_reg_GT_evol_low_nodes_3.rds")
saveRDS(direct_increase_nodes_3, "Cactus_Fabric_reg_GT_direct_increase_nodes_3.rds")
saveRDS(direct_decrease_nodes_3, "Cactus_Fabric_reg_GT_direct_decrease_nodes_3.rds")
saveRDS(time_tree, "Cactus_Fabric_reg_GT_time_tree.rds")
saveRDS(edge_rate_df, "Cactus_Fabric_reg_GT_edge_rate_df.rds")
saveRDS(cactus_height, "Cactus_Fabric_reg_GT_cactus_height.rds")

# Final tree
present_cats_3 <- c(if (length(evol_high_nodes_3) > 0) "High Evolvability (Positive)",
  if (length(evol_low_nodes_3) > 0) "Low Evolvability (Negative)",
  if (length(direct_increase_nodes_3) > 0) "Positive Directionality",
  if (length(direct_decrease_nodes_3) > 0) "Negative Directionality")

shape_vals_3 <- c("High Evolvability (Positive)" = 16, "Low Evolvability (Negative)" = 16,
                  "Positive Directionality" = 17, "Negative Directionality" = 17)

color_vals_3 <- c("High Evolvability (Positive)" = "green", "Low Evolvability (Negative)" = "blue",
                  "Positive Directionality" = "green", "Negative Directionality" = "blue")

tree_plot <- ggtree(time_tree, layout = "rectangular") %<+% edge_rate_df +
  geom_tree(aes(color = log_scalar), linewidth = 0.6) +
  scale_color_viridis_c(option = "plasma", name = "log(rate)")

tree_plot <- tree_plot + ggnewscale::new_scale_color()

tree_plot <- tree_plot + geom_point2(aes(subset = node %in% evol_high_nodes_3, color = "High Evolvability (Positive)",
                  shape = "High Evolvability (Positive)"), size = 2) +
  geom_point2(aes(subset = node %in% evol_low_nodes_3, color = "Low Evolvability (Negative)",
                  shape = "Low Evolvability (Negative)"), size = 2) +
  geom_point2(aes(subset = node %in% direct_increase_nodes_3, color = "Positive Directionality",
                  shape = "Positive Directionality"), size = 2) +
  geom_point2(aes(subset = node %in% direct_decrease_nodes_3, color = "Negative Directionality",
                  shape = "Negative Directionality"), size = 2)

tree_plot <- tree_plot +
  scale_color_manual(values = color_vals_3[present_cats_3],
    breaks = present_cats_3,
    name = "Node Categories") +
  scale_shape_manual(
    values = shape_vals_3[present_cats_3],
    breaks = present_cats_3,
    name = "Node Categories")

tree_plot <- tree_plot +
  geom_fruit(data = cactus_height,
    geom = geom_bar,
    mapping = aes(y = species, x = plant_height),
    orientation = "y",
    stat = "identity",
    width = 0.6,
    fill = "black")

tree_plot <- tree_plot +
  theme_tree2() +
  theme(legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank())

tree_plot