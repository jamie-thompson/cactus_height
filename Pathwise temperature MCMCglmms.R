# Cactus height and pathwise temperature MCMCglmm analyses

library(ape)
library(coda)
library(dplyr)
library(MCMCglmm)
library(phytools)

setwd("C:/Users/hy929145/Documents/Cactus height")

tree <- read.tree("cacti_treepl_THH.tre")
cedb <- read.table("Cactus_height_gf_BT_sample_July_correction_final_prue.txt", sep = "\t")
colnames(cedb) <- c("species", "log_height", "arborescent", "columnar", "cushion_forming", "Epiphytic", "geophytic", "globose_caespitose", "globose_solitary", "shrubby")

pw_beta <- read.csv("Cacti_pathwise_beta_July_log1p.csv")
pw_clim <- read.csv("Cacti_pathwise_climate_changes_July_log1p.csv")
extant_temp_prec <- read.csv("Cacti_extant_temp_prec_geo_Jorge.csv")[c("species", "median_temperature", "median_precipitation")]
lat <- read.csv("CactEcoDB_latitude_summary_stats.csv")[c("species", "median_lat")]
lat$abs_lat <- abs(lat$median_lat)

pw_beta_clim <- pw_beta %>%
  inner_join(pw_clim, by = "species") %>%
  inner_join(cedb, by = "species") %>%
  inner_join(extant_temp_prec, by = "species") %>%
  inner_join(lat, by = "species")

species_keep <- intersect(tree$tip.label, pw_beta_clim$species)
tree_mcmc <- drop.tip(tree, setdiff(tree$tip.label, species_keep))
dat_mcmc <- pw_beta_clim[pw_beta_clim$species %in% species_keep, ]
dat_mcmc$species <- factor(dat_mcmc$species)

# Separate complete-case datasets are retained for the original model fits.
dat_mv <- dat_mcmc[complete.cases(dat_mcmc[, c("species", "log_dir_pos", "log_dir_neg", "log_temp_pos", "log_temp_neg", "median_temperature")]), ]
dat_pos_neg <- dat_mcmc[complete.cases(dat_mcmc[, c("species", "log_height", "log_temp_pos", "log_temp_neg", "median_temperature")]), ]
dat_abs <- dat_mcmc[complete.cases(dat_mcmc[, c("species", "log_height", "log_temp_abs", "median_temperature")]), ]

tree_mv <- drop.tip(tree_mcmc, setdiff(tree_mcmc$tip.label, dat_mv$species))
tree_pos_neg <- drop.tip(tree_mcmc, setdiff(tree_mcmc$tip.label, dat_pos_neg$species))
tree_abs <- drop.tip(tree_mcmc, setdiff(tree_mcmc$tip.label, dat_abs$species))

dat_mv$species <- droplevels(dat_mv$species)
all(tree_mv$tip.label %in% dat_mv$species)
all(dat_mv$species %in% tree_mv$tip.label)

tree_mv <- force.ultrametric(tree_mv, method = "extend")
tree_pos_neg <- force.ultrametric(tree_pos_neg, method = "extend")
tree_abs <- force.ultrametric(tree_abs, method = "extend")

Ainv_mv <- inverseA(tree_mv, nodes = "ALL", scale = TRUE)$Ainv
Ainv_pos_neg <- inverseA(tree_pos_neg, nodes = "ALL", scale = TRUE)$Ainv
Ainv_abs <- inverseA(tree_abs, nodes = "ALL", scale = TRUE)$Ainv

dat_mv$species <- factor(dat_mv$species, levels = rownames(Ainv_mv))
dat_pos_neg$species <- factor(dat_pos_neg$species, levels = rownames(Ainv_pos_neg))
dat_abs$species <- factor(dat_abs$species, levels = rownames(Ainv_abs))

# Multivariate models: positive and negative directional histories
prior_mv <- list(R = list(V = diag(2), nu = 2), G = list(G1 = list(V = diag(2), nu = 2)))

mod_mv_no_ext <- MCMCglmm(cbind(log_dir_pos, log_dir_neg) ~ trait - 1 + trait:log_temp_pos + trait:log_temp_neg + trait:log(median_temperature),
  random = ~ us(trait):species, rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"), ginverse = list(species = Ainv_mv),
  data = dat_mv, prior = prior_mv,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_mv_no_ext, "MCMCglmm_multivariate_dir_temp_5m.rds")
summary(mod_mv_no_ext$Sol)
summary(mod_mv_no_ext$VCV)
summary(mod_mv_no_ext)

mod_mv_no_current <- MCMCglmm(cbind(log_dir_pos, log_dir_neg) ~ trait - 1 + trait:log_temp_pos + trait:log_temp_neg,
  random = ~ us(trait):species, rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"), ginverse = list(species = Ainv_mv),
  data = dat_mv, prior = prior_mv,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_mv_no_current, "MCMCglmm_multivariate_dir_temp_no_current_5m.rds")

mod_mv_interaction <- MCMCglmm(cbind(log_dir_pos, log_dir_neg) ~ trait - 1 + trait:log_temp_pos + trait:log_temp_neg + trait:log_temp_pos:log_temp_neg,
  random = ~ us(trait):species, rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"), ginverse = list(species = Ainv_mv),
  data = dat_mv, prior = prior_mv,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_mv_interaction, "MCMCglmm_multivariate_dir_temp_no_current_INTERACTIONS_5m_RIGHT.rds")

# Alternative multivariate prior
prior_mv_alt1 <- list(R = list(V = diag(2) * (2/3), nu = 3), G = list(G1 = list(V = diag(2) * (2/3), nu = 3)))

mod_mv_no_ext_prior_mv_alt1 <- MCMCglmm(cbind(log_dir_pos, log_dir_neg) ~ trait - 1 + trait:log_temp_pos + trait:log_temp_neg + trait:log(median_temperature),
  random = ~ us(trait):species, rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"), ginverse = list(species = Ainv_mv),
  data = dat_mv, prior = prior_mv_alt1,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_mv_no_ext_prior_mv_alt1, "MCMCglmm_multivariate_dir_temp_5m_prior_mv_alt1.rds")
summary(mod_mv_no_ext_prior_mv_alt1$Sol)
summary(mod_mv_no_ext_prior_mv_alt1$VCV)
summary(mod_mv_no_ext_prior_mv_alt1)

mod_mv_no_current_prior_mv_alt1 <- MCMCglmm(cbind(log_dir_pos, log_dir_neg) ~ trait - 1 + trait:log_temp_pos + trait:log_temp_neg,
  random = ~ us(trait):species, rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"), ginverse = list(species = Ainv_mv),
  data = dat_mv, prior = prior_mv_alt1,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_mv_no_current_prior_mv_alt1, "MCMCglmm_multivariate_dir_temp_no_current_5m_prior_mv_alt1.rds")
summary(mod_mv_no_current_prior_mv_alt1$Sol)
summary(mod_mv_no_current_prior_mv_alt1$VCV)
summary(mod_mv_no_current_prior_mv_alt1)

mod_mv_interaction_prior_mv_alt1 <- MCMCglmm(cbind(log_dir_pos, log_dir_neg) ~ trait - 1 + trait:log_temp_pos + trait:log_temp_neg + trait:log_temp_pos:log_temp_neg,
  random = ~ us(trait):species, rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"), ginverse = list(species = Ainv_mv),
  data = dat_mv, prior = prior_mv_alt1,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_mv_interaction_prior_mv_alt1, "MCMCglmm_multivariate_dir_temp_no_current_INTERACTIONS_5m_RIGHT_prior_mv_alt1.rds")
summary(mod_mv_interaction_prior_mv_alt1$Sol)
summary(mod_mv_interaction_prior_mv_alt1$VCV)
summary(mod_mv_interaction_prior_mv_alt1)

cor.test(dat_mv$log_dir_pos, dat_mv$log_dir_neg)

# Present-day cactus height models
prior_uni <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))

mod_height_pos_neg <- MCMCglmm(log_height ~ log_temp_pos + log_temp_neg + log(median_temperature),
  random = ~ species, rcov = ~ units, family = "gaussian",
  ginverse = list(species = Ainv_pos_neg), data = dat_pos_neg, prior = prior_uni,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_height_pos_neg, "MCMCglmm_height_temp_pos_neg_control_extant_temp.rds")
mod_height_pos_neg <- readRDS("MCMCglmm_height_temp_pos_neg_control_extant_temp.rds")

mod_height_abs <- MCMCglmm(log_height ~ log_temp_abs + log(median_temperature),
  random = ~ species, rcov = ~ units, family = "gaussian",
  ginverse = list(species = Ainv_abs), data = dat_abs, prior = prior_uni,
  nitt = 5500000, burnin = 500000, thin = 5000, verbose = TRUE)
saveRDS(mod_height_abs, "MCMCglmm_height_temp_abs_control_extant_temp.rds")
mod_height_abs <- readRDS("MCMCglmm_height_temp_abs_control_extant_temp.rds")

summary(mod_height_pos_neg)
summary(mod_height_pos_neg$Sol)
summary(mod_height_pos_neg$VCV)
summary(mod_height_abs)
summary(mod_height_abs$Sol)
summary(mod_height_abs$VCV)