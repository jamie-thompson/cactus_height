# Cactus height: precipitation, temperature and latitude MCMCglmm models

library(ape)
library(dplyr)
library(MCMCglmm)
library(phytools)

setwd("C:/Users/hy929145/Documents/Cactus height")

tree <- read.tree("cacti_treepl_THH.tre")
cedb <- read.table("Cactus_height_gf_BT_sample_July_correction_final_prue.txt", sep = "\t")
colnames(cedb) <- c("species", "log_height", "arborescent", "columnar", "cushion_forming",
                    "Epiphytic", "geophytic", "globose_caespitose", "globose_solitary", "shrubby")

extant_temp_prec <- read.csv("Cacti_extant_temp_prec_geo_Jorge.csv")[c("species", "median_temperature", "median_precipitation")]
lat <- read.csv("CactEcoDB_latitude_summary_stats.csv")[c("species", "median_lat")]
lat$abs_lat <- abs(lat$median_lat)

growth_cols <- c("arborescent", "columnar", "cushion_forming", "Epiphytic", "geophytic",
                 "globose_caespitose", "globose_solitary", "shrubby")

dat <- cedb %>%
  inner_join(extant_temp_prec, by = "species") %>%
  inner_join(lat, by = "species")

# One species set and one phylogeny for all three analyses
model_vars <- c("species", "log_height", "median_precipitation", "median_temperature", "abs_lat", growth_cols)
dat_mcmc <- dat[complete.cases(dat[, model_vars]), ]
species_keep <- intersect(tree$tip.label, dat_mcmc$species)
dat_mcmc <- dat_mcmc[dat_mcmc$species %in% species_keep, ]
tree_mcmc <- drop.tip(tree, setdiff(tree$tip.label, species_keep))
tree_mcmc <- force.ultrametric(tree_mcmc, method = "extend")
Ainv <- inverseA(tree_mcmc, nodes = "ALL", scale = TRUE)$Ainv
dat_mcmc$species <- factor(dat_mcmc$species, levels = rownames(Ainv))

prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))

# Precipitation
m0_prec <- MCMCglmm(log_height ~ 1,
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m0_prec, "MCMCglmm_height_precipitation_growthform_m0_null_5m.rds")

m1_prec <- MCMCglmm(log_height ~ log(median_precipitation),
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m1_prec, "MCMCglmm_height_precipitation_growthform_m1_precipitation_5m.rds")

m2_prec <- MCMCglmm(log_height ~ log(median_precipitation) + arborescent + columnar + cushion_forming +
                      Epiphytic + geophytic + globose_caespitose + globose_solitary + shrubby,
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m2_prec, "MCMCglmm_height_precipitation_growthform_m2_additive_5m.rds")

m3_prec <- MCMCglmm(log_height ~ log(median_precipitation) * (arborescent + columnar + cushion_forming +
                       Epiphytic + geophytic + globose_caespitose + globose_solitary + shrubby),
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m3_prec, "MCMCglmm_height_precipitation_growthform_m3_interactions_5m.rds")

# Temperature
m0_temp <- MCMCglmm(log_height ~ 1,
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m0_temp, "MCMCglmm_height_temperature_growthform_m0_null_5m.rds")

m1_temp <- MCMCglmm(log_height ~ log(median_temperature),
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m1_temp, "MCMCglmm_height_temperature_growthform_m1_temperature_5m.rds")

m2_temp <- MCMCglmm(log_height ~ log(median_temperature) + arborescent + columnar + cushion_forming +
                      Epiphytic + geophytic + globose_caespitose + globose_solitary + shrubby,
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m2_temp, "MCMCglmm_height_temperature_growthform_m2_additive_5m.rds")

m3_temp <- MCMCglmm(log_height ~ log(median_temperature) * (arborescent + columnar + cushion_forming +
                       Epiphytic + geophytic + globose_caespitose + globose_solitary + shrubby),
                    random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                    data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                    thin = 5000, verbose = TRUE)
saveRDS(m3_temp, "MCMCglmm_height_temperature_growthform_m3_interactions_5m.rds")

# Latitude
m0_lat <- MCMCglmm(log_height ~ 1,
                   random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                   data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                   thin = 5000, verbose = TRUE)
saveRDS(m0_lat, "MCMCglmm_height_latitude_growthform_m0_null_5m.rds")

m1_lat <- MCMCglmm(log_height ~ log(abs_lat),
                   random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                   data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                   thin = 5000, verbose = TRUE)
saveRDS(m1_lat, "MCMCglmm_height_LOG_latitude_growthform_m1_latitude_5m.rds")

m2_lat <- MCMCglmm(log_height ~ log(abs_lat) + arborescent + columnar + cushion_forming +
                     Epiphytic + geophytic + globose_caespitose + globose_solitary + shrubby,
                   random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                   data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                   thin = 5000, verbose = TRUE)
saveRDS(m2_lat, "MCMCglmm_height_LOG_latitude_growthform_m2_additive_5m.rds")

m3_lat <- MCMCglmm(log_height ~ log(abs_lat) * (arborescent + columnar + cushion_forming +
                      Epiphytic + geophytic + globose_caespitose + globose_solitary + shrubby),
                   random = ~ species, family = "gaussian", ginverse = list(species = Ainv),
                   data = dat_mcmc, prior = prior, nitt = 5500000, burnin = 500000,
                   thin = 5000, verbose = TRUE)
saveRDS(m3_lat, "MCMCglmm_height_LOG_latitude_growthform_m3_interactions_5m.rds")
