# Cactus height: phylogenetic signal and descriptive summaries

library(ape)
library(dplyr)
library(geiger)
library(phytools)

setwd("C:/Users/hy929145/Documents/Cactus height")

# Phylogenetic signal
tree <- read.tree("cacti_treepl_THH.tre")
height <- read.table("Cactus_height_no_gf_BT_sample_July_correction_final_prue.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
tree <- keep.tip(tree, height$V1)
phylosig(tree, setNames(height$V2, height$V1), method = "lambda", test = TRUE)

# Height summaries by taxonomic group
colnames(height) <- c("species", "height")
height$height <- as.numeric(height$height)

node_definitions <- list("Cactoideae" = c("Echinopsis_calorubra", "Blossfeldia_liliputana"),
  "Opuntioideae" = c("Opuntia_lasiacantha", "Pereskiopsis_rotundifolia"),
  "Maihuenioideae" = c("Maihuenia_poeppigii", "Maihuenia_patagonica"),
  "Cylindropuntieae" = c("Cylindropuntia_prolifera", "Pereskiopsis_rotundifolia"),
  "Opuntieae" = c("Opuntia_basilaris", "Opuntia_chlorotica"),
  "Pterocacteae" = c("Cumulopuntia_chichensis", "Pterocactus_valentinii"),
  "Rhipsalideae" = c("Rhipsalis_cuneata", "Hatiora_cylindrica"),
  "Cacteae" = c("Mammillaria_albilanata", "Geohintonia_mexicana"),
  "Cereeae" = c("Uebelmannia_pectinifera", "Echinopsis_calorubra"),
  "Notocacteae" = c("Eriosyce_recondita", "Parodia_haselbergii"),
  "Phyllocacteae" = c("Pfeiffera_paranganiensis", "Echinocereus_bonkerae"),
  "Copiapoeae" = c("Copiapoa_marginata", "Copiapoa_dealbata"))

node_definitions[["Pereskia"]] <- grep("^Pereskia_", tree$tip.label, value = TRUE)
node_definitions[["Leuenbergeria"]] <- grep("^Leuenbergeria_", tree$tip.label, value = TRUE)
node_definitions[["Fraileeae"]] <- grep("^Frailea_", tree$tip.label, value = TRUE)

taxonomy <- data.frame(species = tree$tip.label, stringsAsFactors = FALSE)
taxonomy$genus <- sub("_.*", "", taxonomy$species)
taxonomy$tribe <- NA_character_
taxonomy$subfamily <- NA_character_

for (clade in names(node_definitions)) {
  defining_species <- intersect(node_definitions[[clade]], tree$tip.label)
  if (length(defining_species) < 2) next
  node <- ape::getMRCA(tree, defining_species)
  spp <- geiger::tips(tree, node)
  if (grepl("oideae$", clade, ignore.case = TRUE)) {
    taxonomy$subfamily[taxonomy$species %in% spp] <- clade
  }
  if (grepl("eae$", clade, ignore.case = TRUE) && !grepl("oideae$", clade, ignore.case = TRUE)) {
    taxonomy$tribe[taxonomy$species %in% spp] <- clade
  }
}

taxonomy$subfamily[taxonomy$genus %in% c("Pereskia", "Leuenbergeria")] <- "Pereskioideae"
taxonomy$tribe[taxonomy$genus %in% c("Pereskia")] <- "Pereskia"
taxonomy$tribe[taxonomy$genus %in% c("Leuenbergeria")] <- "Leuenbergeria"
taxonomy$tribe[taxonomy$genus %in% c("Blossfeldia")] <- "Blossfeldieae"
taxonomy$tribe[taxonomy$genus %in% c("Frailea")] <- "Fraileeae"
taxonomy$tribe[taxonomy$genus %in% c("Maihuenia")] <- "Maihuenia"

taxonomy <- left_join(taxonomy, height)
taxonomy$raw_height <- exp(taxonomy$height)

overall_height_summary <- data.frame(
  n_species = sum(!is.na(taxonomy$raw_height)),
  mean = mean(taxonomy$raw_height, na.rm = TRUE),
  median = median(taxonomy$raw_height, na.rm = TRUE),
  sd = sd(taxonomy$raw_height, na.rm = TRUE),
  min = min(taxonomy$raw_height, na.rm = TRUE),
  max = max(taxonomy$raw_height, na.rm = TRUE)
)
overall_height_summary

subfamily_summary <- taxonomy %>%
  group_by(subfamily) %>%
  summarise(n_species = sum(!is.na(raw_height)), mean_height = mean(raw_height, na.rm = TRUE),
            median_height = median(raw_height, na.rm = TRUE), sd_height = sd(raw_height, na.rm = TRUE),
            min_height = min(raw_height, na.rm = TRUE), max_height = max(raw_height, na.rm = TRUE), .groups = "drop")
subfamily_summary

tribe_summary <- taxonomy %>%
  group_by(tribe) %>%
  summarise(n_species = sum(!is.na(raw_height)), mean_height = mean(raw_height, na.rm = TRUE),
            median_height = median(raw_height, na.rm = TRUE), sd_height = sd(raw_height, na.rm = TRUE),
            min_height = min(raw_height, na.rm = TRUE), max_height = max(raw_height, na.rm = TRUE), .groups = "drop")
tribe_summary

# Growth-form summaries
height_gf <- read.table("Cactus_height_gf_BT_sample_July_correction_final_prue.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
height_gf$V2 <- NULL
colnames(height_gf) <- c("species", "arborescent", "columnar", "cushion_forming", "Epiphytic",
                         "geophytic", "globose_caespitose", "globose_solitary", "shrubby")
head(height_gf)

growth_cols <- c("arborescent", "columnar", "cushion_forming", "Epiphytic",
                 "geophytic", "globose_caespitose", "globose_solitary", "shrubby")
growth_summary <- data.frame(growth_form = growth_cols,
                             n_species = colSums(height_gf[, growth_cols]),
                             proportion = colMeans(height_gf[, growth_cols]))
growth_summary

height_gf$n_growth_forms <- rowSums(height_gf[, growth_cols])
table(height_gf$n_growth_forms)

taxonomy_gf <- left_join(taxonomy, height_gf)
subfamily_growth <- taxonomy_gf %>%
  group_by(subfamily) %>%
  summarise(n_species = n(), across(all_of(growth_cols), sum, na.rm = TRUE))
subfamily_growth
