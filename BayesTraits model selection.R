library(coda)
library(tidyverse)

setwd("C:/Users/hy929145/Documents/Cactus height")

# Models
model_folders <- c("cactus_height_reg_may",
  "cactus_height_univ_fabric_gt_may",
  "cactus_height_vr_reg_global_trend_may",
  "cactus_height_vr_July",
  "cactus_height_univ_fabric_may",
  "cactus_height_Fabric_reg_July_gt",
  "cactus_height_Fabric_reg_July")

# Model labels reflect the analyses, rather than the original folder names (minor error in naming)
model_names <- c("VR regression",
  "Univariate Fabric without global trend",
  "VR regression with global trend",
  "VR without regression",
  "Univariate Fabric with global trend",
  "Fabric regression with global trend",
  "Fabric regression without global trend")

# Read BayesTraits logs
log_lists <- setNames(vector("list", length(model_folders)), model_names)

for (i in seq_along(model_folders)) {
  logs <- vector("list", length(log_numbers <- 0:5))
  names(logs) <- paste0("run_", log_numbers)

  for (j in seq_along(log_numbers)) {
    file <- file.path(model_folders[i], sprintf("Data-%07d.txt.Log.txt", log_numbers[j]))
    lines <- readLines(file, warn = FALSE)
    header_position <- grep("Iteration\tLh\tTree No", lines, fixed = TRUE)
    if (length(header_position) == 0) {
      stop("Couldn't find table header in: ", file)
    }
    header_position <- header_position[1]
    header_line <- sub("\t+$", "", lines[header_position])
    data_lines <- lines[(header_position + 1):length(lines)]
    data_lines <- sub("\t+$", "", data_lines)
    data_lines <- data_lines[nzchar(data_lines)]
    column_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
    number_data_fields <- length(strsplit(data_lines[1], "\t", fixed = TRUE)[[1]])
    number_header_fields <- length(column_names)
    if (number_data_fields > number_header_fields) {
      number_missing_names <- number_data_fields - number_header_fields
      if (number_missing_names == 1) {
        column_names <- c(column_names, "Global Trend")
      } else {
        column_names <- c(column_names, paste0("Extra_", seq_len(number_missing_names)))
      }
    }
    logs[[j]] <- read.delim(text = paste(data_lines, collapse = "\n"), sep = "\t",
                            header = FALSE, col.names = column_names, row.names = NULL,
                            stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE,
                            quote = "", comment.char = "")
    names(logs[[j]])[names(logs[[j]]) == "Tree No"] <- "Tree_No"
    logs[[j]]$Iteration <- as.numeric(logs[[j]]$Iteration)
    logs[[j]]$Lh <- as.numeric(logs[[j]]$Lh)
  }
  log_lists[[i]] <- logs
}

# Combine runs and plot likelihood traces
combined_logs <- list()
for (model in names(log_lists)) {
  combined <- do.call(rbind, lapply(names(log_lists[[model]]), function(run) {
    x <- log_lists[[model]][[run]]
    x$Run <- run
    x
  }))
  combined$cumulative_itt <- seq_len(nrow(combined))
  combined_logs[[model]] <- combined
}

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (model in names(combined_logs)) {
  plot(combined_logs[[model]]$cumulative_itt, combined_logs[[model]]$Lh,
       type = "l", col = "steelblue", xlab = "MCMC sample",
       ylab = "Log likelihood (Lh)", main = model)
}
par(mfrow = c(1, 1))

# Stepping-stone comparisons
stones_results <- data.frame(Model = character(), Run = character(), Log_Marginal_Likelihood = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(model_folders)) {
  model <- model_names[i]
  for (j in log_numbers) {
    file <- file.path(model_folders[i], sprintf("Data-%07d.txt.Stones.txt", j))
    lines <- readLines(file)
    ml_line <- grep("Log marginal likelihood", lines, value = TRUE)
    if (length(ml_line) == 0) {
      warning("No marginal likelihood found in: ", file)
      next
    }
    ml <- as.numeric(sub(".*\\t", "", ml_line[length(ml_line)]))
    stones_results <- rbind(stones_results, data.frame(Model = model, Run = paste0("run_", j), Log_Marginal_Likelihood = ml, stringsAsFactors = FALSE))
  }
}

stones_results

stones_wide <- stones_results %>% pivot_wider(names_from = Run, values_from = Log_Marginal_Likelihood) %>%
  mutate(Mean = rowMeans(dplyr::select(., starts_with("run_")), na.rm = TRUE), Median = apply(dplyr::select(., starts_with("run_")), 1, median, na.rm = TRUE)) %>% arrange(desc(Mean))

stones_wide
write.csv(stones_wide, "Cactus_height_BayesTraits_stepping_stones.csv")

stones_plot <- ggplot(stones_results, aes(x = reorder(Model, Log_Marginal_Likelihood, FUN = median), y = Log_Marginal_Likelihood)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_point(aes(shape = Run), size = 3, position = position_jitter(width = 0.08)) +
  coord_flip() + theme_bw() +
  labs(x = NULL, y = "Log marginal likelihood", shape = "Run")

stones_plot

# R2 and global trend for the Fabric regression
median(combined_logs$`Fabric regression with global trend`$`R^2`)

global_trend <- combined_logs$`Fabric regression with global trend`$`Global Trend`
median(global_trend)
HPDinterval(as.mcmc(global_trend), prob = 0.95)

n_negative <- sum(global_trend < 0)
n_positive <- sum(global_trend > 0)
prop_negative <- mean(global_trend < 0)
prop_positive <- mean(global_trend > 0)

data.frame(negative_samples = n_negative, positive_samples = n_positive,proportion_negative = prop_negative, proportion_positive = prop_positive)