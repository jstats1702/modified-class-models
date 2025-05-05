# Clean environment and set working directory
rm(list = ls())

# Load packages
suppressMessages({
     library(igraph)
     library(corrplot)
     library(mclust)
     library(ggplot2)
     library(dplyr)
     library(mcclust.ext)  # Install with: devtools::install_github("muschellij2/mcclust.ext")
})

# Define and set base path
base_path <- "C:/Users/User/Dropbox/PAPERS/projects/eleni/salter/"
setwd(base_path)

dataset <- "salter"
n_thin  <- NA

# Load helper functions
source(paste0(base_path, "helper_functions_data_analysis.R"))

# ==============================================================================
# Get results ==================================================================

# CLASS

set.seed(123)

results_salter_class_dm <- get_results_dm(base_path,
                                          model   = "class", 
                                          prior   = "dm",
                                          dataset,
                                          n_thin)

results_salter_class_dp <- get_results_np(base_path,
                                          model   = "class",
                                          prior   = "dp",
                                          dataset,
                                          n_thin)

results_salter_class_py <- get_results_np(base_path,
                                          model   = "class",
                                          prior   = "py",
                                          dataset,
                                          n_thin)

results_salter_class_gn <- get_results_np(base_path,
                                          model   = "class",
                                          prior   = "gn",
                                          dataset,
                                          n_thin)

# CLIST

results_salter_clist_dm <- get_results_dm(base_path,
                                          model   = "clist",
                                          prior   = "dm",
                                          dataset,
                                          n_thin)

results_salter_clist_dp <- get_results_np(base_path,
                                          model   = "clist",
                                          prior   = "dp",
                                          dataset,
                                          n_thin)

results_salter_clist_py <- get_results_np(base_path,
                                          model   = "clist",
                                          prior   = "py",
                                          dataset,
                                          n_thin)

results_salter_clist_gn <- get_results_np(base_path,
                                          model   = "clist",
                                          prior   = "gn",
                                          dataset,
                                          n_thin)

# CLEIGEN

results_salter_cleigen_dm <- get_results_dm(base_path,
                                            model   = "cleigen",
                                            prior   = "dm",
                                            dataset,
                                            n_thin)

results_salter_cleigen_dp <- get_results_np(base_path,
                                            model   = "cleigen",
                                            prior   = "dp",
                                            dataset,
                                            n_thin)

results_salter_cleigen_py <- get_results_np(base_path,
                                            model   = "cleigen",
                                            prior   = "py",
                                            dataset,
                                            n_thin)

results_salter_cleigen_gn <- get_results_np(base_path,
                                            model   = "cleigen",
                                            prior   = "gn",
                                            dataset,
                                            n_thin)

# save(results_salter_class_dm,
#      results_salter_class_dp,
#      results_salter_class_py,
#      results_salter_class_gn,
#      results_salter_clist_dm,
#      results_salter_clist_dp,
#      results_salter_clist_py,
#      results_salter_clist_gn,
#      results_salter_cleigen_dm,
#      results_salter_cleigen_dp,
#      results_salter_cleigen_py,
#      results_salter_cleigen_gn, 
#      file = "results_salter.RData")

load("results_salter.RData")

# ==============================================================================
# Model/prior comparisson using ARI ============================================

# Define models and priors
models <- c("class", "clist", "cleigen")
priors <- c("dm", "dp", "py", "gn")

# Build all combination names
names_list <- expand.grid(model = models, prior = priors)
names_list <- apply(names_list, 1, function(x) paste(x, collapse = "_"))

# Extract all label vectors
labels_list <- lapply(names_list, function(name) {
     eval(parse(text = paste0("results_salter_", name, "$sorted_network_data$labels")))
})

names(labels_list) <- names_list

# Initialize ARI matrix
n <- length(labels_list)
ari_matrix <- matrix(NA, nrow = n, ncol = n)
rownames(ari_matrix) <- names_list
colnames(ari_matrix) <- names_list

# Fill the ARI matrix
for (i in 1:n) {
     for (j in 1:n) {
          ari_matrix[i, j] <- adjustedRandIndex(labels_list[[i]], labels_list[[j]])
     }
}

# View table
print(round(ari_matrix, 3))

# Plot upper triangle with black labels and no diagonal
pdf(file = paste0(base_path, "ari_comparisson.pdf"))
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

corrplot(ari_matrix,
         method = "color",
         type = "upper",
         col = colorRampPalette(c("white", "red"))(100),
         tl.col = "black",
         tl.cex = 0.8,
         addCoef.col = "black",
         number.cex = 0.6,
         diag = FALSE,
         is.corr = FALSE)

dev.off()

# ==============================================================================
# Number of clusters distribution ==============================================

# Define model–prior combinations
models <- c("class", "clist", "cleigen")
priors <- c("dm", "dp", "py", "gn")
names_list <- expand.grid(model = models, prior = priors)
names_list <- apply(names_list, 1, function(x) paste(x, collapse = "_"))

# Summarize posterior number of clusters
cluster_summary_df <- do.call(rbind, lapply(names_list, function(name) {
     samples <- eval(parse(text = paste0("results_salter_", name, "$post_nclusters")))
     
     data.frame(
          model_prior = name,
          model = strsplit(name, "_")[[1]][1],
          prior = strsplit(name, "_")[[1]][2],
          ymin = min(samples),
          q1 = quantile(samples, 0.25),
          median = median(samples),
          q3 = quantile(samples, 0.75),
          ymax = max(samples)
     )
})) %>%
     mutate(index = as.numeric(factor(model_prior)))

# Plot with model-based color
pdf(file = paste0(base_path, "nclusters_posterior.pdf"), width = 10, height = 7)
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))

ggplot(cluster_summary_df, aes(x = index, color = model)) +
     # Min–max line
     geom_segment(aes(xend = index, y = ymin, yend = ymax), linewidth = 0.8) +
     
     # Quartile points
     geom_point(aes(y = q1), shape = 21, fill = "white", size = 2) +
     geom_point(aes(y = median), shape = 21, fill = "black", size = 3) +
     geom_point(aes(y = q3), shape = 21, fill = "white", size = 2) +
     
     # Axis and styling
     scale_x_continuous(breaks = cluster_summary_df$index,
                        labels = cluster_summary_df$model_prior) +
     scale_color_brewer(palette = "Set1") +
     theme_minimal(base_size = 14) +
     labs(x = "Model and Prior", y = "Posterior Number of Non-Empty Clusters",
          title = "Posterior Clustering Across Models and Priors",
          color = "Model") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Build all combination names
names_list <- expand.grid(model = models, prior = priors)
names_list <- apply(names_list, 1, function(x) paste(x, collapse = "_"))

# Extract posterior nclusters samples
cluster_samples <- lapply(names_list, function(name) {
     eval(parse(text = paste0("results_salter_", name, "$post_nclusters")))
})

# Create tidy dataframe with model and prior split
names(cluster_samples) <- names_list
cluster_df <- do.call(rbind, lapply(seq_along(cluster_samples), function(i) {
     model_prior <- names(cluster_samples)[i]
     parts <- unlist(strsplit(model_prior, "_"))
     data.frame(
          model = parts[1],
          prior = parts[2],
          model_prior = model_prior,
          n_clusters = cluster_samples[[i]]
     )
}))

summary_table <- do.call(rbind, lapply(seq_along(cluster_samples), function(i) {
     x <- cluster_samples[[i]]
     model_prior <- names(cluster_samples)[i]
     parts <- unlist(strsplit(model_prior, "_"))
     
     data.frame(
          model = parts[1],
          prior = parts[2],
          min = min(x),
          Q1 = quantile(x, 0.25),
          median = median(x),
          Q3 = quantile(x, 0.75),
          max = max(x)
     )
}))

rownames(summary_table) <- NULL

# View the table
print(summary_table)

# ==============================================================================
# Test statistics ==============================================================

# Define test statistic names
test_names <- c("Edge Density", "Transitivity", "Assortativity",
                "Mean Distance", "Mean Degree", "SD Degree")

# Models and priors
models <- c("class", "clist", "cleigen")
priors <- c("dm", "dp", "py", "gn")
names_list <- expand.grid(model = models, prior = priors)
names_list <- apply(names_list, 1, function(x) paste(x, collapse = "_"))

# Initialize storage
summary_df <- data.frame()

# Build summary with both 95% and 99% CIs
for (name in names_list) {
     test_post <- eval(parse(text = paste0("results_salter_", name, "$test_stats$test_stats")))
     test_obs  <- eval(parse(text = paste0("results_salter_", name, "$test_stats$test_obs")))
     
     for (i in seq_along(test_names)) {
          stat_samples <- test_post[, i]
          obs_value <- test_obs[i]
          
          quantiles_95 <- quantile(stat_samples, probs = c(0.025, 0.975))
          quantiles_99 <- quantile(stat_samples, probs = c(0.005, 0.995))
          median_value <- median(stat_samples)
          
          summary_df <- rbind(summary_df, data.frame(
               model_prior = name,
               model = strsplit(name, "_")[[1]][1],
               prior = strsplit(name, "_")[[1]][2],
               stat_name = test_names[i],
               lower95 = quantiles_95[1],
               upper95 = quantiles_95[2],
               lower99 = quantiles_99[1],
               upper99 = quantiles_99[2],
               median = median_value,
               observed = obs_value
          ))
     }
}

# Add label columns
summary_df <- summary_df %>%
     mutate(
          model_label = recode(model, "class" = "CM", "clist" = "CDM", "cleigen" = "CBM"),
          prior_label = toupper(prior),
          model_prior_label = prior_label  # <- Only show the prior in axis labels
     )

# Create and save one PDF per test statistic
for (stat in test_names) {
     
     df_stat <- summary_df %>% filter(stat_name == stat)
     df_stat <- df_stat %>% mutate(index = as.numeric(factor(model_prior)))
     
     p <- ggplot(df_stat, aes(x = index, y = median)) +
          # 99% CI
          geom_segment(aes(x = index, xend = index, y = lower99, yend = upper99, color = model_label),
                       linewidth = 1) +
          # 95% CI
          geom_segment(aes(x = index, xend = index, y = lower95, yend = upper95, color = model_label),
                       linewidth = 2) +
          # Median point
          geom_point(aes(color = model_label), size = 3) +
          # Observed value line
          geom_hline(aes(yintercept = observed), linetype = "dashed", color = "gray40") +
          # Labels and theme
          scale_x_continuous(breaks = df_stat$index, labels = df_stat$model_prior_label) +
          theme_minimal(base_size = 30) +
          labs(x = "Prior", y = stat, color = "Model") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
     
     ggsave(filename = paste0("plot_test_stats_", gsub(" ", "_", stat), "_salter.pdf"),
            plot = p, width = 10, height = 6)
}


# Build final summary table
final_table <- summary_df %>%
     select(model, prior, stat_name,
            observed,
            median,
            lower95, upper95,
            lower99, upper99) %>%
     arrange(stat_name, model, prior)  # Optional: organize nicely

rownames(final_table) <- NULL

# View Table
print(final_table)

# ==============================================================================
# Predictive metrics ===========================================================

# Define model and prior combinations
models <- c("class", "clist", "cleigen")
priors <- c("dm", "dp", "py", "gn")
names_list <- expand.grid(model = models, prior = priors)
names_list <- apply(names_list, 1, function(x) paste(x, collapse = "_"))

# Define metric names
metrics_names <- c("MSE", "LogLoss", "AUC")

# Build the performance summary
perf_summary_df <- data.frame()

for (name in names_list) {
     
     perf_stats <- eval(parse(text = paste0("results_salter_", name, "$pred_performance_stats")))
     
     perf_summary_df <- rbind(perf_summary_df, data.frame(
          model_prior = name,
          model = strsplit(name, "_")[[1]][1],
          prior = strsplit(name, "_")[[1]][2],
          
          metric = "MSE",
          mean = perf_stats$mean_mse,
          lower95 = perf_stats$ci_mse[1],
          upper95 = perf_stats$ci_mse[2]
     ))
     
     perf_summary_df <- rbind(perf_summary_df, data.frame(
          model_prior = name,
          model = strsplit(name, "_")[[1]][1],
          prior = strsplit(name, "_")[[1]][2],
          
          metric = "LogLoss",
          mean = perf_stats$mean_logloss,
          lower95 = perf_stats$ci_logloss[1],
          upper95 = perf_stats$ci_logloss[2]
     ))
     
     perf_summary_df <- rbind(perf_summary_df, data.frame(
          model_prior = name,
          model = strsplit(name, "_")[[1]][1],
          prior = strsplit(name, "_")[[1]][2],
          
          metric = "AUC",
          mean = perf_stats$mean_auc,
          lower95 = perf_stats$ci_auc[1],
          upper95 = perf_stats$ci_auc[2]
     ))
}

# Loop over metrics and generate a separate PDF for each
for (m in metrics_names) {
     
     df_metric <- perf_summary_df %>% filter(metric == m)
     df_metric <- df_metric %>% mutate(index = as.numeric(factor(model_prior)))
     
     p <- ggplot(df_metric, aes(x = index, y = mean)) +
          # 95% credible interval as thick vertical segment
          geom_segment(aes(x = index, xend = index,
                           y = lower95, yend = upper95, color = model),
                       linewidth = 1) +
          # Posterior mean as a point
          geom_point(aes(color = model), size = 2) +
          
          # Labeling and formatting
          scale_x_continuous(breaks = df_metric$index, labels = df_metric$model_prior) +
          theme_minimal(base_size = 14) +
          labs(x = "Model and Prior", y = m, title = paste("Predictive Performance:", m), color = "Model") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
     
     # Save as PDF
     ggsave(filename = paste0("predictive_performance_", m, ".pdf"),
            plot = p, width = 10, height = 6)
}

# Create the final predictive performance summary table
predictive_performance_table <- perf_summary_df %>%
     select(model, prior, metric, mean, lower95, upper95) %>%
     arrange(metric, model, prior)  # Optional: organized by metric first

rownames(predictive_performance_table) <- NULL

# View table
print(predictive_performance_table)

# ==============================================================================
# Results summary table ======================================================== 

results_df <- data.frame()

# List of models and priors
models <- c("class", "clist", "cleigen")
priors <- c("dm", "dp", "py", "gn")

# Loop through each model–prior combination
for (model in models) {
     for (prior in priors) {
          obj_name <- paste0("results_salter_", model, "_", prior)
          tmp <- data.frame(
               Model    = toupper(model),
               Prior    = toupper(prior),
               AUC      = round(get(obj_name)$pred_performance_stats$mean_auc, 3),
               MSE      = round(get(obj_name)$pred_performance_stats$mean_mse, 3),
               LogLoss  = round(get(obj_name)$pred_performance_stats$mean_logloss, 3)
          )
          results_df <- bind_rows(results_df, tmp)
     }
}

rownames(results_df) <- NULL

# View the final data frame
print(results_df)
