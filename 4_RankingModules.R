#!/usr/bin/env Rscript
# ==========================================================
# NIRRA — Network-Informed Restricted-set Ridge Analysis
#
# Copyright 2026 William Tower
#
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
# ==========================================================
# ==========================================================
# RankModules_v6.1 — Ranking modules by breadth–depth composite of extreme-tail predictive performance (p99)
# ==========================================================
# Will Tower — 2026-03-16
# ==========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# ==========================================================
# USER PARAMETERS
# ==========================================================

in_file <- "~/path/to/Full_Collections_Clustered_kNN_Jaccard.csv"
out_dir <- "~/path/to/output"

percentile_cutoff <- 0.99

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# LOAD DATA
# ==========================================================

df_all <- read_csv(in_file, show_col_types = FALSE)

cat("Rows loaded:", nrow(df_all), "\n")

# ==========================================================
# DETECT R² COLUMN
# ==========================================================

possible_r2 <- names(df_all)[grepl("ComboR2", names(df_all), ignore.case = TRUE)]

if (length(possible_r2) == 0) {
  stop("❌ ComboR2 column not found")
}

r2_col <- possible_r2[1]

cat("Using R² column:", r2_col, "\n")

# ==========================================================
# CHECK REQUIRED COLUMNS
# ==========================================================

if (!"Cluster" %in% names(df_all)) stop("❌ Cluster column missing")
if (!"SetIndex" %in% names(df_all)) stop("❌ SetIndex column missing")

# ==========================================================
# COLLAPSE TO BEST MODEL PER SETINDEX
# ==========================================================

cat("Selecting best model per SetIndex\n")

df_rank <- df_all %>%
  group_by(SetIndex) %>%
  slice_max(
    order_by = .data[[r2_col]],
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

cat("Sets after collapse:", nrow(df_rank), "\n")

# ==========================================================
# NORMALIZE CLUSTER IDS (DETERMINISTIC)
# ==========================================================

df_rank <- df_rank %>%
  mutate(Cluster = as.numeric(Cluster))

cluster_levels <- sort(unique(df_rank$Cluster[!is.na(df_rank$Cluster)]))

cluster_map <- tibble(
  Cluster = cluster_levels,
  Cluster_det = seq_along(cluster_levels)
)

df_rank <- df_rank %>%
  left_join(cluster_map, by = "Cluster") %>%
  mutate(Cluster = Cluster_det) %>%
  select(-Cluster_det)

cat("Clusters normalized:", length(cluster_levels), "\n")

# ==========================================================
# GLOBAL R² CUTOFF
# ==========================================================

global_cutoff <- quantile(
  df_rank[[r2_col]],
  percentile_cutoff,
  na.rm = TRUE
)

cat("Global cutoff:", round(global_cutoff,4), "\n")

# ==========================================================
# CLUSTER PERFORMANCE
# ==========================================================

cluster_perf <- df_rank %>%
  mutate(IsGlobalTop = .data[[r2_col]] >= global_cutoff) %>%
  group_by(Cluster) %>%
  summarise(
    
    n_sets = n(),
    
    CountTop1_Global = sum(IsGlobalTop, na.rm = TRUE),
    
    MedianTop1_Global =
      ifelse(
        CountTop1_Global > 0,
        median(.data[[r2_col]][IsGlobalTop], na.rm = TRUE),
        NA_real_
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    Score_geom_global =
      sqrt(CountTop1_Global * MedianTop1_Global),
    
    Rank_geom_global =
      if_else(
        is.na(Score_geom_global),
        NA_integer_,
        rank(-Score_geom_global, ties.method = "first")
      )
    
  ) %>%
  arrange(Rank_geom_global)

# ==========================================================
# CREATE MODULE COLUMN
# ==========================================================

ranked_clusters <- cluster_perf %>%
  filter(!is.na(Rank_geom_global)) %>%
  arrange(Rank_geom_global) %>%
  mutate(Module = as.character(Rank_geom_global))

unranked_clusters <- cluster_perf %>%
  filter(is.na(Rank_geom_global)) %>%
  arrange(Cluster)

n_ranked <- nrow(ranked_clusters)

if (nrow(unranked_clusters) > 0) {
  
  unranked_clusters <- unranked_clusters %>%
    mutate(
      Module = paste0(seq_len(n()) + n_ranked, " (unranked)")
    )
  
}

cluster_perf <- bind_rows(
  ranked_clusters,
  unranked_clusters
)

# ==========================================================
# WARNING FOR UNRANKED CLUSTERS
# ==========================================================

if (nrow(unranked_clusters) > 0) {
  
  warning(
    paste0(
      nrow(unranked_clusters),
      " clusters had no sets in the global top percentile."
    )
  )
  
} else {
  
  cat("All clusters ranked\n")
  
}

# ==========================================================
# MERGE MODULES BACK TO FULL DATASET
# ==========================================================

df_out <- df_all %>%
  left_join(
    cluster_perf %>%
      select(
        Cluster,
        Module,
        Rank_geom_global,
        CountTop1_Global,
        MedianTop1_Global,
        Score_geom_global
      ),
    by = "Cluster"
  )

# ==========================================================
# SAVE OUTPUTS
# ==========================================================

summary_file <- file.path(
  out_dir,
  "Cluster_Performance_GlobalTop1pct.csv"
)

full_file <- file.path(
  out_dir,
  "Combo_RidgeResults_AllModels_WithModules.csv"
)

write_csv(cluster_perf, summary_file)
write_csv(df_out, full_file)

cat("Saved summary:", summary_file, "\n")
cat("Saved full dataset:", full_file, "\n")

# ==========================================================
# PLOT ALL RANKED MODULES
# ==========================================================

plot_df <- cluster_perf %>%
  filter(!grepl("unranked", Module)) %>%
  mutate(
    Module = factor(Module, levels = rev(sort(as.numeric(Module))))
  )

plot_file <- file.path(
  out_dir,
  "Cluster_Performance_AllRanked_Modules.pdf"
)

pdf(
  plot_file,
  width = 10,
  height = 6 + nrow(plot_df) * 0.15
)

ggplot(plot_df, aes(x = Module, y = Score_geom_global)) +
  geom_col(fill = "darkblue") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Module Ranking by Global Top-1% Performance",
    subtitle = paste0(
      "Score = √(count × median R²), cutoff ≥ ",
      round(global_cutoff,3)
    ),
    x = "Module",
    y = "Composite Score"
  )

dev.off()

cat("Saved plot:", plot_file, "\n")

# ==========================================================
# PRINT TOP MODULES
# ==========================================================

cat("\nTop modules:\n")

print(
  head(
    cluster_perf %>%
      select(
        Module,
        Cluster,
        CountTop1_Global,
        MedianTop1_Global,
        Score_geom_global
      ),
    10
  )
)