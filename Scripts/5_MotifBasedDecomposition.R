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
# MotifDecomposition+HierarchicalEmpiricalBayesEvaluation_v2.5 — IntraModule motif generation with global evaluation (HEB)
# ==========================================================
# Will Tower — 2026-03-16
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(lme4)
  library(tibble)
})

# =========================
# USER INPUTS
# =========================
in_csv      <- "~/path/to/Combo_RidgeResults_AllModels_WithModules.csv"
min_frac    <- 0.05
min_carrier <- 5

out_csv <- "~/path/to/output/motif_global_shrinkage_ModuleSpecific.csv"

# =========================
# LOAD
# =========================
df <- read.csv(in_csv, stringsAsFactors = FALSE)

required_cols <- c("SetIndex", "Module", "ComboR2", "Coefficients")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0)
  stop("Missing required columns: ", paste(missing_cols, collapse=", "))

df <- df %>%
  dplyr::mutate(
    SetID   = as.character(SetIndex),
    Module = as.factor(Module),
    ComboR2 = as.numeric(ComboR2)
  )

Modules <- sort(unique(df$Module))

# =========================
# FUNCTIONS
# =========================
parse_coefficients <- function(coeff_string, min_frac = 0.05) {
  
  if (is.na(coeff_string) || !nzchar(coeff_string)) return(NULL)
  
  parts <- stringr::str_split(coeff_string, ",\\s*")[[1]]
  prot_beta <- stringr::str_split(parts, "=")
  
  proteins <- sapply(prot_beta, `[`, 1)
  betas    <- suppressWarnings(as.numeric(sapply(prot_beta, `[`, 2)))
  
  valid <- is.finite(betas)
  proteins <- proteins[valid]
  betas    <- betas[valid]
  
  if (length(betas) < 2) return(NULL)
  
  abs_w <- abs(betas)
  total_w <- sum(abs_w)
  if (!is.finite(total_w) || total_w == 0) return(NULL)
  
  keep <- (abs_w / total_w) >= min_frac
  if (sum(keep) < 2) return(NULL)
  
  tibble(
    Protein = proteins[keep],
    Beta    = betas[keep],
    Sign    = ifelse(betas[keep] > 0, "+", "-")
  )
}

generate_pairs <- function(prot_table) {
  
  if (is.null(prot_table) || nrow(prot_table) < 2) return(NULL)
  
  combs <- combn(nrow(prot_table), 2)
  
  motifs <- apply(combs, 2, function(idx) {
    p1 <- paste0(prot_table$Sign[idx[1]], prot_table$Protein[idx[1]])
    p2 <- paste0(prot_table$Sign[idx[2]], prot_table$Protein[idx[2]])
    paste(sort(c(p1, p2)), collapse="|")
  })
  
  unique(motifs)
}

# =========================
# GENERATE INTRAModule MOTIFS
# =========================
all_rows <- list()

for (cl in Modules) {
  
  df_Module <- df %>%
    dplyr::filter(Module == cl)
  
  motif_df <- df_Module %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      parsed = list(parse_coefficients(Coefficients, min_frac)),
      motifs = list(generate_pairs(parsed))
    ) %>%
    dplyr::ungroup()
  
  motif_long <- motif_df %>%
    dplyr::select(SetID, ComboR2, motifs) %>%
    tidyr::unnest(motifs) %>%
    dplyr::rename(Motif = motifs)
  
  if (nrow(motif_long) == 0) next
  
  motif_long <- motif_long %>%
    dplyr::mutate(
      Motif   = as.character(Motif),
      SetID   = as.character(SetID),
      ComboR2 = as.numeric(ComboR2),
      Module = cl
    ) %>%
    dplyr::filter(is.finite(ComboR2)) %>%
    dplyr::distinct(Motif, SetID, .keep_all=TRUE)
  
  # intraModule carrier filter
  motif_counts <- motif_long %>%
    dplyr::group_by(Motif) %>%
    dplyr::summarise(
      n_carriers = dplyr::n_distinct(SetID),
      .groups="drop"
    ) %>%
    dplyr::filter(n_carriers >= min_carrier)
  
  motif_long_f <- motif_long %>%
    dplyr::filter(Motif %in% motif_counts$Motif)
  
  if (nrow(motif_long_f) == 0) next
  
  all_rows[[as.character(cl)]] <- motif_long_f
}

# =========================
# BIND ALL ModuleS
# =========================
global_motif_long <- dplyr::bind_rows(all_rows)

# Module-specific motif identity
global_motif_long <- global_motif_long %>%
  dplyr::mutate(
    MotifKey = interaction(Module, Motif, drop = TRUE)
  )

# =========================
# COMPUTE INTRAModule CARRIER COUNTS
# =========================
carrier_counts <- global_motif_long %>%
  dplyr::group_by(MotifKey) %>%
  dplyr::summarise(
    n_carriers = dplyr::n_distinct(SetID),
    .groups = "drop"
  )

# =========================
# GLOBAL SHRINKAGE MODEL
# =========================
fit <- lmer(
  ComboR2 ~ 1 + (1 | MotifKey),
  data = global_motif_long,
  REML = TRUE
)

# =========================
# EXTRACT EFFECTS
# =========================
re <- ranef(fit, condVar = TRUE)$MotifKey
pv <- attr(re, "postVar")

motif_sd <- sqrt(sapply(seq_len(dim(pv)[3]), function(i) pv[1,1,i]))
grand_mean <- fixef(fit)[1]

motif_results <- tibble(
  MotifKey  = rownames(re),
  alpha_hat = re[,1],
  alpha_sd  = motif_sd
) %>%
  dplyr::mutate(
    MotifEstimatedPredictiveness = grand_mean + alpha_hat,
    Lower95 = MotifEstimatedPredictiveness - 1.96 * alpha_sd,
    Upper95 = MotifEstimatedPredictiveness + 1.96 * alpha_sd
  ) %>%
  # split back into Module and Motif
  dplyr::mutate(
    Module = sub("\\..*$", "", MotifKey),
    Motif   = sub("^[^\\.]+\\.", "", MotifKey)
  ) %>%
  dplyr::left_join(carrier_counts, by = "MotifKey") %>%
  dplyr::arrange(desc(MotifEstimatedPredictiveness))

# =========================
# SAVE
# =========================
write.csv(motif_results, out_csv, row.names = FALSE)

cat("Saved to:", out_csv, "\n")