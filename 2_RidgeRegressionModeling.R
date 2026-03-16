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
# RunRidge_v5.2 — Ridge Regression Modeling (unconstrained and sign-constrained)
# ==========================================================
# Will Tower — 2026-03-16
# ==========================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

set.seed(123)

eps_beta <- 1e-6
tol_yhat <- 1e-4
tol_beta <- 1e-4

# ==========================================================
# PATHS
# ==========================================================
in_file  <- "~/path/to/inputmatrix.csv"
ppi_file <- "~/path/to/PPI_HighDensity_Collections_AllHubs.csv"
out_dir  <- "~/path/to/output"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_models_csv <- file.path(out_dir, "ALL_Set_Models.csv")
out_scores_csv <- file.path(out_dir, "SubjectScores_AllSets.csv")
out_check_csv  <- file.path(out_dir, "SetScore_R2_Equals_ComboR2.csv")
out_metrics    <- file.path(out_dir, "RunMetrics.txt")
out_equiv_csv  <- file.path(out_dir, "Equivalence_Audit.csv")

# ==========================================================
# LOAD DATA
# ==========================================================
df <- read_csv(in_file, show_col_types = FALSE)

if (!"Outcome" %in% names(df)) stop("❌ Column 'Outcome' not found.")
if (!"Subject" %in% names(df)) stop("❌ Column 'Subject' not found.")

df <- df[is.finite(df$Outcome), ]
y  <- df$Outcome
n_subj <- length(y)

prot_start <- which(names(df) == "FirstProtein")
prot_end   <- which(names(df) == "LastProtein")

X_all <- df[, prot_start:prot_end]

X_all <- as.data.frame(lapply(X_all, function(v) {
  v <- as.numeric(v)
  bad <- !is.finite(v)
  if (any(bad)) {
    med <- median(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    v[bad] <- med
  }
  v
}))

# ==========================================================
# LOAD PPI SETS
# ==========================================================
ppi <- read_csv(ppi_file, show_col_types = FALSE) |>
  filter(!is.na(PPI_p)) |>
  mutate(
    Proteins = str_replace_all(Proteins, "\\s*,\\s*", ","),
    Proteins = str_trim(Proteins)
  )

n_sets_total <- nrow(ppi)

# ==========================================================
# RIDGE FITTER
# ==========================================================
ridge_fit <- function(X, y, mode) {
  
  lower <- switch(mode,
                  free   = rep(-Inf, ncol(X)),
                  nonneg = rep(0,    ncol(X)),
                  nonpos = rep(-Inf, ncol(X)))
  
  upper <- switch(mode,
                  free   = rep( Inf, ncol(X)),
                  nonneg = rep( Inf, ncol(X)),
                  nonpos = rep(0,    ncol(X)))
  
  cv <- tryCatch(
    cv.glmnet(
      X, y,
      alpha = 0,
      standardize = FALSE,
      lower.limits = lower,
      upper.limits = upper
    ),
    error = function(e) NULL
  )
  
  if (is.null(cv)) return(NULL)
  
  lam  <- cv$lambda.min
  yhat <- as.numeric(predict(cv$glmnet.fit, X, s = lam))
  
  sst <- sum((y - mean(y))^2)
  sse <- sum((y - yhat)^2)
  R2  <- max(0, 1 - sse / sst)
  
  beta <- as.numeric(coef(cv, s = lam))[-1]
  names(beta) <- colnames(X)
  
  list(R2 = R2, yhat = yhat, beta = beta, lambda = lam, cv = cv)
}

# ==========================================================
# EQUIVALENCE AUDIT
# ==========================================================
equiv_audit <- function(X, free, con, which) {
  
  if (is.null(free) || is.null(con)) return(NULL)
  
  beta_free <- free$beta
  
  violates <- FALSE
  if (which == "nonneg" && any(beta_free < -eps_beta)) violates <- TRUE
  if (which == "nonpos" && any(beta_free >  eps_beta)) violates <- TRUE
  
  if (violates) {
    return(tibble(
      Which = which,
      FreeSatisfiesConstraint = FALSE
    ))
  }
  
  lam <- free$lambda
  
  y1 <- as.numeric(predict(free$cv$glmnet.fit, X, s = lam))
  y2 <- as.numeric(predict(con$cv$glmnet.fit,  X, s = lam))
  
  max_dy <- max(abs(y1 - y2))
  
  b1 <- as.numeric(coef(free$cv$glmnet.fit, s = lam))[-1]
  b2 <- as.numeric(coef(con$cv$glmnet.fit,  s = lam))[-1]
  
  max_db <- max(abs(b1 - b2))
  
  tibble(
    Which = which,
    FreeSatisfiesConstraint = TRUE,
    MaxAbsDeltaYhat = max_dy,
    MaxAbsDeltaBeta = max_db
  )
}

# ==========================================================
# RUN MODELS
# ==========================================================
results  <- list()
scores   <- list()
equivlog <- list()

skipped_small <- 0L
skipped_fit   <- 0L

start_time <- Sys.time()

for (i in seq_len(n_sets_total)) {
  
  genes <- intersect(
    strsplit(ppi$Proteins[i], ",")[[1]] |> str_trim(),
    colnames(X_all)
  )
  
  if (length(genes) < 2) {
    skipped_small <- skipped_small + 1
    next
  }
  
  X <- scale(as.matrix(X_all[, genes, drop = FALSE]))
  colnames(X) <- genes
  
  fit_r <- ridge_fit(X, y, "free")
  fit_p <- ridge_fit(X, y, "nonneg")
  fit_n <- ridge_fit(X, y, "nonpos")
  
  if (is.null(fit_r) && is.null(fit_p) && is.null(fit_n)) {
    skipped_fit <- skipped_fit + 1
    next
  }
  
  a_p <- equiv_audit(X, fit_r, fit_p, "nonneg")
  a_n <- equiv_audit(X, fit_r, fit_n, "nonpos")
  
  if (!is.null(a_p)) equivlog[[length(equivlog)+1]] <- a_p |> mutate(SetIndex=i)
  if (!is.null(a_n)) equivlog[[length(equivlog)+1]] <- a_n |> mutate(SetIndex=i)
  
  model_list <- list(
    Ridge  = fit_r,
    NonNeg = fit_p,
    NonPos = fit_n
  )
  
  for (m in names(model_list)) {
    
    fit <- model_list[[m]]
    if (is.null(fit)) next
    
    coef_str <- paste0(
      names(fit$beta),"=",sprintf("%.5f",fit$beta),
      collapse=", "
    )
    
    results[[length(results)+1]] <- tibble(
      SetIndex=i,
      Model=m,
      Proteins=ppi$Proteins[i],
      NumVars=length(genes),
      ComboR2=fit$R2,
      Lambda=fit$lambda,
      Coefficients=coef_str
    )
    
    scores[[length(scores)+1]] <- fit$yhat
  }
  
  if (i %% 25 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units="mins")
    message("Processed ",i," / ",n_sets_total," | elapsed ",round(elapsed,2)," min")
  }
}

# ==========================================================
# OUTPUT
# ==========================================================
final_df <- bind_rows(results)

score_mat <- do.call(cbind, scores)

colnames(score_mat) <- paste0(
  "Set_",final_df$SetIndex,"_",final_df$Model
)

rownames(score_mat) <- df$Subject

write_csv(final_df,out_models_csv)

scores_df <- as.data.frame(t(score_mat))
scores_df <- tibble::rownames_to_column(scores_df,"SetID")

write_csv(scores_df,out_scores_csv)

# ==========================================================
# VERIFY R²
# ==========================================================
sst <- sum((y - mean(y))^2)

check_df <- lapply(seq_len(ncol(score_mat)), function(j){
  
  sc <- score_mat[,j]
  sse <- sum((y-sc)^2)
  
  r2_score <- max(0,1-sse/sst)
  
  tibble(
    SetID=colnames(score_mat)[j],
    ComboR2=final_df$ComboR2[j],
    R2_score=r2_score,
    Diff=r2_score-final_df$ComboR2[j]
  )
  
}) |> bind_rows()

write_csv(check_df,out_check_csv)

# ==========================================================
# METRICS
# ==========================================================
cat(
  "N subjects:",n_subj,"\n",
  "N PPI sets:",n_sets_total,"\n",
  "N models:",nrow(final_df),"\n",
  "Skipped small:",skipped_small,"\n",
  "Skipped fit:",skipped_fit,"\n",
  file=out_metrics
)

equiv_df <- bind_rows(equivlog)
write_csv(equiv_df,out_equiv_csv)

message("✔ COMPLETE — All three models retained")
