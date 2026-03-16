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
# ClusterSets_v2.1 — k-NN Jaccard Graph + Leiden Clustering for Large Collections
# ==========================================================
# Will Tower — 2026-03-16
# ==========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(Matrix)
  library(data.table)
  library(igraph)
})

# ----------------------------
# 1) File paths
# ----------------------------
in_file <- "~/path/to/ALL_Set_Models.csv"
out_dir <- "~/path/to/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_csv        <- file.path(out_dir, "Full_Collections_Clustered_kNN_Jaccard.csv")
out_unique_csv <- file.path(out_dir, "Unique_Collections_Clustered_kNN_Jaccard.csv")
out_edges      <- file.path(out_dir, "OverlapGraph_kNN_Jaccard_Edges.csv")
out_metrics    <- file.path(out_dir, "OverlapGraph_kNN_Jaccard_RunMetrics.txt")

# ----------------------------
# 2) Parameters
# ----------------------------
k_neighbors  <- 200
min_overlap  <- 1
resolution   <- 0.01
use_leiden   <- TRUE

# Optional safety caps
max_edges_pre_knn  <- Inf
max_edges_post_knn <- Inf

# ----------------------------
# 3) Helpers
# ----------------------------
stop_if_missing_col <- function(df, col) {
  if (!col %in% names(df)) stop("Missing required column: ", col)
}

clean_token <- function(x) {
  x <- gsub("\\.\\.\\.[0-9]+$", "", x)
  x <- gsub("\\.[0-9]+$", "", x)
  x
}

normalize_protein_string <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  toks <- strsplit(x, ",", fixed = TRUE)[[1]]
  toks <- toks[nchar(toks) > 0]
  toks <- clean_token(toks)
  toks <- unique(toks)
  paste(toks, collapse = ",")
}

write_metrics <- function(lines, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, con = path)
}

# ----------------------------
# 4) Load data
# ----------------------------
df <- readr::read_csv(in_file, show_col_types = FALSE)
message("✅ Loaded ", nrow(df), " total rows")
stop_if_missing_col(df, "Proteins")
df <- df %>% dplyr::mutate(Proteins = as.character(Proteins))

n_rows_original <- nrow(df)

# ----------------------------
# 5) Normalize Proteins and collapse to unique rows
# ----------------------------
message("🧬 Normalizing Proteins strings and collapsing duplicates...")

df$Proteins_Normalized <- vapply(df$Proteins, normalize_protein_string, character(1))

if (any(df$Proteins_Normalized == "" | is.na(df$Proteins_Normalized))) {
  bad <- which(df$Proteins_Normalized == "" | is.na(df$Proteins_Normalized))
  stop("Found ", length(bad), " rows with empty Proteins after normalization. Example rows: ",
       paste(head(bad, 20), collapse = ", "),
       if (length(bad) > 20) " ..." else "")
}

# Keep lookup so clusters can be mapped back later
dup_map <- data.frame(
  OriginalRow = seq_len(nrow(df)),
  Proteins_Normalized = df$Proteins_Normalized,
  stringsAsFactors = FALSE
)

# Unique rows only for graph construction
df_unique <- df %>%
  dplyr::distinct(Proteins_Normalized, .keep_all = TRUE)

n_unique <- nrow(df_unique)
n_dup_removed <- n_rows_original - n_unique

message("✅ Unique Protein rows: ", n_unique)
message("♻️ Duplicate rows collapsed: ", n_dup_removed)

# Use normalized string as Proteins for downstream clustering
df_unique$Proteins <- df_unique$Proteins_Normalized
n_sets <- nrow(df_unique)

# ----------------------------
# 6) Parse + normalize protein tokens
# ----------------------------
message("🧬 Parsing unique Proteins tokens...")

prot_str <- df_unique$Proteins
proteins_list <- strsplit(prot_str, ",", fixed = TRUE)

proteins_list <- lapply(proteins_list, function(v) {
  v <- v[nchar(v) > 0]
  v <- unique(v)
  v
})

if (any(lengths(proteins_list) == 0)) {
  bad <- which(lengths(proteins_list) == 0)
  stop("Found ", length(bad), " empty Protein lists after parsing unique rows. Example rows: ",
       paste(head(bad, 20), collapse = ", "),
       if (length(bad) > 20) " ..." else "")
}

all_proteins <- unique(unlist(proteins_list, use.names = FALSE))
message("🧬 Unique proteins: ", length(all_proteins))

# ----------------------------
# 7) Build sparse incidence matrix
# ----------------------------
message("🧮 Building sparse incidence matrix...")
lens <- lengths(proteins_list)
i_idx <- rep.int(seq_along(proteins_list), times = lens)
j_idx <- match(unlist(proteins_list, use.names = FALSE), all_proteins)

mat <- Matrix::sparseMatrix(
  i = i_idx,
  j = j_idx,
  x = 1L,
  dims = c(n_sets, length(all_proteins)),
  giveCsparse = TRUE
)

set_sizes <- Matrix::rowSums(mat)
message("📏 Unique set sizes: min=", min(set_sizes),
        " median=", stats::median(set_sizes),
        " max=", max(set_sizes))

# ----------------------------
# 8) Overlap matrix via tcrossprod
# ----------------------------
message("⚙️ Computing sparse overlap matrix (tcrossprod) on unique rows...")
overlap <- Matrix::tcrossprod(mat)
Matrix::diag(overlap) <- 0

overlap@x[overlap@x < min_overlap] <- 0
overlap <- Matrix::drop0(overlap)

nnz_ov <- length(overlap@x)
message("🔩 Nonzero overlaps (>= ", min_overlap, "): ", format(nnz_ov, big.mark=","))

if (nnz_ov == 0) {
  stop("No overlaps found at min_overlap = ", min_overlap,
       ". This usually indicates token mismatch or delimiter issues in Proteins.")
}

# ----------------------------
# 9) Extract edges
# ----------------------------
message("🔗 Extracting overlap edges...")
edges <- data.table::as.data.table(Matrix::summary(overlap))
data.table::setnames(edges, c("i", "j", "intersect"))

# Optional pre-cap by intersect strength
if (is.finite(max_edges_pre_knn) && nrow(edges) > max_edges_pre_knn) {
  message("⚠️ Pre-kNN edge cap applied: keeping top ", max_edges_pre_knn, " by intersect")
  data.table::setorder(edges, -intersect)
  edges <- edges[1:max_edges_pre_knn]
}

n_edges_pre_knn <- nrow(edges)
message("🔗 Pre-kNN edges: ", format(n_edges_pre_knn, big.mark=","))

# ----------------------------
# 10) Compute Jaccard and build symmetric k-NN graph
# ----------------------------
message("🧠 Computing Jaccard and building symmetric k-NN graph (k = ", k_neighbors, ")...")

edges[, jaccard := intersect / (set_sizes[i] + set_sizes[j] - intersect)]

# For each i, keep top-k by Jaccard
data.table::setorder(edges, i, -jaccard, -intersect, j)
edges_knn_i <- edges[, head(.SD, k_neighbors), by = i]

# For each j, keep top-k by Jaccard
data.table::setorder(edges, j, -jaccard, -intersect, i)
edges_knn_j <- edges[, head(.SD, k_neighbors), by = j]

edges_knn <- unique(rbind(
  edges_knn_i[, .(i, j, intersect, jaccard)],
  edges_knn_j[, .(i, j, intersect, jaccard)]
))

rm(edges, overlap)
gc()

# Optional post-cap by Jaccard strength
if (is.finite(max_edges_post_knn) && nrow(edges_knn) > max_edges_post_knn) {
  message("⚠️ Post-kNN edge cap applied: keeping top ", max_edges_post_knn, " by Jaccard")
  data.table::setorder(edges_knn, -jaccard, -intersect)
  edges_knn <- edges_knn[1:max_edges_post_knn]
}

n_edges_post_knn <- nrow(edges_knn)
message("🔗 Post-kNN edges: ", format(n_edges_post_knn, big.mark=","))

readr::write_csv(as.data.frame(edges_knn), out_edges)
message("💾 k-NN edges saved to: ", out_edges)

# ----------------------------
# 11) Build graph and cluster unique rows
# ----------------------------
message("🧩 Building igraph object...")
edge_df <- data.frame(
  from   = as.character(edges_knn$i),
  to     = as.character(edges_knn$j),
  weight = as.numeric(edges_knn$jaccard)
)

vert_df <- data.frame(name = as.character(seq_len(n_sets)))

g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = vert_df)

if (igraph::any_multiple(g) || igraph::any_loop(g)) {
  g <- igraph::simplify(
    g,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "max")
  )
}

deg <- igraph::degree(g)
comp <- igraph::components(g)

message("📈 Graph: ", igraph::gorder(g), " vertices; ", igraph::gsize(g), " edges")
message("🔌 Isolated vertices (degree 0): ", sum(deg == 0))
message("🧩 Components: ", comp$no, " | largest component size: ", max(comp$csize))
message("📊 Degree: min=", min(deg), " median=", stats::median(deg), " max=", max(deg))

do_leiden <- FALSE
if (use_leiden) do_leiden <- "cluster_leiden" %in% getNamespaceExports("igraph")

if (do_leiden) {
  message("🧩 Running Leiden (resolution = ", resolution, ")...")
  cl <- igraph::cluster_leiden(
    g,
    weights = igraph::E(g)$weight,
    resolution_parameter = resolution
  )
} else {
  message("🧩 Leiden not available; running Louvain (resolution ignored)...")
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
}

membership_vec <- igraph::membership(cl)
df_unique$Cluster <- as.integer(membership_vec[as.character(seq_len(n_sets))])

n_clusters <- length(unique(df_unique$Cluster))
message("✅ Identified ", n_clusters, " clusters on unique Protein rows")

# ----------------------------
# 12) Map clusters back to all original rows
# ----------------------------
message("🔁 Mapping unique-row clusters back to all original rows...")

cluster_lookup <- df_unique %>%
  dplyr::select(Proteins_Normalized, Cluster)

df_out <- df %>%
  dplyr::left_join(cluster_lookup, by = "Proteins_Normalized")

if (any(is.na(df_out$Cluster))) {
  stop("Cluster mapping failed: some original rows did not receive a cluster.")
}

# ----------------------------
# 13) Save outputs
# ----------------------------
readr::write_csv(df_unique, out_unique_csv)
message("💾 Unique clustered collections saved to: ", out_unique_csv)

readr::write_csv(df_out, out_csv)
message("💾 Full clustered collections saved to: ", out_csv)

# ----------------------------
# 14) Write run metrics
# ----------------------------
metrics_lines <- c(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Input file: ", in_file),
  paste0("Original rows: ", n_rows_original),
  paste0("Unique Proteins rows clustered: ", n_unique),
  paste0("Duplicate rows collapsed before clustering: ", n_dup_removed),
  paste0("Unique proteins: ", length(all_proteins)),
  paste0("Unique set size min/median/max: ", min(set_sizes), "/", stats::median(set_sizes), "/", max(set_sizes)),
  paste0("Similarity metric: Jaccard = intersect / (|A| + |B| - intersect)"),
  paste0("min_overlap: ", min_overlap),
  paste0("k_neighbors: ", k_neighbors),
  paste0("resolution: ", resolution),
  paste0("Edges pre-kNN: ", n_edges_pre_knn),
  paste0("Edges post-kNN: ", n_edges_post_knn),
  paste0("Graph edges: ", igraph::gsize(g)),
  paste0("Isolated vertices: ", sum(deg == 0)),
  paste0("Components: ", comp$no),
  paste0("Largest component size: ", max(comp$csize)),
  paste0("Degree min/median/max: ", min(deg), "/", stats::median(deg), "/", max(deg)),
  paste0("Clusters on unique rows: ", n_clusters),
  "",
  "Workflow:",
  " 1) Normalize Proteins strings",
  " 2) Collapse identical Proteins rows",
  " 3) Cluster only unique rows",
  " 4) Map cluster labels back to all original rows",
  "",
  "Tuning guidance:",
  " - If clusters >> 300: increase k_neighbors and/or decrease resolution.",
  " - If clusters << 100: decrease k_neighbors and/or increase resolution.",
  " - If many components persist: increase k_neighbors."
)
write_metrics(metrics_lines, out_metrics)
message("📝 Run metrics saved to: ", out_metrics)