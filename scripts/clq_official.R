# Part IV, stage 2: co-location quotients per image with the Plevritis lab's own CLQ code
# (plevritis-lab/Spatial_Permutation_and_Normalization, Spatial_Perm_and_Normalization.R, pinned in
# results/genotype/environment.json), run unmodified.
#
#   Rscript scripts/clq_official.R --labels labels.csv.gz --cells cells.csv.gz --types types.csv \
#       --out clq.csv.gz [--k 10,20,40] [--perms 500] [--seed 20261019] [--workers 18] \
#       [--label-col cell_type_number] [--keep-unknown 0] [--seeds-out seeds.csv]
#
# The lab script is a driver file: its top-level loops glob *cell_type_assignment.csv in the working
# directory and hard-code 18 cell types and k = 20. We therefore evaluate only its function definitions
# (get_CLQ, find_cell_type_neighbors, KNN_neighbors, significance_matrix_gen), byte-for-byte from the
# original file, and drive them here for any number of cell types and k. The permutation step repeats
# the lab driver's own loop (sample() of the whole label vector, coordinates fixed).
#
# CLQ as implemented by get_CLQ(cell_A = a, cell_B = b):
#   CLQ = (C_ab / N_a) / (N_b / (N - 1))      [(N_b - 1) for a == b]
# where C_ab = number of type-b cells among the k nearest neighbours of all type-a cells, and N counts
# typed cells (label not 0 / NA). C_ab / N_a is a COUNT out of k, so this CLQ equals k x the
# Leslie-Kronenfeld fraction-based CLQ (expected value k, not 1, under random labelling). We report
# both: clq (lab scale) and clq_k = clq / k (classic scale, null mean 1). get_CLQ returns 0 when either
# type has <= 5 cells; those pairs are set to NA here (as the lab's CLQ_normalization_by_sample does).
#
# Inputs (CSV): labels = cell_id, image, <label-col> (integer type number, 0 = Unknown, NA = filtered;
#   or class names, mapped to type numbers via --types, other names -> 0);
#   cells = cell_id, image, X, Y; types = cell_type (row i = type number i).
# Output (CSV, one row per image x k x ordered pair): image, k, cell_A, cell_B, n_A, n_B, n_typed, n_all,
#   clq, clq_k, perm_mean, perm_sd, percentile (lab rule), significant (lab rule: percentile <= 0.05 or
#   >= 0.95), p_two_sided (new: (1 + #{|perm - perm_mean| >= |clq - perm_mean|}) / (1 + perms)).

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

LAB_FUNCTIONS <- c("get_CLQ", "find_cell_type_neighbors", "KNN_neighbors", "significance_matrix_gen")
LAB <- new.env(parent = globalenv())  # the lab's functions, evaluated from the original file

lab_script <- function(root) {
  file.path(root, "third_party", "Spatial_Permutation_and_Normalization", "Spatial_Perm_and_Normalization.R")
}

# Evaluate the first definition of each named function from the lab file into `env`, unmodified.
load_lab_functions <- function(path, env = LAB) {
  exprs <- parse(path, keep.source = FALSE)
  done <- character()
  for (e in exprs) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% setdiff(LAB_FUNCTIONS, done)) {
      eval(e, envir = env)
      done <- c(done, as.character(e[[2]]))
    }
  }
  stopifnot(setequal(done, LAB_FUNCTIONS))
  invisible(done)
}

parse_args <- function(args) {
  out <- list(k = "10,20,40", perms = "500", seed = "20261019", workers = "1", label_col = "cell_type_number",
              keep_unknown = "0", seeds_out = "")
  i <- 1
  while (i <= length(args)) {
    out[[gsub("-", "_", sub("^--", "", args[i]))]] <- args[i + 1]
    i <- i + 2
  }
  for (k in c("labels", "cells", "types", "out")) if (is.null(out[[k]])) stop("missing --", k)
  out
}

# CLQ matrix (types x types) with the lab's get_CLQ. get_CLQ fails when cell_A has no cells
# (find_cell_type_neighbors returns NULL); the lab's permutation driver guards this by writing 0 for
# every pair of an absent type (`if (is.na(count_a)) CLQ_array = append(CLQ_array, rep(0, 18))`),
# which is repeated here.
clq_matrix <- function(lab, coords, nb_list, n_types, counts) {
  m <- matrix(0, n_types, n_types)
  for (a in seq_len(n_types)) {
    if (counts[a] == 0) next
    for (b in seq_len(n_types)) {
      m[a, b] <- LAB$get_CLQ(lab, coords, nb_list, a, b, number_of_neighbor = NA, as.character(b))
    }
  }
  m
}

# One image: CLQ at every k plus `perms` label permutations (the same permutations for every k).
clq_one_image <- function(d, types, ks, perms, seed, keep_unknown, script) {
  suppressPackageStartupMessages({
    library(data.table); library(spdep); library(dplyr)
  })
  if (!exists("get_CLQ", envir = LAB, inherits = FALSE)) load_lab_functions(script, LAB)
  n_types <- length(types)
  if (!keep_unknown) d <- d[!is.na(lab) & lab != 0]
  lab <- d$lab
  coords <- cbind(X = d$X, Y = d$Y)
  n_typed <- sum(!is.na(lab) & lab != 0)
  counts <- tabulate(lab[!is.na(lab) & lab != 0], nbins = n_types)
  out <- list()
  if (nrow(d) <= max(ks)) return(NULL)  # kNN needs more cells than k
  set.seed(seed)
  perm_labels <- lapply(seq_len(perms), function(r) sample(lab, length(lab), replace = FALSE))
  for (k in ks) {
    nb_list <- LAB$KNN_neighbors(coords, number_of_neighbors = k)
    obs <- clq_matrix(lab, coords, nb_list, n_types, counts)
    perm <- vapply(perm_labels, function(pl) as.vector(t(clq_matrix(pl, coords, nb_list, n_types, counts))),
                   numeric(n_types^2))  # rows ordered cellA-major, as in the lab's CLQ_matrix_R
    # Lab significance rule, run unmodified in a scratch directory (it reads cell_types_celesta.csv and
    # writes <sample>_CLQ_data_full.csv into the working directory).
    wd <- tempfile("clq_")
    dir.create(wd)
    old <- setwd(wd)
    write.csv(data.frame(cell_types = types), "cell_types_celesta.csv", row.names = FALSE)
    perm_df <- data.frame(cellA = rep(seq_len(n_types), each = n_types), cellB = rep(seq_len(n_types), n_types))
    perm_df <- cbind(perm_df, setNames(as.data.frame(perm), paste0("P", seq_len(perms))))
    obs_named <- obs
    dimnames(obs_named) <- list(types, types)
    df_c <- data.frame(cell_type = types, cell_type_assignment = seq_len(n_types),
                       Freq = ifelse(counts > 0, counts, NA))
    sig <- LAB$significance_matrix_gen(iternum = perms, filename = "img_cell_type_assignment.csv",
                                   list_of_matrices = perm_df, CLQ_matrix_original = obs_named, df_c = df_c)
    setwd(old)
    unlink(wd, recursive = TRUE)
    obs_v <- as.vector(t(obs))
    pm <- rowMeans(perm)
    p2 <- (1 + rowSums(abs(perm - pm) >= abs(obs_v - pm) - 1e-12)) / (1 + perms)
    a <- rep(seq_len(n_types), each = n_types)
    b <- rep(seq_len(n_types), n_types)
    rare <- counts[a] <= 5 | counts[b] <= 5
    out[[as.character(k)]] <- data.table(
      image = d$image[1], k = k, cell_A = types[a], cell_B = types[b], n_A = counts[a], n_B = counts[b],
      n_typed = n_typed, n_all = nrow(d),
      clq = ifelse(rare, NA_real_, obs_v), clq_k = ifelse(rare, NA_real_, obs_v / k),
      perm_mean = ifelse(rare, NA_real_, pm), perm_sd = ifelse(rare, NA_real_, apply(perm, 1, sd)),
      percentile = ifelse(rare, NA_real_, as.numeric(sig$percentile)),
      significant = ifelse(rare, NA, as.logical(sig$significant)),
      p_two_sided = ifelse(rare, NA_real_, p2))
  }
  rbindlist(out)
}

run_index <- function(i) clq_one_image(groups[[i]], types, ks, perms, seeds$seed[i], keep_unknown, script)

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
  script <- lab_script(root)
  types <- read.csv(opt$types, stringsAsFactors = FALSE)$cell_type
  labels <- fread(opt$labels, colClasses = list(character = "image"))
  cells <- fread(opt$cells, select = c("cell_id", "image", "X", "Y"), colClasses = list(character = "image"))
  d <- merge(cells, labels[, c("cell_id", "image", opt$label_col), with = FALSE], by = c("cell_id", "image"))
  setnames(d, opt$label_col, "lab")
  if (is.character(d$lab)) {  # class names (e.g. published labels) -> type numbers; other names -> 0
    d[, lab := match(lab, types, nomatch = 0L)]
  }
  d[, lab := as.integer(lab)]
  ks <- as.integer(strsplit(opt$k, ",")[[1]])
  perms <- as.integer(opt$perms)
  base_seed <- as.integer(opt$seed)
  images <- sort(unique(d$image))
  seeds <- data.table(image = images, seed = base_seed + seq_along(images))
  # Shared state lives in the global environment (and in each worker's), so the task function
  # run_index() is a top-level function and its closure carries no data when serialised.
  shared <- list(groups = split(d, by = "image", keep.by = TRUE)[images], types = types, ks = ks, perms = perms,
                 seeds = seeds, keep_unknown = opt$keep_unknown == "1", script = script)
  list2env(shared, .GlobalEnv)
  t0 <- Sys.time()
  workers <- as.integer(opt$workers)
  if (workers > 1) {
    cl <- makeCluster(workers)
    on.exit(stopCluster(cl))
    clusterExport(cl, c("clq_one_image", "clq_matrix", "load_lab_functions", "LAB_FUNCTIONS", "LAB", names(shared)),
                  envir = .GlobalEnv)
    res <- parLapplyLB(cl, seq_along(images), run_index)
  } else {
    res <- lapply(seq_along(images), run_index)
  }
  out <- rbindlist(res)
  fwrite(out, opt$out)
  if (nzchar(opt$seeds_out)) fwrite(seeds[, .(image, seed, perms = perms, k = opt$k)], opt$seeds_out)
  cat(sprintf("clq_official: %d images, k = %s, %d permutations, %.1f min, %d workers\n", length(images), opt$k,
              perms, as.numeric(difftime(Sys.time(), t0, units = "mins")), workers))
}

if (sys.nframe() == 0) main()
