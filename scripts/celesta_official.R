# Part IV: cell typing with the official CELESTA R package (plevritis-lab/CELESTA), one image at a time.
#
#   Rscript scripts/celesta_official.R --cells cells.csv.gz --signature sig.csv --thresholds thr.csv \
#       --out labels.csv.gz [--cofactor 5] [--norm none|zscore|p99] [--workers 18] [--images a,b,...] \
#       [--filter 1] [--ep-out ep.csv.gz] [--status-out status.csv]
#
# Inputs (CSV, explicit columns):
#   --cells       cell_id:int, image:str, X:num, Y:num, then one numeric column per marker (raw IMC mean
#                 ion counts). Extra columns are ignored.
#   --signature   CELESTA prior_marker_info layout: column 1 = cell type, column 2 = Lineage_level
#                 ("round_parent_number"), then one column per marker with 1 / 0 / NA.
#   --thresholds  cell_type, high_anchor, low_anchor, high_index, low_index; one row per signature row,
#                 in the same order (row i = cell type number i).
# Output (CSV): cell_id:int, image:str, celesta_type:str ("Unknown" = unassigned, "Filtered" = removed by
#   FilterCells, "Skipped" = image not typed, see status), cell_type_number:int (0 Unknown, NA otherwise),
#   round_1..round_R:str, anchor_round1:bool.
#
# Scaling: CreateCelestaObject() takes a `cofactor` argument but passes a hard-coded cofactor of 10 to
# GetMarkerExpMatrix() (CELESTA_functions.R, commit e5e92ff), so the transform is done here and CELESTA
# is called with transform_type = 0. Per-channel normalisation statistics are computed over all cells in
# the --cells file (cohort level), never per image, so they cannot differ between images of one cohort.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

parse_args <- function(args) {
  out <- list(cofactor = "5", norm = "none", workers = "1", filter = "1", images = "", ep_out = "",
              status_out = "", number_of_neighbors = "5", bandwidth = "100", max_iteration = "10",
              cell_change_threshold = "0.01", filter_high = "0.9", filter_low = "0.4")
  i <- 1
  while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[i]))
    out[[key]] <- args[i + 1]
    i <- i + 2
  }
  for (k in c("cells", "signature", "thresholds", "out")) if (is.null(out[[k]])) stop("missing --", k)
  out
}

read_signature <- function(path) {
  sig <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
  sig[[2]] <- as.character(sig[[2]])
  sig
}

threshold_vectors <- function(path, sig) {
  thr <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(identical(thr$cell_type, sig[[1]]))
  num <- as.integer(sapply(strsplit(sig[[2]], "_"), `[`, 3))
  vec <- function(col) {
    v <- rep(NA_real_, max(num, 50))
    v[num] <- thr[[col]]
    v
  }
  list(high_anchor = vec("high_anchor"), low_anchor = vec("low_anchor"),
       high_index = vec("high_index"), low_index = vec("low_index"))
}

# Cohort-level scaling of the raw counts; returns a function applied per image.
scale_cohort <- function(cells, markers, cofactor, norm) {
  x <- asinh(as.matrix(cells[, ..markers]) / cofactor)
  if (norm == "zscore") {
    mu <- colMeans(x)
    s <- apply(x, 2, sd)
    x <- sweep(sweep(x, 2, mu), 2, s, "/")
  } else if (norm == "p99") {
    # Divide each channel by its cohort 99th percentile and multiply by 10, so every channel spans
    # roughly the arcsinh(x/10) range of the package's CODEX demo data (about 0-10).
    p99 <- apply(x, 2, quantile, 0.99)
    x <- sweep(x, 2, p99 / 10, "/")
  } else if (norm != "none") {
    stop("unknown --norm ", norm)
  }
  x
}

type_one_image <- function(d, sig, thr, opt) {
  suppressPackageStartupMessages({
    library(data.table); library(CELESTA); library(Rmixmod); library(spdep); library(zeallot)
  })
  markers <- colnames(sig)[3:ncol(sig)]
  img <- data.frame(X = d$X, Y = d$Y, d[, markers, with = FALSE], check.names = FALSE)
  rounds <- max(as.integer(sapply(strsplit(sig[[2]], "_"), `[`, 1)))
  empty <- function(status) list(
    labels = data.table(cell_id = d$cell_id, image = d$image, celesta_type = "Skipped",
                        cell_type_number = NA_integer_, anchor_round1 = FALSE),
    ep = NULL, status = status)
  log <- character()
  res <- tryCatch({
    log <- capture.output({
      obj <- CreateCelestaObject(project_title = "unused", prior_marker_info = sig, imaging_data_file = img,
                                 transform_type = 0,
                                 number_of_neighbors = as.integer(opt$number_of_neighbors),
                                 bandwidth = as.numeric(opt$bandwidth))
    })
    if (is.null(obj)) return(empty("too_few_cells"))
    if (opt$filter == "1") {
      log <- c(log, capture.output(obj <- FilterCells(obj, high_marker_threshold = as.numeric(opt$filter_high),
                                                      low_marker_threshold = as.numeric(opt$filter_low))))
    }
    log <- c(log, capture.output(obj <- AssignCells(
      obj, max_iteration = as.integer(opt$max_iteration),
      cell_change_threshold = as.numeric(opt$cell_change_threshold),
      high_expression_threshold_anchor = thr$high_anchor, low_expression_threshold_anchor = thr$low_anchor,
      high_expression_threshold_index = thr$high_index, low_expression_threshold_index = thr$low_index,
      bandwidth = as.numeric(opt$bandwidth), save_result = FALSE)))
    fa <- obj@final_cell_type_assignment
    num <- suppressWarnings(as.integer(fa[, "Cell type number"]))
    lab <- fa[, "Final cell type"]
    filtered <- is.na(obj@current_cell_type_assignment[, 1])
    lab[filtered] <- "Filtered"
    num[filtered] <- NA_integer_
    out <- data.table(cell_id = d$cell_id, image = d$image, celesta_type = lab, cell_type_number = num,
                      anchor_round1 = !filtered & obj@anchor_cell_type_assignment[, 1] != 0)
    for (r in seq_len(rounds)) {
      rl <- fa[, r]
      rl[filtered] <- "Filtered"
      out[[paste0("round_", r)]] <- rl
    }
    ep <- data.table(cell_id = d$cell_id, image = d$image, obj@marker_exp_prob)
    list(labels = out, ep = ep, status = "ok")
  }, error = function(e) empty(paste("error:", conditionMessage(e))))
  res
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  sig <- read_signature(opt$signature)
  thr <- threshold_vectors(opt$thresholds, sig)
  markers <- colnames(sig)[3:ncol(sig)]
  cells <- fread(opt$cells, colClasses = list(character = "image"))
  need <- c("cell_id", "image", "X", "Y", markers)
  miss <- setdiff(need, names(cells))
  if (length(miss)) stop("cells file lacks columns: ", paste(miss, collapse = ", "))
  cells <- cells[, ..need]
  scaled <- scale_cohort(cells, markers, as.numeric(opt$cofactor), opt$norm)
  cells[, (markers) := as.data.table(scaled)]
  if (nzchar(opt$images)) cells <- cells[image %in% strsplit(opt$images, ",")[[1]]]
  groups <- split(cells, by = "image", keep.by = TRUE)
  workers <- as.integer(opt$workers)
  t0 <- Sys.time()
  if (workers > 1) {
    cl <- makeCluster(workers)
    on.exit(stopCluster(cl))
    clusterExport(cl, c("type_one_image"), envir = environment())
    res <- parLapplyLB(cl, groups, type_one_image, sig = sig, thr = thr, opt = opt)
  } else {
    res <- lapply(groups, type_one_image, sig = sig, thr = thr, opt = opt)
  }
  labels <- rbindlist(lapply(res, `[[`, "labels"), fill = TRUE)
  setcolorder(labels, c("cell_id", "image", "celesta_type", "cell_type_number"))
  fwrite(labels, opt$out)
  if (nzchar(opt$ep_out)) fwrite(rbindlist(lapply(res, `[[`, "ep")), opt$ep_out)
  status <- data.table(image = names(groups), n_cells = sapply(groups, nrow), status = sapply(res, `[[`, "status"))
  if (nzchar(opt$status_out)) fwrite(status, opt$status_out)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("celesta_official: %d images (%d ok), %d cells, %.1f s, %d workers\n", nrow(status),
              sum(status$status == "ok"), nrow(labels), secs, workers))
}

if (sys.nframe() == 0) main()
