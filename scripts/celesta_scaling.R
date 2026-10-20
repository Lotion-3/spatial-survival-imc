# Part IV, stage 1: choose the IMC input scaling for CELESTA on 20 seeded images per cohort.
#
#   Rscript scripts/celesta_scaling.R
#
# Grid: arcsinh cofactor {1, 5} x per-channel normalisation {none, zscore, p99} (see celesta_official.R;
# statistics are cohort-level). CELESTA (official package, default thresholds) is run on the same 20
# images per cohort under each setting. Only typing diagnostics are used -- no published labels, no
# subtype, genotype or outcome.
#
# Decision rule (fixed before running, AMBIGUITIES_GENOTYPE.md S1.2):
#   eligible = filtered fraction < 10% in both cohorts AND every round-1 class is assigned (>= 1 cell)
#              in >= 80% of the 20 images in both cohorts;
#   choose the eligible setting with the highest mean assigned fraction (cells not Unknown / Filtered)
#   over the two cohorts; ties -> fewer transforms (none < p99 < zscore), then cofactor 5 before 1.
#   If no setting is eligible, the class-presence bar is lowered to the best value any setting reaches
#   and the same ordering is applied (reported).
# Writes results/genotype/scaling/{grid.csv, per_image.csv, choice.json, ep_*.png}.

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2) })

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(root, "scripts", "celesta_official.R"))
out <- file.path(root, "results", "genotype", "scaling")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
SEED <- 20261019L
N_IMAGES <- 20L
ROUND1 <- c("Tumour", "T cell", "B cell", "Macrophage", "Endothelial", "Stroma")
GRID <- CJ(cofactor = c(1, 5), norm = c("none", "zscore", "p99"))
NORM_ORDER <- c(none = 1, p99 = 2, zscore = 3)

cohort_run <- function(cohort, cof, norm, images, opt_base) {
  sig <- read_signature(file.path(root, "results", "genotype", paste0("signature_", cohort, ".csv")))
  thr <- threshold_vectors(file.path(root, "results", "genotype", paste0("thresholds_", cohort, ".csv")), sig)
  markers <- colnames(sig)[3:ncol(sig)]
  cells <- CELLS[[cohort]][, c("cell_id", "image", "X", "Y", markers), with = FALSE]
  scaled <- scale_cohort(cells, markers, cof, norm)
  cells[, (markers) := as.data.table(scaled)]
  groups <- split(cells[image %in% images], by = "image", keep.by = TRUE)
  res <- parLapplyLB(CL, groups, type_one_image, sig = sig, thr = thr, opt = opt_base)
  labels <- rbindlist(lapply(res, `[[`, "labels"), fill = TRUE)
  ep <- rbindlist(lapply(res, `[[`, "ep"), fill = TRUE)
  list(labels = labels, ep = ep, status = sapply(res, `[[`, "status"), sig = sig)
}

main <- function() {
  opt_base <- parse_args(c("--cells", "x", "--signature", "x", "--thresholds", "x", "--out", "x"))
  CELLS <<- list(metabric = fread(file.path(root, "data", "processed", "genotype", "metabric_cells.csv.gz")),
                 basel = fread(file.path(root, "data", "processed", "genotype", "basel_cells.csv.gz")))
  set.seed(SEED)
  picks <- lapply(CELLS, function(d) sort(sample(sort(unique(d$image)), N_IMAGES)))
  fwrite(rbindlist(lapply(names(picks), function(c) data.table(cohort = c, image = picks[[c]], seed = SEED))),
         file.path(out, "images.csv"))
  CL <<- makeCluster(min(detectCores() - 2, 18))
  on.exit(stopCluster(CL))
  clusterExport(CL, "type_one_image")
  per_image <- list()
  grid <- list()
  eps <- list()
  for (g in seq_len(nrow(GRID))) {
    cof <- GRID$cofactor[g]
    norm <- GRID$norm[g]
    for (cohort in names(CELLS)) {
      r <- cohort_run(cohort, cof, norm, picks[[cohort]], opt_base)
      lab <- r$labels
      lab[, round1 := ifelse(celesta_type %in% c("CD8 T cell", "CD4 T cell", "Treg"), "T cell", celesta_type)]
      pres <- lab[, .(cls = ROUND1, present = ROUND1 %in% round1), by = image]
      pres <- pres[, .(images_with_class = mean(present)), by = cls]
      amb <- r$ep[, lapply(.SD, function(v) mean(v > 0.3 & v < 0.7)), .SDcols = !c("cell_id", "image")]
      grid[[length(grid) + 1]] <- data.table(
        cofactor = cof, norm = norm, cohort = cohort, n_cells = nrow(lab),
        assigned = mean(!lab$celesta_type %in% c("Unknown", "Filtered", "Skipped")),
        filtered = mean(lab$celesta_type == "Filtered"), unknown = mean(lab$celesta_type == "Unknown"),
        images_failed = sum(r$status != "ok"),
        min_class_presence = min(pres$images_with_class),
        mean_ambiguous_ep = mean(unlist(amb)))
      per_image[[length(per_image) + 1]] <- cbind(
        data.table(cofactor = cof, norm = norm, cohort = cohort),
        dcast(lab[, .N, by = .(image, celesta_type)], image ~ celesta_type, value.var = "N", fill = 0))
      for (cl_name in ROUND1) grid[[length(grid)]][[paste0("presence_", gsub(" ", "_", cl_name))]] <-
        pres$images_with_class[pres$cls == cl_name]
      eps[[paste(cof, norm, cohort)]] <- r
    }
  }
  grid <- rbindlist(grid, fill = TRUE)
  fwrite(grid, file.path(out, "grid.csv"))
  fwrite(rbindlist(per_image, fill = TRUE), file.path(out, "per_image.csv"))

  # ---- decision rule
  s <- grid[, .(assigned = mean(assigned), max_filtered = max(filtered), presence = min(min_class_presence)),
            by = .(cofactor, norm)]
  bar <- 0.8
  eligible <- s[max_filtered < 0.10 & presence >= bar]
  relaxed <- FALSE
  if (nrow(eligible) == 0) {
    bar <- max(s[max_filtered < 0.10, presence], -Inf)
    eligible <- s[max_filtered < 0.10 & presence >= bar]
    relaxed <- TRUE
  }
  if (nrow(eligible) == 0) stop("no setting has < 10% filtered cells in both cohorts")
  eligible[, norm_rank := NORM_ORDER[norm]]
  setorder(eligible, -assigned, norm_rank, -cofactor)
  choice <- eligible[1]
  write_json(list(rule = "filtered < 10% both cohorts; every round-1 class in >= 80% of images both cohorts; max mean assigned",
                  presence_bar_used = bar, relaxed = relaxed, seed = SEED, n_images_per_cohort = N_IMAGES,
                  summary = s, chosen = list(cofactor = choice$cofactor, norm = choice$norm,
                                             mean_assigned = choice$assigned, max_filtered = choice$max_filtered,
                                             min_class_presence = choice$presence)),
             file.path(out, "choice.json"), auto_unbox = TRUE, pretty = TRUE, digits = 4)
  print(s)
  cat(sprintf("chosen: cofactor %g, norm %s (presence bar %.2f%s)\n", choice$cofactor, choice$norm, bar,
              if (relaxed) ", relaxed" else ""))

  # ---- marker-probability plots (package's PlotExpProb) for the chosen and the package-like setting
  suppressPackageStartupMessages(library(CELESTA))
  for (key in unique(c(paste(choice$cofactor, choice$norm), "5 none"))) {
    for (cohort in names(CELLS)) {
      r <- eps[[paste(key, cohort)]]
      img <- picks[[cohort]][1]
      d <- CELLS[[cohort]][image == img]
      e <- r$ep[image == img]
      pdir <- file.path(out, paste0("ep_", cohort, "_", gsub(" ", "_", key)))
      dir.create(pdir, showWarnings = FALSE)
      PlotExpProb(coords = cbind(X = d$X, Y = d$Y), marker_exp_prob = as.matrix(e[, !c("cell_id", "image")]),
                  prior_marker_info = r$sig, save_plot = TRUE, output_dir = pdir)
    }
  }
}

if (sys.nframe() == 0) main()
