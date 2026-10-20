# Part IV, stage 2: per-patient spatial features from the per-image CLQ tables, and the stage-2
# acceptance summary (permutation null centred on 1).
#
#   Rscript scripts/clq_features.R
#
# For every run results/genotype/clq/<run>.csv.gz (written by clq_official.R):
#   features/<run>_clq.csv.gz      PID, k, cell_A, cell_B, clq_k (typed-cell-weighted mean over the patient's
#                                  images where the pair is defined), n_images
#   features/<run>_composition.csv PID, n_images, n_typed, log_n_typed, count_<class>, clr_<class>
# and clq/acceptance.csv: per run and k, the distribution of the permutation mean CLQ / k over
# image x pair entries (should be centred on 1), plus the share of pairs flagged by the lab rule.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
gen <- file.path(root, "results", "genotype")
clq_dir <- file.path(gen, "clq")
feat_dir <- file.path(gen, "features")
dir.create(feat_dir, showWarnings = FALSE)
PSEUDO <- 0.5  # as spatialsurv.features.clr

clr <- function(counts) {
  l <- log(counts + PSEUDO)
  l - rowMeans(l)
}

image_pid <- function(cohort) {
  unique(fread(file.path(root, "data", "processed", "genotype", paste0(cohort, "_cells.csv.gz")),
               select = c("image", "PID"), colClasses = list(character = c("image", "PID"))))
}

features_for_run <- function(run) {
  cohort <- sub("_.*", "", run)
  d <- fread(file.path(clq_dir, paste0(run, ".csv.gz")))
  d <- merge(d, image_pid(cohort), by = "image")
  types <- unique(d$cell_A)
  # per-image class counts (identical for every k and partner; take one k and the diagonal)
  cnt <- unique(d[k == min(k) & cell_A == cell_B, .(image, PID, cell = cell_A, n = n_A, n_typed)])
  comp <- dcast(cnt, PID ~ cell, value.var = "n", fun.aggregate = sum)
  setcolorder(comp, c("PID", types))
  imgs <- cnt[, .(n_images = uniqueN(image), n_typed = sum(n) ), by = PID]
  comp <- merge(imgs, comp, by = "PID")
  cm <- as.matrix(comp[, ..types])
  cl <- clr(cm)
  comp[, log_n_typed := log(n_typed)]
  setnames(comp, types, paste0("count_", types))
  comp <- cbind(comp, setNames(as.data.table(cl), paste0("clr_", types)))
  fwrite(comp, file.path(feat_dir, paste0(run, "_composition.csv")))

  p <- d[!is.na(clq_k), .(clq_k = weighted.mean(clq_k, n_typed), n_images = .N), by = .(PID, k, cell_A, cell_B)]
  fwrite(p, file.path(feat_dir, paste0(run, "_clq.csv.gz")))

  acc <- d[!is.na(perm_mean), .(
    entries = .N, images = uniqueN(image),
    perm_mean_over_k_median = median(perm_mean / k), perm_mean_over_k_q025 = quantile(perm_mean / k, 0.025),
    perm_mean_over_k_q975 = quantile(perm_mean / k, 0.975), perm_mean_over_k_mean = mean(perm_mean / k),
    observed_clq_k_median = median(clq_k), share_flagged_lab_rule = mean(significant),
    share_p_two_sided_lt_005 = mean(p_two_sided < 0.05)), by = k]
  acc[, run := run]
  cat(sprintf("%s: %d patients, %d images\n", run, nrow(comp), uniqueN(d$image)))
  acc
}

main <- function() {
  runs <- sub("\\.csv\\.gz$", "", list.files(clq_dir, pattern = "^(metabric|basel)_.*\\.csv\\.gz$"))
  acc <- rbindlist(lapply(runs, features_for_run))
  setcolorder(acc, "run")
  acc[, centred_on_1 := abs(perm_mean_over_k_median - 1) <= 0.05]
  fwrite(acc, file.path(clq_dir, "acceptance.csv"))
  print(acc[, .(run, k, images, perm_med = round(perm_mean_over_k_median, 4), q025 = round(perm_mean_over_k_q025, 3),
                q975 = round(perm_mean_over_k_q975, 3), centred_on_1)])
}

if (sys.nframe() == 0) main()
