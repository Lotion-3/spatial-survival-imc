# Part IV, stage 1 report: CELESTA cell typing on METABRIC and Basel.
#
#   Rscript scripts/celesta_report.R
#
# Writes results/genotype/celesta/:
#   class_share.csv          share of cells per assigned class (incl. Unknown / Filtered), per cohort
#   composition_<cohort>.csv per-image cell counts per class (input to stage 2 covariates)
#   acceptance.csv           per class: share of images assigned >= 1 cell, among images where the
#                            class's defining marker is present (>= 5 cells above the cohort-level
#                            crossing point from CELESTA's own BuildSigmoidFunction GMM)
#   marker_models.csv        that crossing point per marker and cohort
#   concordance_<cohort>.csv cross-tabulation vs. the published labels (concordance only)
#   concordance_summary.csv  agreement and Cohen's kappa on cells both methods assign
#   marker_means_<cohort>.csv mean scaled expression per marker by assigned class
# and figures/genotype_marker_heatmap.png, figures/genotype_celesta_example_<cohort>.png (package
# PlotCellsAnyCombination).

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(CELESTA); library(Rmixmod); library(zeallot); library(jsonlite)
})

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
source(file.path(root, "scripts", "celesta_official.R"))  # scale_cohort(), read_signature()
gen <- file.path(root, "results", "genotype")
out <- file.path(gen, "celesta")
fig <- file.path(root, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
SEED <- 20261019L
SCALING <- fromJSON(file.path(gen, "scaling", "choice.json"))$chosen
DEFINING <- c(Tumour = "panCK", `T cell` = "CD3", `B cell` = "CD20", Macrophage = "CD68",
              Endothelial = "CD31-vWF", Stroma = "SMA", `CD8 T cell` = "CD8", `CD4 T cell` = "CD4", Treg = "FOXP3")
T_SUBSETS <- c("CD8 T cell", "CD4 T cell", "Treg")
ORDER <- c(names(DEFINING), "Unknown", "Filtered", "Skipped")
# Basel metacluster coarse labels / METABRIC published coarse labels use these names
PUB <- c(Tumor = "Tumour", T = "T cell", B = "B cell", Macrophage = "Macrophage", Endothelial = "Endothelial",
         Stroma = "Stroma")
EXTRA <- list(metabric = c("CK8-18", "CK5", "CD45RA", "CD45RO", "PDGFRB", "FSP1", "Podoplanin", "HER2", "Ki-67"),
              basel = c("CK8-18", "CK5", "CD45", "Vimentin", "Fibronectin", "HER2", "Ki-67"))

kappa <- function(a, b) {
  lv <- union(a, b)
  t <- table(factor(a, lv), factor(b, lv))
  po <- sum(diag(t)) / sum(t)
  pe <- sum(rowSums(t) * colSums(t)) / sum(t)^2
  (po - pe) / (1 - pe)
}

report_cohort <- function(cohort) {
  sig <- read_signature(file.path(gen, paste0("signature_", cohort, ".csv")))
  markers <- colnames(sig)[3:ncol(sig)]
  cells <- fread(file.path(root, "data", "processed", "genotype", paste0(cohort, "_cells.csv.gz")))
  lab <- fread(file.path(gen, paste0("cell_labels_", cohort, ".csv.gz")))
  d <- merge(cells, lab[, .(cell_id, celesta_type)], by = "cell_id")
  shown <- c(markers, EXTRA[[cohort]])
  scaled <- scale_cohort(d, shown, SCALING$cofactor, SCALING$norm)
  d[, (shown) := as.data.table(scaled)]

  # share of cells per class
  share <- d[, .N, by = celesta_type][, share := N / sum(N)][, cohort := cohort]
  # per-image composition
  comp <- dcast(d[, .N, by = .(image, PID, celesta_type)], image + PID ~ celesta_type, value.var = "N", fill = 0)
  fwrite(comp, file.path(out, paste0("composition_", cohort, ".csv")))

  # cohort-level marker positivity: CELESTA's own GMM crossing point, on a seeded subsample
  set.seed(SEED)
  sub <- as.matrix(d[sample(.N, min(.N, 100000L)), ..markers])
  invisible(capture.output(sp <- BuildSigmoidFunction(sub)))
  xroot <- setNames(sp[1, ], markers)
  models <- data.table(cohort = cohort, marker = markers, crossing_point = xroot,
                       frac_cells_above = sapply(markers, function(m) mean(d[[m]] > xroot[m])))
  # acceptance: class assigned in >= 80% of images where its defining marker is present
  acc <- rbindlist(lapply(intersect(names(DEFINING), sig[[1]]), function(cls) {
    m <- DEFINING[[cls]]
    per_img <- d[, .(present = sum(get(m) > xroot[m]) >= 5, assigned = sum(celesta_type == cls)), by = image]
    data.table(cohort = cohort, cell_type = cls, defining_marker = m, images = nrow(per_img),
               images_marker_present = per_img[, sum(present)],
               images_assigned = per_img[present == TRUE, sum(assigned >= 1)],
               share_assigned = per_img[present == TRUE, mean(assigned >= 1)])
  }))
  acc[, passes_80 := share_assigned >= 0.8]

  # marker means by class
  mm <- d[!celesta_type %in% c("Filtered", "Skipped"), lapply(.SD, mean), by = celesta_type, .SDcols = shown]
  fwrite(mm, file.path(out, paste0("marker_means_", cohort, ".csv")))

  # concordance with published labels (round-1 level; T subsets -> T cell)
  d[, celesta_r1 := fifelse(celesta_type %in% T_SUBSETS, "T cell", celesta_type)]
  d[, published_r1 := PUB[published]]
  ct <- dcast(d[, .N, by = .(published_r1, celesta_r1)], published_r1 ~ celesta_r1, value.var = "N", fill = 0)
  fwrite(ct, file.path(out, paste0("concordance_", cohort, ".csv")))
  if (cohort == "metabric") {  # T-cell subclasses vs. the published phenotypes (concordance only)
    pub_t <- c("CD8 T cell" = "CD8^{+} T cells", "CD4 T cell" = "CD4^{+} T cells", "Treg" = "T_{Reg} & T_{Ex}")
    sub <- rbindlist(lapply(names(pub_t), function(cls) data.table(
      celesta_class = cls, published_phenotype = pub_t[[cls]], celesta_cells = d[celesta_type == cls, .N],
      published_cells = d[published_phenotype == pub_t[[cls]], .N],
      both = d[celesta_type == cls & published_phenotype == pub_t[[cls]], .N],
      celesta_cells_published_epithelial = d[celesta_type == cls & published == "Tumor", .N])))
    sub[, `:=`(share_of_celesta_in_published = both / celesta_cells, share_of_published_in_celesta = both / published_cells)]
    fwrite(sub, file.path(out, "concordance_tsubsets_metabric.csv"))
  }
  both <- d[celesta_r1 %in% PUB]
  conc <- data.table(cohort = cohort, cells = nrow(d), cells_celesta_assigned = nrow(both),
                     agreement_assigned = both[, mean(celesta_r1 == published_r1)],
                     kappa_assigned = kappa(both$celesta_r1, both$published_r1))

  # example core: the image with the median assigned fraction, plotted with the package function
  af <- d[, .(assigned = mean(!celesta_type %in% c("Unknown", "Filtered", "Skipped")), n = .N), by = image]
  af <- af[n >= 500]
  ex <- af[order(abs(assigned - median(assigned)))][1, image]
  e <- d[image == ex]
  num <- match(e$celesta_type, sig[[1]])
  num[e$celesta_type == "Unknown"] <- 0
  num[is.na(num)] <- 0  # Filtered shown with Unknown
  types_shown <- seq_len(nrow(sig))
  cols <- c("#c0392b", "#e67e22", "#2e86c1", "#8e44ad", "#16a085", "#95a5a6", "#f1c40f", "#a3e635", "#ff4fd8")
  g <- PlotCellsAnyCombination(cell_type_assignment_to_plot = num, coords = cbind(e$X, e$Y), prior_info = sig,
                               cell_number_to_use = c(types_shown, 0),
                               cell_type_colors = c(cols[types_shown], "grey40"), test_size = 0.8,
                               save_plot = FALSE)
  ggsave(file.path(fig, paste0("genotype_celesta_example_", cohort, ".png")),
         g + ggtitle(paste0(cohort, " ", ex, ": official CELESTA labels")), width = 7, height = 6, dpi = 150)

  list(share = share, acc = acc, models = models, conc = conc, mm = mm[, cohort := cohort], example = ex)
}

main <- function() {
  res <- lapply(c(metabric = "metabric", basel = "basel"), report_cohort)
  share <- rbindlist(lapply(res, `[[`, "share"))
  share[, celesta_type := factor(celesta_type, ORDER)]
  setorder(share, cohort, celesta_type)
  fwrite(share, file.path(out, "class_share.csv"))
  acc <- rbindlist(lapply(res, `[[`, "acc"))
  fwrite(acc, file.path(out, "acceptance.csv"))
  fwrite(rbindlist(lapply(res, `[[`, "models")), file.path(out, "marker_models.csv"))
  fwrite(rbindlist(lapply(res, `[[`, "conc")), file.path(out, "concordance_summary.csv"))
  write_json(lapply(res, `[[`, "example"), file.path(out, "example_cores.json"), auto_unbox = TRUE)

  # heatmap: mean scaled expression by class, each marker z-scored across classes, per cohort
  mm <- rbindlist(lapply(res, function(r) melt(r$mm, id.vars = c("celesta_type", "cohort"), variable.name = "marker")),
                  fill = TRUE)
  mm[, z := (value - mean(value)) / sd(value), by = .(cohort, marker)]
  mm[, celesta_type := factor(celesta_type, rev(ORDER))]
  g <- ggplot(mm, aes(marker, celesta_type, fill = z)) + geom_tile(colour = "white") +
    facet_wrap(~cohort, scales = "free_x") +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", name = "z (across\nclasses)") +
    labs(x = NULL, y = "CELESTA class",
         title = "Mean scaled marker expression by official-CELESTA class",
         subtitle = sprintf("arcsinh(x / %g), divided by cohort p99 x 10; each marker standardised across classes",
                            SCALING$cofactor)) +
    theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(fig, "genotype_marker_heatmap.png"), g, width = 11, height = 4.8, dpi = 160)
  print(acc)
  print(share[, .(cohort, celesta_type, N, share = round(share, 4))])
  print(rbindlist(lapply(res, `[[`, "conc")))
}

if (sys.nframe() == 0) main()
