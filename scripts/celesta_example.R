# Part IV, stage 0: run the official CELESTA package's bundled example (README "Usage" block,
# demo data = reg009 of Schurch et al. 2020 CODEX) and compare the final labels with the
# reference assignment shipped in the package's tests/testthat/.
#
#   Rscript scripts/celesta_example.R
#
# Writes results/genotype/example/{example_summary.json, *_final_cell_type_assignment.csv, plots}.

suppressPackageStartupMessages({
  library(CELESTA)
  library(Rmixmod)
  library(spdep)
  library(ggplot2)
  library(reshape2)
  library(zeallot)
  library(jsonlite)
})

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
out <- file.path(root, "results", "genotype", "example")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ref_dir <- file.path(root, "third_party", "CELESTA", "tests", "testthat")

# The package's lazy-loaded demo objects (data/*.rda)
data("prior_marker_info", "imaging_data", "high_marker_threshold_anchor", "low_marker_threshold_anchor",
     "high_marker_threshold_iteration", "low_marker_threshold_iteration", package = "CELESTA")

t0 <- Sys.time()
obj <- CreateCelestaObject(project_title = file.path(out, "example"), prior_marker_info, imaging_data)
obj <- FilterCells(obj, high_marker_threshold = 0.9, low_marker_threshold = 0.4)
obj <- AssignCells(obj, max_iteration = 10, cell_change_threshold = 0.01,
                   high_expression_threshold_anchor = high_marker_threshold_anchor,
                   low_expression_threshold_anchor = low_marker_threshold_anchor,
                   high_expression_threshold_index = high_marker_threshold_iteration,
                   low_expression_threshold_index = low_marker_threshold_iteration)
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

final <- obj@final_cell_type_assignment[, "Final cell type"]
PlotCellsAnyCombination(cell_type_assignment_to_plot = obj@final_cell_type_assignment[, obj@total_rounds + 1],
                        coords = obj@coords, prior_info = prior_marker_info,
                        cell_number_to_use = c(0, 1, 2, 3), cell_type_colors = c("yellow", "red", "blue"),
                        output_dir = out)
# The demo marker names contain ":", which is not allowed in Windows file names (ggsave would
# write NTFS alternate data streams), so the plot copy uses names with ":" -> "_".
safe <- function(x) gsub(":", "_", x, fixed = TRUE)
ep <- obj@marker_exp_prob
colnames(ep) <- safe(colnames(ep))
pmi <- prior_marker_info
colnames(pmi) <- safe(colnames(pmi))
PlotExpProb(coords = obj@coords, marker_exp_prob = ep, prior_marker_info = pmi, save_plot = TRUE, output_dir = out)

ref <- read.csv(file.path(ref_dir, "project_title_final_cell_type_assignment.csv"), check.names = FALSE)
ref_final <- ref[["Final cell type"]]
same <- (is.na(final) & is.na(ref_final)) | (!is.na(final) & !is.na(ref_final) & final == ref_final)
summary <- list(
  n_cells = length(final),
  runtime_seconds = round(runtime, 1),
  assigned_fraction = mean(!is.na(final) & final != "Unknown"),
  filtered_fraction = mean(is.na(final)),
  label_counts = as.list(table(final, useNA = "ifany")),
  reference_file = "third_party/CELESTA/tests/testthat/project_title_final_cell_type_assignment.csv",
  agreement_with_reference = mean(same),
  completed = TRUE
)
write_json(summary, file.path(out, "example_summary.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("example: %d cells, %.1f s, assigned %.3f, agreement with shipped reference %.4f\n",
            summary$n_cells, runtime, summary$assigned_fraction, summary$agreement_with_reference))
