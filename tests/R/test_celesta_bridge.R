# Output format of the R bridge scripts/celesta_official.R on a small synthetic input.

suppressPackageStartupMessages({ library(testthat); library(data.table) })

root <- normalizePath(file.path(Sys.getenv("SPATIALSURV_ROOT", ".")))
rscript <- file.path(R.home("bin"), "Rscript")

test_that("celesta_official.R writes one typed row per input cell with the documented columns", {
  set.seed(11)
  dir <- tempfile("bridge_")
  dir.create(dir)
  n <- 600
  type <- sample(c("Tumour", "T cell", "Other"), n, replace = TRUE, prob = c(0.6, 0.25, 0.15))
  # raw counts: bimodal per marker, as in IMC mean ion counts
  panck <- ifelse(type == "Tumour", rlnorm(n, log(8), 0.3), rlnorm(n, log(0.3), 0.5))
  cd3 <- ifelse(type == "T cell", rlnorm(n, log(5), 0.3), rlnorm(n, log(0.2), 0.5))
  cells <- data.table(cell_id = 1000L + seq_len(n), image = rep(c("imgA", "imgB"), each = n / 2),
                      X = runif(n, 0, 400), Y = runif(n, 0, 400), panCK = panck, CD3 = cd3, extra = 1)
  fwrite(cells, file.path(dir, "cells.csv"))
  fwrite(data.table(cell_type = c("Tumour", "T cell"), Lineage_level = c("1_0_1", "1_0_2"),
                    panCK = c(1, 0), CD3 = c(0, 1)), file.path(dir, "sig.csv"))
  fwrite(data.table(cell_type = c("Tumour", "T cell"), high_anchor = 0.7, low_anchor = 0.9,
                    high_index = 0.5, low_index = 1), file.path(dir, "thr.csv"))
  status <- system2(rscript, c(file.path(root, "scripts", "celesta_official.R"),
                               "--cells", file.path(dir, "cells.csv"), "--signature", file.path(dir, "sig.csv"),
                               "--thresholds", file.path(dir, "thr.csv"), "--out", file.path(dir, "labels.csv"),
                               "--cofactor", "1", "--norm", "p99", "--status-out", file.path(dir, "status.csv")),
                    stdout = TRUE, stderr = TRUE)
  expect_true(file.exists(file.path(dir, "labels.csv")), info = paste(status, collapse = "\n"))
  lab <- fread(file.path(dir, "labels.csv"))
  expect_equal(names(lab), c("cell_id", "image", "celesta_type", "cell_type_number", "anchor_round1", "round_1"))
  expect_type(lab$cell_id, "integer")
  expect_type(lab$cell_type_number, "integer")
  expect_type(lab$anchor_round1, "logical")
  expect_setequal(lab$cell_id, cells$cell_id)
  expect_true(all(lab$celesta_type %in% c("Tumour", "T cell", "Unknown", "Filtered")))
  expect_true(all(lab[celesta_type == "Unknown", cell_type_number] == 0L))
  expect_true(all(is.na(lab[celesta_type == "Filtered", cell_type_number])))
  st <- fread(file.path(dir, "status.csv"))
  expect_equal(st$status, c("ok", "ok"))
  # the synthetic types are easy: most tumour cells come back as Tumour
  m <- merge(lab, cells[, .(cell_id, truth = type)], by = "cell_id")
  expect_gt(m[truth == "Tumour", mean(celesta_type == "Tumour")], 0.8)
})
