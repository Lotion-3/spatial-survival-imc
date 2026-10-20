# CLQ tests for the Plevritis lab code as driven by scripts/clq_official.R.
# Run: Rscript tests/R/run_tests.R  (or via pytest: tests/test_r_scripts.py)

suppressPackageStartupMessages({ library(testthat); library(data.table); library(spdep); library(dplyr) })

root <- normalizePath(file.path(Sys.getenv("SPATIALSURV_ROOT", ".")))
source(file.path(root, "scripts", "clq_official.R"))
script <- lab_script(root)
load_lab_functions(script, LAB)

# Ten pairs of cells 100 um apart, partners 1 um apart: with k = 1 each cell's only neighbour is its
# partner. 6 A-B pairs, 2 A-A pairs, 2 B-B pairs -> N_A = N_B = 10, N = 20.
#   C_AB = 6 (A cells with a B neighbour):  CLQ(A->B) = (6/10) / (10/19)       = 1.14
#   C_AA = 4:                               CLQ(A->A) = (4/10) / ((10-1)/19)  = 0.8444...
pairs <- c(rep("AB", 6), rep("AA", 2), rep("BB", 2))
toy_coords <- cbind(X = rep(seq(0, by = 100, length.out = 10), each = 2) + rep(c(0, 1), 10), Y = 0)
toy_lab <- as.integer(unlist(lapply(pairs, function(p) c(AB = c(1, 2), AA = c(1, 1), BB = c(2, 2))[
  paste0(p, 1:2)])))

test_that("lab get_CLQ matches the hand-computed toy example", {
  expect_equal(toy_lab, c(rep(c(1L, 2L), 6), 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))
  nb <- LAB$KNN_neighbors(toy_coords, number_of_neighbors = 1)
  m <- clq_matrix(toy_lab, toy_coords, nb, 2, tabulate(toy_lab, 2))
  expect_equal(m[1, 2], (6 / 10) / (10 / 19))
  expect_equal(m[1, 1], (4 / 10) / (9 / 19))
  expect_equal(m[2, 1], (6 / 10) / (10 / 19))
  expect_equal(m[2, 2], (4 / 10) / (9 / 19))
})

test_that("lab CLQ is a neighbour count: k x the fraction-based CLQ, about k under random labels", {
  set.seed(1)
  n <- 3000
  xy <- cbind(X = runif(n, 0, 1000), Y = runif(n, 0, 1000))
  lab <- sample(1:4, n, replace = TRUE, prob = c(0.55, 0.15, 0.15, 0.15))
  for (k in c(10, 20)) {
    nb <- LAB$KNN_neighbors(xy, number_of_neighbors = k)
    m <- clq_matrix(lab, xy, nb, 4, tabulate(lab, 4))
    expect_lt(abs(mean(m / k) - 1), 0.03)  # classic scale centred on 1
    expect_lt(max(abs(m / k - 1)), 0.15)
  }
})

test_that("rare types (<= 5 cells) give 0 in the lab code and NA in our output; absent types do not crash", {
  set.seed(2)
  n <- 400
  d <- data.table(image = "toy", X = runif(n, 0, 300), Y = runif(n, 0, 300),
                  lab = c(rep(3L, 4), sample(1:2, n - 4, replace = TRUE)))
  res <- clq_one_image(d, types = c("A", "B", "C", "D"), ks = 10, perms = 20, seed = 1, keep_unknown = FALSE,
                       script = script)
  expect_true(all(is.na(res[cell_A == "C" | cell_B == "C", clq])))
  expect_true(all(is.na(res[cell_A == "D" | cell_B == "D", clq])))  # type D absent from the image
  expect_false(any(is.na(res[cell_A %in% c("A", "B") & cell_B %in% c("A", "B"), clq])))
})

test_that("within-image label permutations are centred on 1 (classic scale)", {
  set.seed(3)
  n <- 1500
  d <- data.table(image = "toy", X = runif(n, 0, 800), Y = runif(n, 0, 800),
                  lab = sample(c(0L, 1:3), n, replace = TRUE, prob = c(0.1, 0.5, 0.2, 0.2)))
  # spatial structure: type 2 clustered in one corner -> observed self-CLQ far above 1
  d[X < 250 & Y < 250 & lab != 0, lab := 2L]
  res <- clq_one_image(d, types = c("A", "B", "C"), ks = 20, perms = 100, seed = 7, keep_unknown = FALSE,
                       script = script)
  expect_lt(abs(mean(res$perm_mean / 20) - 1), 0.02)
  expect_gt(res[cell_A == "B" & cell_B == "B", clq_k], 1.3)
  expect_lt(res[cell_A == "B" & cell_B == "B", p_two_sided], 0.05)
  expect_true(res[cell_A == "B" & cell_B == "B", significant])
  expect_equal(sort(unique(res$k)), 20)
  # the same seed reproduces the permutations exactly
  res2 <- clq_one_image(d, types = c("A", "B", "C"), ks = 20, perms = 100, seed = 7, keep_unknown = FALSE,
                        script = script)
  expect_identical(res$perm_mean, res2$perm_mean)
})

test_that("k values outside perm_ks get the identical observed CLQ and NA permutation columns", {
  set.seed(5)
  n <- 600
  d <- data.table(image = "toy", X = runif(n, 0, 400), Y = runif(n, 0, 400), lab = sample(1:3, n, replace = TRUE))
  full <- clq_one_image(d, types = c("A", "B", "C"), ks = c(10, 20), perms = 30, seed = 3, keep_unknown = FALSE,
                        script = script)
  part <- clq_one_image(d, types = c("A", "B", "C"), ks = c(10, 20), perms = 30, seed = 3, keep_unknown = FALSE,
                        script = script, perm_ks = 20)
  expect_equal(part$clq, full$clq)
  expect_true(all(is.na(part[k == 10, perm_mean])))
  expect_equal(part[k == 20, perm_mean], full[k == 20, perm_mean])
})
