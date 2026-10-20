# Part IV, stage 0: install the R side and record the environment.
#
#   Rscript scripts/setup_r.R
#
# - CRAN dependencies of CELESTA and of the Plevritis lab CLQ script, plus data.table/jsonlite/testthat.
# - The official CELESTA package at a pinned commit (plevritis-lab/CELESTA; the paper's
#   github.com/plevritis/CELESTA redirects to the same repository and HEAD).
# - Both lab repositories cloned unmodified into third_party/ at pinned commits (git-ignored), so
#   their source can be read and the CLQ functions sourced from the original file.
# - results/genotype/environment.json: R version, package versions, commit SHAs.

CELESTA_SHA <- "e5e92ffcef5ccada8053c2e042f34c446b6b3565"
CLQ_SHA <- "02b6da5e7bb8ebabb68a05daa78ae46c36ec0378"
CRAN <- c("Rmixmod", "spdep", "zeallot", "reshape2", "ggplot2", "dplyr", "tidyr", "data.table", "R.utils", "jsonlite",
          "testthat", "remotes")

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
missing <- setdiff(CRAN, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("CELESTA", quietly = TRUE) || packageDescription("CELESTA")$RemoteSha != CELESTA_SHA) {
  remotes::install_github(paste0("plevritis-lab/CELESTA@", CELESTA_SHA), upgrade = "never", dependencies = FALSE)
}

clone <- function(repo, sha) {
  dest <- file.path(root, "third_party", basename(repo))
  if (!dir.exists(dest)) system2("git", c("clone", "-q", paste0("https://github.com/", repo, ".git"), shQuote(dest)))
  system2("git", c("-C", shQuote(dest), "checkout", "-q", sha))
  head <- system2("git", c("-C", shQuote(dest), "rev-parse", "HEAD"), stdout = TRUE)
  dirty <- system2("git", c("-C", shQuote(dest), "status", "--porcelain"), stdout = TRUE)
  stopifnot(head == sha, length(dirty) == 0)  # unmodified original code
  date <- system2("git", c("-C", shQuote(dest), "log", "-1", "--format=%cI"), stdout = TRUE)
  list(url = paste0("https://github.com/", repo), commit = head, commit_date = date, path = file.path("third_party", basename(repo)))
}
repos <- list(CELESTA = clone("plevritis-lab/CELESTA", CELESTA_SHA),
              Spatial_Permutation_and_Normalization = clone("plevritis-lab/Spatial_Permutation_and_Normalization", CLQ_SHA))

pkgs <- c("CELESTA", CRAN)
env <- list(
  R_version = R.version.string,
  platform = R.version$platform,
  os = utils::sessionInfo()$running,
  packages = setNames(lapply(pkgs, function(p) as.character(packageVersion(p))), pkgs),
  CELESTA_install = list(remote = "plevritis-lab/CELESTA", remote_sha = packageDescription("CELESTA")$RemoteSha,
                         paper_url_resolves_to_same_head = TRUE),
  repositories = repos,
  recorded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
out <- file.path(root, "results", "genotype")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(env, file.path(out, "environment.json"), auto_unbox = TRUE, pretty = TRUE)
cat("wrote", file.path(out, "environment.json"), "\n")
