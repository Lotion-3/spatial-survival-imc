# Part IV, stage 3: patient-level associations between subtype / genotype / node status and CLQ,
# exactly as pre-registered in HYPOTHESES_GENOTYPE.md (frozen at commit e6e9389).
#
#   Rscript scripts/genotype_models.R             # real labels  -> results/genotype/models/
#   Rscript scripts/genotype_models.R --dry-run   # every patient label shuffled (seed 1) -> models_dryrun/
#
# The dry run validates the code end to end without looking at any real association.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
gen <- file.path(root, "results", "genotype")
feat_dir <- file.path(gen, "features")
DRY <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(gen, if (DRY) "models_dryrun" else "models")
dir.create(out, showWarnings = FALSE)

MIN_TYPED <- 500      # repo METABRIC MIN_CELLS
MIN_PATIENTS <- 50
K_PRIMARY <- 20
PERM_SEED <- 20261020L
N_PERM <- 1000L
PAM50_LEVELS <- c("LumA", "LumB", "Her2", "Basal", "claudin-low")
SURR_LEVELS <- c("ER+/HER2-", "HER2+", "TN")
IMMUNE <- c("T cell", "CD8 T cell", "CD4 T cell", "Treg", "B cell", "Macrophage")
seeds_log <- list()

# ------------------------------------------------------------------ data

read_patients <- function(cohort) {
  p <- fread(file.path(gen, paste0("patients_", cohort, ".csv")), colClasses = list(character = "PID"),
             na.strings = c("", "NA"))  # pandas writes missing values as empty fields
  if (DRY) {  # shuffle every patient-level label jointly across patients (keeps their correlations)
    set.seed(1)
    lab <- intersect(c("pam50", "pam50_raw", "intclust", "ER", "PR", "HER2", "TP53", "PIK3CA", "BRCA1", "CASP8",
                       "sequenced", "node_positive", "pN", "surrogate", "clinical_type_published"), names(p))
    idx <- sample(nrow(p))
    p[, (lab) := .SD[idx], .SDcols = lab]
    seeds_log[[length(seeds_log) + 1]] <<- data.table(step = paste0("dry-run label shuffle ", cohort), seed = 1L)
  }
  p
}

# Patient-level analysis table for one CLQ run: one row per PID x k x pair.
load_run <- function(run, patients) {
  comp <- fread(file.path(feat_dir, paste0(run, "_composition.csv")), colClasses = list(character = "PID"))
  clr_cols <- grep("^clr_", names(comp), value = TRUE)
  keep_clr <- head(clr_cols, -1)  # CLR sums to zero: drop the last component
  cmp <- comp[, c("PID", "n_typed", "log_n_typed", keep_clr, grep("^count_", names(comp), value = TRUE)), with = FALSE]
  setnames(cmp, keep_clr, paste0("clr", seq_along(keep_clr)))
  clq <- fread(file.path(feat_dir, paste0(run, "_clq.csv.gz")), colClasses = list(character = "PID"))
  d <- merge(clq, cmp, by = "PID")
  d <- merge(d, patients, by = "PID")
  d <- d[n_typed >= MIN_TYPED & !is.na(age) & !is.na(grade)]
  attr(d, "comp") <- merge(cmp, patients, by = "PID")[n_typed >= MIN_TYPED & !is.na(age) & !is.na(grade)]
  d
}

covars <- function(d, adjusted) c(if (adjusted) grep("^clr[0-9]+$", names(d), value = TRUE), "log_n_typed", "age", "grade")

pair_label <- function(a, b) paste0(a, " -> ", b)

# ------------------------------------------------------------------ model helpers

coef_rows <- function(m, pattern, extra = list()) {
  s <- summary(m)$coefficients
  rows <- grep(pattern, rownames(s), value = TRUE)
  if (!length(rows)) return(NULL)
  ci <- suppressMessages(confint(m)[rows, , drop = FALSE])
  data.table(term = rows, estimate = s[rows, 1], se = s[rows, 2], lo = ci[, 1], hi = ci[, 2], p = s[rows, 4],
             n = nobs(m), as.data.table(extra))
}

fit_group <- function(df, group, cov, family) {
  f_full <- reformulate(c(group, cov), "clq_k")
  m <- lm(f_full, data = df)
  m
}

# ------------------------------------------------------------------ Q1: subtype (F-test of the PAM50 block)

q1_family <- function(d, label = "") {
  d <- d[!is.na(pam50)]
  d[, pam50 := factor(pam50, PAM50_LEVELS)]
  splits <- split(d, by = c("k", "cell_A", "cell_B"), keep.by = TRUE)
  tests <- list(); contrasts <- list()
  for (s in splits) {
    if (nrow(s) < MIN_PATIENTS) next
    for (adj in c(TRUE, FALSE)) {
      cv <- covars(d, adj)
      full <- lm(reformulate(c("pam50", cv), "clq_k"), data = s)
      red <- lm(reformulate(cv, "clq_k"), data = s)
      a <- anova(red, full)
      key <- list(k = s$k[1], cell_A = s$cell_A[1], cell_B = s$cell_B[1], adjusted = adj)
      tests[[length(tests) + 1]] <- data.table(as.data.table(key), n = nrow(s), F = a$F[2], df1 = a$Df[2],
                                               p = a$`Pr(>F)`[2])
      contrasts[[length(contrasts) + 1]] <- coef_rows(full, "^pam50", key)
    }
  }
  tests <- rbindlist(tests); contrasts <- rbindlist(contrasts)
  tests[, q := p.adjust(p, "BH"), by = adjusted]
  contrasts[, contrast := paste0(sub("^pam50", "", term), " - LumA")]
  list(tests = tests, contrasts = contrasts)
}

# ------------------------------------------------------------------ Q2: one binary group, PAM50 (+ other) as covariate

binary_family <- function(d, var, subset_expr = NULL, add_subtype = TRUE) {
  d <- d[!is.na(get(var))]
  if (!is.null(subset_expr)) d <- d[eval(subset_expr)]
  d[, pam50_cov := factor(fifelse(is.na(pam50), "other", pam50), c(PAM50_LEVELS, "other"))]
  d[, grp := as.numeric(get(var))]
  splits <- split(d, by = c("k", "cell_A", "cell_B"), keep.by = TRUE)
  res <- list()
  for (s in splits) {
    if (nrow(s) < MIN_PATIENTS || length(unique(s$grp)) < 2) next
    for (adj in c(TRUE, FALSE)) {
      cv <- c(if (add_subtype && length(unique(s$pam50_cov)) > 1) "pam50_cov", covars(d, adj))
      m <- lm(reformulate(c("grp", cv), "clq_k"), data = s)
      key <- list(k = s$k[1], cell_A = s$cell_A[1], cell_B = s$cell_B[1], adjusted = adj,
                  n_group1 = sum(s$grp == 1), n_group0 = sum(s$grp == 0))
      res[[length(res) + 1]] <- coef_rows(m, "^grp$", key)
    }
  }
  res <- rbindlist(res)
  if (!nrow(res)) return(data.table())  # no pair with >= MIN_PATIENTS and both groups present
  res[, q := p.adjust(p, "BH"), by = adjusted]
  res[, variable := var]
  res
}

# ------------------------------------------------------------------ permutation calibration of Q1

perm_null <- function(d) {
  d <- d[!is.na(pam50)]
  pids <- sort(unique(d$PID))
  lab <- unique(d[, .(PID, pam50)])[match(pids, PID), pam50]
  splits <- split(d, by = c("k", "cell_A", "cell_B"), keep.by = TRUE)
  splits <- splits[sapply(splits, nrow) >= MIN_PATIENTS]
  cv <- covars(d, TRUE)
  pre <- lapply(splits, function(s) {
    Xr <- model.matrix(reformulate(cv), data = s)
    list(y = s$clq_k, Xr = Xr, idx = match(s$PID, pids), rss_r = sum(lm.fit(Xr, s$clq_k)$residuals^2),
         rank_r = qr(Xr)$rank)
  })
  ftest <- function(pr, g) {
    G <- model.matrix(~ factor(g, PAM50_LEVELS))[, -1, drop = FALSE]
    X <- cbind(pr$Xr, G)
    f <- lm.fit(X, pr$y)
    rss <- sum(f$residuals^2); rk <- f$rank
    df1 <- rk - pr$rank_r; df2 <- length(pr$y) - rk
    Fv <- ((pr$rss_r - rss) / df1) / (rss / df2)
    c(F = Fv, p = pf(Fv, df1, df2, lower.tail = FALSE))
  }
  obs <- t(sapply(pre, function(pr) ftest(pr, lab[pr$idx])))
  set.seed(PERM_SEED)
  seeds_log[[length(seeds_log) + 1]] <<- data.table(step = "Q1 subtype-label permutations (1,000)", seed = PERM_SEED)
  perms <- replicate(N_PERM, sample(length(pids)))
  Fp <- matrix(NA_real_, N_PERM, length(pre)); Pp <- Fp
  for (r in seq_len(N_PERM)) {
    lr <- lab[perms[, r]]
    v <- sapply(pre, function(pr) ftest(pr, lr[pr$idx]))
    Fp[r, ] <- v["F", ]; Pp[r, ] <- v["p", ]
  }
  per_perm <- data.table(perm = seq_len(N_PERM), n_discoveries = apply(Pp, 1, function(p) sum(p.adjust(p, "BH") < 0.05)),
                         min_p = apply(Pp, 1, min), mean_p = rowMeans(Pp))
  keys <- rbindlist(lapply(splits, function(s) data.table(k = s$k[1], cell_A = s$cell_A[1], cell_B = s$cell_B[1])))
  per_test <- cbind(keys, F_obs = obs[, "F"], p_param = obs[, "p"],
                    p_perm = (1 + colSums(Fp >= matrix(obs[, "F"], N_PERM, length(pre), byrow = TRUE))) / (1 + N_PERM))
  summ <- list(n_perm = N_PERM, seed = PERM_SEED, n_tests = length(pre),
               share_perms_any_discovery = mean(per_perm$n_discoveries > 0),
               mean_null_p = mean(Pp), share_null_p_lt_005 = mean(Pp < 0.05),
               calibrated = mean(per_perm$n_discoveries > 0) <= 0.10 && abs(mean(Pp) - 0.5) <= 0.02)
  list(per_perm = per_perm, per_test = per_test, summary = summ)
}

# ------------------------------------------------------------------ main

main <- function() {
  pm <- read_patients("metabric")
  pb <- read_patients("basel")
  met <- load_run("metabric_celesta", pm)
  met1 <- load_run("metabric_round1", pm)
  bas <- load_run("basel_celesta", pb)
  met_pub <- load_run("metabric_published_k20", pm)
  bas_pub <- load_run("basel_published_k20", pb)
  met_unk <- load_run("metabric_celesta_keepunknown_k20", pm)

  # ---- Q1
  q1 <- q1_family(met)
  fwrite(q1$tests, file.path(out, "q1_subtype_tests.csv"))
  fwrite(q1$contrasts, file.path(out, "q1_subtype_contrasts.csv"))

  # ---- Q2 (TP53, PIK3CA) + within-subtype TP53 (descriptive)
  q2t <- binary_family(met, "TP53")
  q2p <- binary_family(met, "PIK3CA")
  fwrite(q2t, file.path(out, "q2_tp53.csv"))
  fwrite(q2p, file.path(out, "q2_pik3ca.csv"))
  within <- rbindlist(lapply(c("LumA", "LumB", "Her2", "claudin-low"), function(st) {
    r <- binary_family(met[k == K_PRIMARY], "TP53", subset_expr = bquote(pam50 == .(st)), add_subtype = FALSE)
    if (nrow(r)) r[, subtype := st]
    r
  }), fill = TRUE)
  fwrite(within, file.path(out, "q2_tp53_within_subtype.csv"))

  # ---- H1, H2, (d), P2, P3 from the families / dedicated fits
  get_c <- function(tab, a, b, kk, ctr, adj = TRUE) tab[cell_A == a & cell_B == b & k == kk & adjusted == adj &
                                                          contrast == ctr]
  h1 <- get_c(q1$contrasts, "Tumour", "CD8 T cell", K_PRIMARY, "Basal - LumA")
  h1rev <- get_c(q1$contrasts, "CD8 T cell", "Tumour", K_PRIMARY, "Basal - LumA")
  h2 <- binary_family(met[k == K_PRIMARY & cell_A == "Tumour" & cell_B == "Macrophage"], "TP53",
                      subset_expr = quote(ER == 1))
  node_fit <- function(d, label) {
    r <- binary_family(d[k == K_PRIMARY & ((cell_A == "Tumour" & cell_B == "Treg") | (cell_A == "Treg" & cell_B == "Tumour"))],
                       "node_positive")
    if (nrow(r)) r[, labels := label]
    r
  }
  qd <- rbind(node_fit(met, "CELESTA"), node_fit(met_pub, "published"), fill = TRUE)
  fwrite(qd, file.path(out, "question_d_node_treg.csv"))

  q1_pub <- q1_family(met_pub)
  q1_r1 <- q1_family(met1)
  fwrite(q1_pub$tests, file.path(out, "alt_published_q1_tests.csv"))
  fwrite(q1_pub$contrasts, file.path(out, "alt_published_q1_contrasts.csv"))
  fwrite(q1_r1$tests, file.path(out, "round1_q1_tests.csv"))
  fwrite(q1_r1$contrasts, file.path(out, "round1_q1_contrasts.csv"))

  p2 <- rbind(
    get_c(q1$contrasts, "Tumour", "Tumour", K_PRIMARY, "Basal - LumA")[, control := "P2a"],
    get_c(q1$contrasts, "Stroma", "Stroma", K_PRIMARY, "Basal - LumA")[, control := "P2a"],
    get_c(q1_r1$contrasts, "B cell", "T cell", K_PRIMARY, "Basal - LumA")[, control := "P2b"],
    get_c(q1_r1$contrasts, "T cell", "B cell", K_PRIMARY, "Basal - LumA")[, control := "P2b"],
    get_c(q1$contrasts, "Treg", "CD8 T cell", K_PRIMARY, "Basal - LumA")[, control := "P2c"],
    get_c(q1$contrasts, "Treg", "CD4 T cell", K_PRIMARY, "Basal - LumA")[, control := "P2c"],
    get_c(q1_pub$contrasts, "Treg", "CD8 T cell", K_PRIMARY, "Basal - LumA")[, control := "P2c (published labels)"],
    get_c(q1_pub$contrasts, "Treg", "CD4 T cell", K_PRIMARY, "Basal - LumA")[, control := "P2c (published labels)"],
    fill = TRUE)
  fwrite(p2, file.path(out, "p2_controls.csv"))

  p3 <- rbind(
    binary_family(met[k == K_PRIMARY & cell_A == "Treg" & cell_B %in% c("CD8 T cell", "CD4 T cell")], "BRCA1")[, labels := "CELESTA"],
    binary_family(met_pub[k == K_PRIMARY & cell_A == "Treg" & cell_B %in% c("CD8 T cell", "CD4 T cell")], "BRCA1")[, labels := "published"],
    fill = TRUE)
  p3_desc <- rbind(
    met[k == K_PRIMARY & cell_A == "Treg" & cell_B %in% c("CD8 T cell", "CD4 T cell") & !is.na(BRCA1),
        .(labels = "CELESTA", n = .N, mean_clq = mean(clq_k), sd = sd(clq_k)), by = .(cell_A, cell_B, BRCA1)],
    met_pub[k == K_PRIMARY & cell_A == "Treg" & cell_B %in% c("CD8 T cell", "CD4 T cell") & !is.na(BRCA1),
            .(labels = "published", n = .N, mean_clq = mean(clq_k), sd = sd(clq_k)), by = .(cell_A, cell_B, BRCA1)])
  fwrite(p3, file.path(out, "p3_brca1.csv"))
  fwrite(p3_desc, file.path(out, "p3_brca1_descriptive.csv"))

  # ---- P1 composition (METABRIC; Basel replication with surrogate subtypes)
  comp_m <- attr(met, "comp")[!is.na(pam50)]
  imm_cols <- paste0("count_", intersect(IMMUNE, sub("^count_", "", grep("^count_", names(comp_m), value = TRUE))))
  comp_m[, immune_fraction := rowSums(.SD) / n_typed, .SDcols = imm_cols]
  comp_m[, pam50 := factor(pam50, PAM50_LEVELS)]
  p1m <- coef_rows(lm(immune_fraction ~ pam50 + log_n_typed + age + grade, data = comp_m), "^pam50",
                   list(cohort = "METABRIC", labels = "CELESTA"))
  comp_mp <- attr(met_pub, "comp")[!is.na(pam50)]
  comp_mp[, immune_fraction := rowSums(.SD) / n_typed, .SDcols = imm_cols]
  comp_mp[, pam50 := factor(pam50, PAM50_LEVELS)]
  p1mp <- coef_rows(lm(immune_fraction ~ pam50 + log_n_typed + age + grade, data = comp_mp), "^pam50",
                    list(cohort = "METABRIC", labels = "published"))
  comp_b <- attr(bas, "comp")[!is.na(surrogate)]
  imm_b <- paste0("count_", intersect(IMMUNE, sub("^count_", "", grep("^count_", names(comp_b), value = TRUE))))
  comp_b[, immune_fraction := rowSums(.SD) / n_typed, .SDcols = imm_b]
  comp_b[, surrogate := factor(surrogate, SURR_LEVELS)]
  p1b <- coef_rows(lm(immune_fraction ~ surrogate + log_n_typed + age + grade, data = comp_b), "^surrogate",
                   list(cohort = "Basel", labels = "CELESTA"))
  p1 <- rbind(p1m, p1mp, p1b)
  fwrite(p1, file.path(out, "p1_immune_fraction.csv"))
  fwrite(rbind(comp_m[, .(cohort = "METABRIC", group = pam50, immune_fraction)],
               comp_b[, .(cohort = "Basel", group = surrogate, immune_fraction)])[
    , .(n = .N, mean = mean(immune_fraction), median = median(immune_fraction)), by = .(cohort, group)],
    file.path(out, "p1_immune_fraction_by_group.csv"))

  # ---- H3 replication: METABRIC round-1 Basal / Her2 contrasts -> Basel TN / HER2+
  rep_m <- q1_r1$contrasts[adjusted == TRUE & contrast %in% c("Basal - LumA", "Her2 - LumA")]
  rep_m[, q := p.adjust(p, "BH")]
  fwrite(rep_m, file.path(out, "replication_metabric_family.csv"))
  sig <- rep_m[q < 0.05]
  bas_s <- bas[!is.na(surrogate)]
  bas_s[, surrogate := factor(surrogate, SURR_LEVELS)]
  rep <- rbindlist(lapply(seq_len(nrow(sig)), function(i) {
    r <- sig[i]
    s <- bas_s[k == r$k & cell_A == r$cell_A & cell_B == r$cell_B]
    if (nrow(s) < MIN_PATIENTS) return(NULL)
    m <- lm(reformulate(c("surrogate", covars(bas, TRUE)), "clq_k"), data = s)
    term <- if (r$contrast == "Basal - LumA") "surrogateTN" else "surrogateHER2+"
    b <- coef_rows(m, paste0("^", gsub("\\+", "\\\\+", term), "$"))
    data.table(k = r$k, cell_A = r$cell_A, cell_B = r$cell_B, metabric_contrast = r$contrast,
               metabric_estimate = r$estimate, metabric_lo = r$lo, metabric_hi = r$hi, metabric_q = r$q,
               basel_contrast = sub("surrogate", "", term), basel_estimate = b$estimate, basel_lo = b$lo,
               basel_hi = b$hi, basel_p = b$p, basel_n = b$n)
  }))
  if (nrow(rep)) rep[, same_sign := sign(metabric_estimate) == sign(basel_estimate)]
  fwrite(rep, file.path(out, "replication.csv"))
  # all METABRIC round-1 contrasts with their Basel counterpart (for the scatter plot)
  allrep <- rbindlist(lapply(split(rep_m, by = c("k", "cell_A", "cell_B", "contrast")), function(r) {
    s <- bas_s[k == r$k & cell_A == r$cell_A & cell_B == r$cell_B]
    if (nrow(s) < MIN_PATIENTS) return(NULL)
    m <- lm(reformulate(c("surrogate", covars(bas, TRUE)), "clq_k"), data = s)
    term <- if (r$contrast == "Basal - LumA") "^surrogateTN$" else "^surrogateHER2\\+$"
    b <- coef_rows(m, term)
    data.table(k = r$k, cell_A = r$cell_A, cell_B = r$cell_B, contrast = r$contrast, metabric_estimate = r$estimate,
               metabric_q = r$q, basel_estimate = b$estimate, basel_p = b$p)
  }))
  fwrite(allrep, file.path(out, "replication_all_contrasts.csv"))

  # ---- robustness: permutation calibration, leave-one-k-out, alternative typing, keep-unknown
  pn <- perm_null(met)
  fwrite(pn$per_perm, file.path(out, "permutation_null_per_perm.csv"))
  fwrite(pn$per_test, file.path(out, "permutation_null_per_test.csv"))
  write_json(pn$summary, file.path(out, "permutation_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)

  disc <- q1$tests[adjusted == TRUE & q < 0.05]
  loko <- rbindlist(lapply(seq_len(nrow(disc)), function(i) {
    r <- disc[i]
    o <- q1$tests[adjusted == TRUE & cell_A == r$cell_A & cell_B == r$cell_B & k != r$k]
    data.table(k = r$k, cell_A = r$cell_A, cell_B = r$cell_B, q = r$q, other_k = paste(o$k, collapse = ","),
               other_p = paste(signif(o$p, 3), collapse = ","), consistent = all(o$p < 0.05) && nrow(o) == 2)
  }))
  fwrite(loko, file.path(out, "leave_one_k_out.csv"))

  alt <- merge(q1$tests[adjusted == TRUE & k == K_PRIMARY, .(cell_A, cell_B, F_celesta = F, p_celesta = p, q_celesta = q)],
               q1_pub$tests[adjusted == TRUE, .(cell_A, cell_B, F_published = F, p_published = p, q_published = q)],
               by = c("cell_A", "cell_B"))
  bl <- merge(q1$contrasts[adjusted == TRUE & k == K_PRIMARY & contrast == "Basal - LumA", .(cell_A, cell_B, est_celesta = estimate)],
              q1_pub$contrasts[adjusted == TRUE & contrast == "Basal - LumA", .(cell_A, cell_B, est_published = estimate)],
              by = c("cell_A", "cell_B"))
  alt <- merge(alt, bl, by = c("cell_A", "cell_B"), all = TRUE)
  fwrite(alt, file.path(out, "alt_typing_comparison.csv"))
  q1_unk <- q1_family(met_unk)
  fwrite(q1_unk$tests, file.path(out, "sensitivity_keepunknown_q1_tests.csv"))
  fwrite(q1_unk$contrasts, file.path(out, "sensitivity_keepunknown_q1_contrasts.csv"))

  # ---- descriptive (Keren, Jackson)
  desc_m <- met1[k == K_PRIMARY & !is.na(pam50), .(cohort = "METABRIC", group = pam50, PID, cell_A, cell_B, clq_k)]
  desc_b <- bas[k == K_PRIMARY & !is.na(surrogate), .(cohort = "Basel", group = surrogate, PID, cell_A, cell_B, clq_k)]
  desc <- rbind(desc_m, desc_b)[(cell_A == "Tumour" & cell_B %in% c("T cell", "Macrophage", "B cell")) |
                                  (cell_A == "T cell" & cell_B == "Endothelial")]
  fwrite(desc[, .(n = .N, median = median(clq_k), q25 = quantile(clq_k, 0.25), q75 = quantile(clq_k, 0.75),
                  iqr = IQR(clq_k)), by = .(cohort, group, cell_A, cell_B)], file.path(out, "descriptive_keren_jackson.csv"))

  # ---- hypothesis / control outcomes
  pass <- function(x) if (length(x) == 0 || is.na(x)) "not estimable" else if (x) "pass" else "fail"
  # positive and p < 0.05; "not estimable" when the model could not be fitted
  pos_sig <- function(r) if (!nrow(r) || is.na(r$estimate[1])) "not estimable" else pass(r$estimate[1] > 0 && r$p[1] < 0.05)
  all_pos <- function(x) if (!length(x) || anyNA(x)) "not estimable" else pass(all(x > 0))
  na1 <- function(x) if (length(x)) x[1] else NA_real_
  h2a <- h2[adjusted == TRUE]
  qd_c <- qd[adjusted == TRUE & labels == "CELESTA" & cell_A == "Tumour" & cell_B == "Treg"]
  n_rep <- nrow(rep); n_same <- if (n_rep) sum(rep$same_sign) else 0
  p1_get <- function(tm) p1m[term == tm]
  outcomes <- rbind(
    data.table(id = "H1", test = "METABRIC Tumour -> CD8 T cell, k = 20, Basal - LumA (adjusted)",
               estimate = na1(h1$estimate), lo = na1(h1$lo), hi = na1(h1$hi), p = na1(h1$p), n = na1(h1$n),
               outcome = pos_sig(h1)),
    data.table(id = "H1-rev", test = "METABRIC CD8 T cell -> Tumour, k = 20, Basal - LumA (adjusted)",
               estimate = na1(h1rev$estimate), lo = na1(h1rev$lo), hi = na1(h1rev$hi), p = na1(h1rev$p), n = na1(h1rev$n),
               outcome = pos_sig(h1rev)),
    data.table(id = "H2", test = "METABRIC ER+, Tumour -> Macrophage, k = 20, TP53 mutant - wild-type (adjusted)",
               estimate = na1(h2a$estimate), lo = na1(h2a$lo), hi = na1(h2a$hi), p = na1(h2a$p), n = na1(h2a$n),
               outcome = pos_sig(h2a)),
    data.table(id = "H3", test = sprintf("METABRIC round-1 contrasts with q < 0.05 (%d) with the same sign in Basel", n_rep),
               estimate = n_same, lo = NA, hi = NA, p = NA, n = n_rep, outcome = pass(n_same >= 1)),
    data.table(id = "H3-strict", test = "share of replicated signs > 0.5 (one-sided binomial; needs >= 5 discoveries)",
               estimate = if (n_rep) n_same / n_rep else NA, lo = NA, hi = NA,
               p = if (n_rep) binom.test(n_same, n_rep, 0.5, alternative = "greater")$p.value else NA, n = n_rep,
               outcome = if (n_rep < 5) "not estimable" else pass(binom.test(n_same, n_rep, 0.5, alternative = "greater")$p.value < 0.05)),
    data.table(id = "P1", test = "METABRIC immune fraction: Basal - LumA and Her2 - LumA",
               estimate = na1(p1_get("pam50Basal")$estimate), lo = na1(p1_get("pam50Basal")$lo),
               hi = na1(p1_get("pam50Basal")$hi), p = na1(p1_get("pam50Basal")$p), n = na1(p1_get("pam50Basal")$n),
               outcome = all_pos(c(p1_get("pam50Basal")$lo, p1_get("pam50Her2")$lo))),
    data.table(id = "P2a", test = "Tumour -> Tumour and Stroma -> Stroma, Basal - LumA > 0",
               estimate = NA, lo = NA, hi = NA, p = NA, n = NA, outcome = all_pos(p2[control == "P2a", estimate])),
    data.table(id = "P2b", test = "B -> T and T -> B (round-1), Basal - LumA > 0",
               estimate = NA, lo = NA, hi = NA, p = NA, n = NA, outcome = all_pos(p2[control == "P2b", estimate])),
    data.table(id = "P2c", test = "Treg -> CD8 T and Treg -> CD4 T, Basal - LumA > 0 (CELESTA)",
               estimate = NA, lo = NA, hi = NA, p = NA, n = NA, outcome = all_pos(p2[control == "P2c", estimate])),
    data.table(id = "d", test = "METABRIC Tumour -> Treg, k = 20, node-positive - node-negative (adjusted, CELESTA)",
               estimate = na1(qd_c$estimate), lo = na1(qd_c$lo), hi = na1(qd_c$hi), p = na1(qd_c$p), n = na1(qd_c$n),
               outcome = pos_sig(qd_c)),
    fill = TRUE)
  p2_n <- sum(outcomes[id %in% c("P2a", "P2b", "P2c"), outcome == "pass"])
  outcomes <- rbind(outcomes, data.table(id = "P2", test = "at least 2 of P2a-c pass", estimate = p2_n,
                                         outcome = pass(p2_n >= 2)), fill = TRUE)
  fwrite(outcomes, file.path(out, "hypothesis_outcomes.csv"))
  fwrite(rbindlist(seeds_log), file.path(out, "seeds.csv"))
  print(outcomes[, .(id, estimate = signif(estimate, 3), lo = signif(lo, 3), hi = signif(hi, 3), p = signif(p, 3), n, outcome)])
  str(pn$summary)
  cat(sprintf("Q1 discoveries (adjusted, q < 0.05): %d of %d; unadjusted: %d\n", q1$tests[adjusted == TRUE, sum(q < 0.05)],
              q1$tests[adjusted == TRUE, .N], q1$tests[adjusted == FALSE, sum(q < 0.05)]))
  cat(sprintf("TP53 discoveries: %d; PIK3CA: %d\n", q2t[adjusted == TRUE, sum(q < 0.05)], q2p[adjusted == TRUE, sum(q < 0.05)]))
}

if (sys.nframe() == 0) main()
