# =============================================================================
# exactness_ties.R   [IMPROVEMENT-5]  Is the fast scoring EXACT? What about ties?
#
# Part A. TCPD (42 datasets x 2 costs): the fast scores are compared with the
#   scores of the original code (saved by the baseline in
#   project0_reproduction/output-data, read only), for several values of
#   exact_max_len (L0): segments with <= L0 points use the bit-identical replica of
#   the original cost function, longer ones the O(1) prefix-sum formula.
#     L0 = 0   : pure prefix sums (mathematically the same algorithm)
#     L0 = 64  : default of the package
#     L0 = Inf : replica costs everywhere (bit-identical scores expected)
#   For each L0: max |score difference|, number of bit-identical score vectors, and
#   number of identical selections among the 24024 (dataset, cost, threshold, level)
#   rows of the pipeline (Algorithm 2 run on the fast scores with the exact fast
#   selection path, compared with the baseline selections).
# Part B. Small examples that show WHY pure prefix sums are not enough: in the
#   original code, gains that are mathematically equal (or equal to 0) are compared
#   through their floating-point rounding noise.
#
# Part C. Time of one segment-cost call (original / replica / prefix sums).
# Part D. Longest constant and collinear runs of the TCPD series (where exact ties occur).
#
# Usage (from the project folder):  Rscript experiments/exactness_ties.R
# Outputs: results/exactness_by_exact_max_len.csv, results/exactness_details.csv,
#          results/tie_examples.txt, results/cost_call_microbenchmark.csv,
#          results/tcpd_constant_and_collinear_runs.csv
# Run time: about 1-3 minutes.
# =============================================================================
suppressMessages(library(dplyr))
stopifnot(requireNamespace("scoth.fast", quietly = TRUE))
ns = asNamespace("scoth.fast")

RESULTS = file.path(here::here(), "results")
dir.create(RESULTS, showWarnings = FALSE)
BASE = file.path(dirname(here::here()), "project0_reproduction", "output-data")   # baseline (read only)
stopifnot(dir.exists(file.path(BASE, "algorithm_selection_results")))
thresholds = c(seq(from = 0.03, to = 0.2, by = 0.01), 0.3, 0.4, 0.5, 0.6, 0.7, 1)
L0_GRID = c(0, 2, 8, 64, Inf)

## Part A -----------------------------------------------------------------------
options(scoth.fast_selection = TRUE)
files = list.files(file.path(BASE, "algorithm_scoring_results"))
details = list()
for (L0 in L0_GRID) {
  for (f in files) {
    d = readRDS(file.path(BASE, "algorithm_scoring_results", f))
    sb = readRDS(file.path(BASE, "algorithm_selection_results", f))
    for (r in seq_len(nrow(d))) {
      x = d$x[[r]]
      cn = d$cost_method[r]
      cn2 = if (is.data.frame(x)) paste0(cn, "_multivariate") else cn
      t0 = Sys.time()
      sc = ns$scoth_changes_fast(x, cn2, exact_max_len = L0)
      sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
      base_score = d$scoring[[r]]$score
      scoring = data.frame(index = 2:(length(sc) + 1), score = sc)
      invisible(capture.output(
        sel <- ns$pipeline_selection_single_scoring(scoring, thresholds, 10, d$cost_initial[r],
                                                    ns$all_cost_functions()[[cn2]], x)))
      base = sb %>% filter(cost_method == cn)
      stopifnot(identical(base$selection_threshold_included, sel$selection_threshold_included),
                identical(base$selection_level, sel$selection_level))
      same = mapply(function(a, b) identical(a$predictions, b$predictions), base$selection, sel$selection)
      details[[length(details) + 1]] = data.frame(
        exact_max_len = L0, dataset = d$dataset_name[r], cost = cn, T = length(sc) + 1,
        max_abs_score_diff = max(abs(sc - base_score)),
        n_scores_diff_gt_1e9 = sum(abs(sc - base_score) > 1e-9),
        bit_identical = identical(sc, base_score),
        n_selection_rows = length(same), n_identical_selection_rows = sum(same),
        scoring_seconds = sec)
    }
  }
  message(sprintf("exact_max_len = %s done", L0))
}
details = bind_rows(details)
write.csv(details, file.path(RESULTS, "exactness_details.csv"), row.names = FALSE)

by_L0 = details %>% group_by(exact_max_len) %>%
  summarize(n_runs = n(),
            n_runs_bit_identical = sum(bit_identical),
            n_runs_max_diff_le_1e12 = sum(max_abs_score_diff <= 1e-12),
            n_runs_all_selections_identical = sum(n_identical_selection_rows == n_selection_rows),
            # (the next columns re-use the names of per-run columns: they come after the counts)
            max_abs_score_diff = max(max_abs_score_diff),
            n_selection_rows = sum(n_selection_rows),
            n_identical_selection_rows = sum(n_identical_selection_rows),
            pct_identical_selection_rows = round(100 * n_identical_selection_rows / n_selection_rows, 2),
            total_scoring_seconds = round(sum(scoring_seconds), 2), .groups = "drop")
write.csv(by_L0, file.path(RESULTS, "exactness_by_exact_max_len.csv"), row.names = FALSE)
print(as.data.frame(by_L0))

## Part B -----------------------------------------------------------------------
# The original loop of scoth_changes(), with the order of the removed cuts recorded.
original_with_order = function(x, cost_func) {
  e = length(x) + 1
  k = e - 2
  cost_initial = cost_func(x, 1, e)
  l = ns$build_df_df2(x, 1:e, cost_func)
  gain = rep(0, k); removed = integer(0); first_gains = NULL
  for (loop in 1:k) {
    diff_gains = ns$compute_diff_costs(l) / cost_initial
    if (loop == 1) first_gains = diff_gains
    alive = l$df$e - 1
    alive = alive[-length(alive)]
    gain[alive] = pmax(diff_gains, gain[alive])
    i = which.min(gain[alive])
    removed = c(removed, alive[i] + 1)          # time index of the removed cut
    l = ns$update_costs(l, i, x, cost_func)
  }
  list(score = gain, removed = removed, first_gains = first_gains)
}
# the same algorithm where a line through 2 points has cost EXACTLY 0 (exact arithmetic)
c_reg_exact2 = function(x, s, e) if (e - s <= 2) 0 else ns$c_reg_cusum_L2(x, s, e)

out = c("Floating-point ties in the original scoring algorithm (Algorithm 1)", strrep("=", 70), "")
set.seed(1)
x = round(rnorm(12), 2)
a = original_with_order(x, ns$c_reg_cusum_L2)
b = original_with_order(x, c_reg_exact2)
out = c(out,
  "Example 1 - linear cost, x = round(rnorm(12), 2) with set.seed(1):",
  paste("  x =", paste(x, collapse = ", ")),
  "  At the first iteration every gain is the cost of a line through 2 points,",
  "  i.e. exactly 0 in exact arithmetic. The original code (lm) returns rounding noise:",
  paste("  first normalized gains =", paste(sprintf("%.1e", a$first_gains), collapse = ", ")),
  paste("  order of removal (original, decided by the noise):     ", paste(a$removed, collapse = ", ")),
  paste("  order of removal (exact zeros, ties -> smallest index):", paste(b$removed, collapse = ", ")),
  sprintf("  max |score difference| between the two runs = %.3f", max(abs(a$score - b$score))),
  sprintf("  fast algorithm (default exact_max_len = 64) identical to the original: %s",
          identical(ns$scoth_changes_fast(x, "reg_cusum_L2"), a$score)),
  sprintf("  fast algorithm with pure prefix sums (exact_max_len = 0), max |diff| = %.3f",
          max(abs(ns$scoth_changes_fast(x, "reg_cusum_L2", exact_max_len = 0) - a$score))),
  "")

x2 = rep(c(0, 0, 1, 1, 3), 6)
a2 = original_with_order(x2, ns$c_cusum_L2)
out = c(out,
  "Example 2 - quadratic cost, integer series with exact ties: x = rep(c(0,0,1,1,3), 6):",
  paste("  first normalized gains =", paste(sprintf("%.4f", a2$first_gains), collapse = ", ")),
  "  Many gains are EXACTLY equal; which.min() removes the smallest index first.",
  paste("  order of removal (original):", paste(head(a2$removed, 15), collapse = ", "), "..."),
  sprintf("  fast algorithm (heap ordered by (score, index)) identical to the original: %s",
          identical(ns$scoth_changes_fast(x2, "cusum_L2"), a2$score)),
  "")

x3 = rep(2.5, 10)
out = c(out,
  "Example 3 - constant series rep(2.5, 10):",
  sprintf("  quadratic cost: V(empty) = %g (exactly 0, all scores are 0)", ns$c_cusum_L2(x3, 1, 11)),
  sprintf("  linear cost:    V(empty) = %.3e (lm rounding noise, not 0!)", ns$c_reg_cusum_L2(x3, 1, 11)),
  paste("  original linear-cost scores = noise / noise =",
        paste(sprintf("%.3f", ns$scoth_changes(x3, ns$c_reg_cusum_L2, 1, 11, verbose = FALSE)), collapse = ", ")),
  sprintf("  fast algorithm identical to the original: %s",
          identical(ns$scoth_changes_fast(x3, "reg_cusum_L2"), ns$scoth_changes(x3, ns$c_reg_cusum_L2, 1, 11, verbose = FALSE))),
  "",
  "Conclusion: the original results depend on rounding noise whenever gains are tied",
  "(mathematically equal or equal to 0). To reproduce them EXACTLY, short segments",
  "must be computed with the same floating-point operations as the original",
  "(bit-identical replicas, exact_max_len); long segments use O(1) prefix sums.")
writeLines(out, file.path(RESULTS, "tie_examples.txt"))
cat(out, sep = "\n")

## Part C -----------------------------------------------------------------------
# Time of ONE segment-cost call: original function (lm() for the linear cost),
# bit-identical replica, and O(1) prefix-sum formula (x = rnorm(10000), seed 1).
time_per_call = function(f) {   # repeat until the measurement lasts >= 0.5 s (clock resolution ~10 ms)
  n = 100
  repeat {
    sec = system.time(for (i in seq_len(n)) f())[["elapsed"]]
    if (sec >= 0.5 || n >= 1e6) return(c(microseconds = 1e6 * sec / n, n_calls = n))
    n = n * 4
  }
}
set.seed(1)
x = rnorm(10000)
micro = list()
for (cn in c("cusum_L2", "reg_cusum_L2")) {
  oracle = ns$fast_segment_cost(x, cn, exact_max_len = 64)
  orig_f = ns$all_cost_functions()[[cn]]
  for (len in c(3, 30, 300, 3000)) {
    s = 1001; e = s + len
    for (impl in c("original", "replica", "prefix_sum")) {
      f = switch(impl, original = function() orig_f(x, s, e), replica = function() oracle$cost_exact(s, e),
                 prefix_sum = function() oracle$cost_prefix(s, e))
      tm = time_per_call(f)
      micro[[length(micro) + 1]] = data.frame(cost = cn, segment_length = len, implementation = impl,
                                              microseconds_per_call = round(tm[["microseconds"]], 2),
                                              n_calls_timed = tm[["n_calls"]], value = f())
    }
  }
}
micro = bind_rows(micro)
write.csv(micro, file.path(RESULTS, "cost_call_microbenchmark.csv"), row.names = FALSE)
print(micro)

## Part D -----------------------------------------------------------------------
# Where can exact ties occur in TCPD? Longest run of identical rows (constant
# segment: exact 0 gains) and longest exactly collinear run (univariate only:
# zero cost for the linear cost) of each dataset.
kpis = ns$load_kpis(file.path(here::here(), "public-data"))
runs = bind_rows(lapply(seq_len(nrow(kpis)), function(i) {
  M = as.matrix(kpis$kpis[[i]])
  same_row = c(FALSE, rowSums(M[-1, , drop = FALSE] != M[-nrow(M), , drop = FALSE]) == 0)
  const_runs = rle(cumsum(!same_row))$lengths
  collinear = if (ncol(M) == 1) {
    r = rle(abs(diff(M[, 1], differences = 2)) <= 1e-12 * max(1, abs(M[, 1])))  # 2nd differences = 0
    if (any(r$values)) max(r$lengths[r$values]) + 2 else 2
  } else NA
  data.frame(dataset = kpis$dataset_name[i], T = nrow(M), n_features = ncol(M),
             longest_constant_run = max(const_runs), longest_collinear_run = collinear)
}))
write.csv(runs, file.path(RESULTS, "tcpd_constant_and_collinear_runs.csv"), row.names = FALSE)
cat(sprintf("TCPD: longest constant run = %d points (%s); longest collinear run = %d points (%s)\n",
            max(runs$longest_constant_run), paste(runs$dataset[runs$longest_constant_run == max(runs$longest_constant_run)], collapse = ", "),
            max(runs$longest_collinear_run, na.rm = TRUE),
            paste(runs$dataset[which(runs$longest_collinear_run == max(runs$longest_collinear_run, na.rm = TRUE))], collapse = ", ")))
