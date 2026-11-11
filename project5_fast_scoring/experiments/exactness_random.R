# =============================================================================
# exactness_random.R   [IMPROVEMENT-5]  exactness on many random series
#
# Complement of the targeted experiment (targeted_runtime.R) with >= 50
# replications per condition: is the fast scoring exact on random series, also on
# tie-heavy data? For every series, the scores of the ORIGINAL algorithm
# (scoth::scoth_scoring_func, baseline package) are compared with the fast ones
# for exact_max_len = 0 (pure prefix sums), 64 (default) and Inf (replicas only),
# and the selections of Algorithm 2 (24 thresholds of the paper x 10 levels) are
# compared row by row.
#
# Conditions: 4 data types x 2 lengths (T = 100, 400) x 2 costs, 50 series each (seeds 1..50)
#   continuous : piecewise-constant mean (3 change points) + N(0,1) noise
#   integer    : Poisson counts with 3 level shifts (many exactly equal values)
#   decimal    : random walk rounded to 1 decimal (many near ties)
#   const_runs : piecewise-constant mean, noise only on 50% of the points (exact constant runs)
#
# Usage (from the project folder):  Rscript experiments/exactness_random.R
# Outputs: results/exactness_random_series.csv (one row per series and exact_max_len)
#          results/exactness_random_summary.csv
# Run time: about 10 minutes (mostly the original algorithm).
# =============================================================================
suppressMessages(library(dplyr))
stopifnot(requireNamespace("scoth", quietly = TRUE), requireNamespace("scoth.fast", quietly = TRUE))
ns = asNamespace("scoth.fast")
RESULTS = file.path(here::here(), "results")
dir.create(RESULTS, showWarnings = FALSE)
N_REP = 50
if (nzchar(Sys.getenv("SCOTH_QUICK"))) {            # quick check of the script
  N_REP = 2; RESULTS = file.path(tempdir(), "results_quick"); dir.create(RESULTS, showWarnings = FALSE)
}
TYPES = c("continuous", "integer", "decimal", "const_runs")
LENGTHS = c(100, 400)
COSTS = c("cusum_L2", "reg_cusum_L2")
L0_GRID = c(0, 64, Inf)
thresholds = c(seq(from = 0.03, to = 0.2, by = 0.01), 0.3, 0.4, 0.5, 0.6, 0.7, 1)

make_series = function(type, T, seed) {
  set.seed(seed)
  cps = sort(sample(10:(T - 10), 3))
  len = diff(c(1, cps, T + 1))
  switch(type,
    continuous = rep(rnorm(4, sd = 2), len) + rnorm(T),
    integer    = as.numeric(rpois(T, rep(c(3, 8, 4, 10)[sample(4)], len))),
    decimal    = round(cumsum(rnorm(T, sd = 0.3)) + rep(c(0, 3, -2, 1), len), 1),
    const_runs = rep(round(rnorm(4, sd = 2), 2), len) + ifelse(runif(T) < 0.5, 0, round(rnorm(T), 2)))
}

# selections of Algorithm 2 for all thresholds and levels (exact fast selection path)
selections = function(score, cost_initial, x, cost) {
  options(scoth.fast_selection = TRUE)
  cf = ns$all_cost_functions()[[cost]]
  lapply(thresholds, function(thr) {
    if (thr > max(score)) return(NULL)   # the pipeline also skips these thresholds (no change point)
    lapply(ns$scoth_recursive_partitions(thr, 10, score, cost_initial, x, cf)$partitions, function(d) d$s)
  })
}

rows = list()
for (type in TYPES) for (T in LENGTHS) for (cost in COSTS) {
  for (rep in seq_len(N_REP)) {
    x = make_series(type, T, seed = rep)
    orig = scoth::scoth_scoring_func(x, cost, verbose = FALSE)
    sel_orig = selections(orig$changepoint_gain, orig$cost_initial, x, cost)
    for (L0 in L0_GRID) {
      fast = ns$scoth_changes_fast(x, cost, exact_max_len = L0)
      same_scores = identical(fast, orig$changepoint_gain)
      n_same = if (same_scores) {
        length(thresholds) * 10           # identical scores -> identical selections (deterministic)
      } else {
        sel_fast = selections(fast, orig$cost_initial, x, cost)
        sum(unlist(mapply(function(a, b) {
          if (is.null(a) && is.null(b)) return(rep(TRUE, 10))
          if (is.null(a) || is.null(b)) return(rep(FALSE, 10))
          mapply(identical, a, b)
        }, sel_orig, sel_fast, SIMPLIFY = FALSE)))
      }
      rows[[length(rows) + 1]] = data.frame(
        data_type = type, T = T, cost = cost, seed = rep, exact_max_len = L0,
        bit_identical = same_scores,
        max_abs_score_diff = max(abs(fast - orig$changepoint_gain)),
        n_selection_settings = length(thresholds) * 10,
        n_identical_selections = n_same)
    }
  }
  message(sprintf("%-10s T=%d %-12s done", type, T, cost))
}
rows = bind_rows(rows)
write.csv(rows, file.path(RESULTS, "exactness_random_series.csv"), row.names = FALSE)

summ = rows %>% group_by(exact_max_len, data_type, T, cost) %>%
  summarize(n_series = n(),
            n_bit_identical = sum(bit_identical),
            n_max_diff_le_1e12 = sum(max_abs_score_diff <= 1e-12),
            n_all_selections_identical = sum(n_identical_selections == n_selection_settings),
            pct_identical_selection_settings = round(100 * sum(n_identical_selections) / sum(n_selection_settings), 2),
            max_abs_score_diff = max(max_abs_score_diff), .groups = "drop")
write.csv(summ, file.path(RESULTS, "exactness_random_summary.csv"), row.names = FALSE)
options(width = 200)
print(as.data.frame(summ))
overall = rows %>% group_by(exact_max_len) %>%
  summarize(n_series = n(), n_bit_identical = sum(bit_identical),
            n_all_selections_identical = sum(n_identical_selections == n_selection_settings),
            max_abs_score_diff = max(max_abs_score_diff), .groups = "drop")
print(as.data.frame(overall))
