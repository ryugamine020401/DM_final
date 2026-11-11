# [IMPROVEMENT-5] Tests of the fast exact scoring algorithm (R/scoth_fast_scoring.R)
# and of the fast exact selection path (options(scoth.fast_selection = TRUE)).
#
# The reference is the original implementation of this package (scoth_changes,
# build_df_df2, ...), which is unchanged, so the tests do not need the baseline package.

test_series_univariate = function() {
  set.seed(5)
  list(
    random     = c(rnorm(60), rnorm(60, 2)),                               # no ties
    integers   = as.numeric(sample(0:3, 120, replace = TRUE)),             # many exactly tied gains
    const_runs = rep(c(1, 1.5, 4, 2, 0.1, 2, 7, 0.3),                     # constant runs (exact 0 gains)
                     times = c(10, 5, 3, 12, 1, 8, 6, 15)),
    decimals   = round(cumsum(rnorm(150, sd = 0.3)), 1),                   # 1 decimal: near ties
    pattern    = rep(c(0, 0, 1, 1, 5), 25),                                # repeated pattern
    trend      = c(1:40, 40 - 0.5 * (1:40), rep(20, 30)) + 0               # collinear runs (linear cost)
  )
}

test_series_multivariate = function() {
  set.seed(7)
  list(
    random  = data.frame(a = c(rnorm(50), rnorm(50, 3)), b = rnorm(100, 0, 5), c = c(rnorm(70), rnorm(30, -2))),
    integer = data.frame(a = as.numeric(sample(0:2, 90, replace = TRUE)), b = as.numeric(rep(c(1, 1, 2), 30))),
    # an integer column, like "Volume" in the apple dataset (mean() of integers is computed differently)
    mixed   = data.frame(price = round(c(rnorm(60, 10), rnorm(60, 12)), 2),
                         volume = as.integer(round(c(rlnorm(60, 15), rlnorm(60, 16)))))
  )
}

cost_fun_of = function(name) all_cost_functions()[[name]]

test_that("[IMPROVEMENT-5] replicas of the cost functions are bit-identical to the original ones", {
  set.seed(11)
  for (x in test_series_univariate()) {
    T = length(x)
    for (i in 1:60) {
      len = sample(1:min(40, T), 1); s = sample(1:(T - len + 1), 1); e = s + len
      expect_identical(replica_c_cusum_L2(x, s, e), c_cusum_L2(x, s, e))
      expect_identical(replica_c_reg_cusum_L2(x, s, e), c_reg_cusum_L2(x, s, e))
    }
  }
  for (x in test_series_multivariate()) {
    M = as.matrix(x); storage.mode(M) = "double"; T = nrow(x)
    cols = lapply(seq_len(ncol(x)), function(f) x[[f]])
    for (i in 1:40) {
      len = sample(1:min(40, T), 1); s = sample(1:(T - len + 1), 1); e = s + len
      expect_identical(replica_c_cusum_L2_multivariate(cols, s, e), c_cusum_L2_multivariate(x, s, e))
      expect_identical(replica_c_reg_cusum_L2_multivariate(M, s, e), c_reg_cusum_L2_multivariate(x, s, e))
    }
  }
})

test_that("[IMPROVEMENT-5] fast scoring with exact costs (exact_max_len = Inf) is bit-identical to scoth_changes", {
  for (x in test_series_univariate()) {
    for (cn in c("cusum_L2", "reg_cusum_L2")) {
      orig = scoth_changes(x, cost_fun_of(cn), s = 1, e = length(x) + 1, verbose = FALSE)
      fast = scoth_changes_fast(x, cn, exact_max_len = Inf)
      expect_identical(fast, orig)
    }
  }
  for (x in test_series_multivariate()) {
    for (cn in c("cusum_L2_multivariate", "reg_cusum_L2_multivariate")) {
      orig = scoth_changes(x, cost_fun_of(cn), s = 1, e = nrow(x) + 1, verbose = FALSE)
      expect_identical(scoth_changes_fast(x, cn, exact_max_len = Inf), orig)
    }
  }
})

same_selections = function(x, cn2, orig, fast, thresholds = c(0.03, 0.05, 0.1, 0.2)) {
  T = if (is.data.frame(x)) nrow(x) else length(x)
  cost_initial = cost_fun_of(cn2)(x, 1, T + 1)
  for (thr in thresholds[thresholds <= max(orig)]) {
    p_orig = scoth_recursive_partitions(thr, 10, orig, cost_initial, x, cost_fun_of(cn2))$partitions
    p_fast = scoth_recursive_partitions(thr, 10, fast, cost_initial, x, cost_fun_of(cn2))$partitions
    if (!identical(lapply(p_fast, function(p) p$s), lapply(p_orig, function(p) p$s))) return(FALSE)
  }
  TRUE
}

test_that("[IMPROVEMENT-5] fast scoring with prefix sums for long segments gives the same scores and selections", {
  series = c(test_series_univariate(), test_series_multivariate())
  for (nm in names(series)) {
    x = series[[nm]]
    multi = is.data.frame(x)
    T = if (multi) nrow(x) else length(x)
    for (cn in c("cusum_L2", "reg_cusum_L2")) {
      cn2 = if (multi) paste0(cn, "_multivariate") else cn
      orig = scoth_changes(x, cost_fun_of(cn2), s = 1, e = T + 1, verbose = FALSE)
      # default exact_max_len = 64: tie-heavy short segments are computed exactly
      fast = scoth_changes_fast(x, cn2)
      expect_lt(max(abs(fast - orig)), 1e-12)
      expect_true(same_selections(x, cn2, orig, fast))
      # exact_max_len = 20: many more prefix-sum costs. Not used for the noise-free
      # piecewise-linear "trend" series, see the next test.
      if (nm != "trend") {
        fast20 = scoth_changes_fast(x, cn2, exact_max_len = 20)
        expect_lt(max(abs(fast20 - orig)), 1e-12)
        expect_true(same_selections(x, cn2, orig, fast20))
      }
    }
  }
})

test_that("[IMPROVEMENT-5] degenerate noise-free series: exact copies are needed to follow the original", {
  # x is EXACTLY piecewise linear: with the linear cost the level-2 segmentation
  # has V(S) ~ 1e-30 (rounding noise of lm), so the zoom factor V(empty)/V(S) ~ 1e30
  # turns scores that are mathematically 0 (computed as 1e-32 by the original) into
  # O(1) numbers. The original selection is then driven by rounding noise; it is
  # reproduced only when these segments are computed with the same arithmetic.
  x = test_series_univariate()$trend
  orig = scoth_changes(x, c_reg_cusum_L2, s = 1, e = length(x) + 1, verbose = FALSE)
  fast_inf = scoth_changes_fast(x, "reg_cusum_L2", exact_max_len = Inf)
  expect_identical(fast_inf, orig)
  expect_true(same_selections(x, "reg_cusum_L2", orig, fast_inf))
  fast20 = scoth_changes_fast(x, "reg_cusum_L2", exact_max_len = 20)
  expect_lt(max(abs(fast20 - orig)), 1e-12)                       # same scores up to rounding ...
  expect_false(same_selections(x, "reg_cusum_L2", orig, fast20))  # ... but noise is amplified 1e30 times
})

test_that("[IMPROVEMENT-5] ties are broken like which.min (smallest index first)", {
  # all initial gains are exactly equal: the original removes the cuts from left to right
  x = rep(c(0, 1), 20)
  orig = scoth_changes(x, c_cusum_L2, s = 1, e = length(x) + 1, verbose = FALSE)
  expect_identical(scoth_changes_fast(x, "cusum_L2"), orig)
  # constant runs: exact zero gains inside the runs, for both costs
  x = rep(c(3, 3.7, 0.1), each = 30)
  for (cn in c("cusum_L2", "reg_cusum_L2")) {
    orig = scoth_changes(x, cost_fun_of(cn), s = 1, e = length(x) + 1, verbose = FALSE)
    expect_identical(scoth_changes_fast(x, cn, exact_max_len = Inf), orig)
    expect_lt(max(abs(scoth_changes_fast(x, cn, exact_max_len = 10) - orig)), 1e-12)
  }
})

test_that("[IMPROVEMENT-5] scoth_scoring_func(algorithm = 'fast') returns the same list as the original", {
  set.seed(3)
  x = c(rnorm(80), rnorm(80, 1.5))
  for (cn in c("cusum_L2", "reg_cusum_L2")) {
    a = scoth_scoring_func(x, cn, verbose = FALSE, algorithm = "original")
    b = scoth_scoring_func(x, cn, verbose = FALSE, algorithm = "fast")
    expect_identical(b$cost_initial, a$cost_initial)
    expect_lt(max(abs(b$changepoint_gain - a$changepoint_gain)), 1e-12)
  }
  # the global option selects the algorithm (default = original)
  old = options(scoth.scoring_algorithm = "fast")
  on.exit(options(old))
  expect_lt(max(abs(scoth_scoring_func(x, "cusum_L2", verbose = FALSE)$changepoint_gain -
                      scoth_changes(x, c_cusum_L2, 1, length(x) + 1, verbose = FALSE))), 1e-12)
  # unknown cost function (c_box, used by the original tests): the original algorithm runs
  y = c(1, 1, 10, 1)
  expect_equal(scoth_scoring_func(y, c_box, verbose = FALSE, algorithm = "fast")$changepoint_gain, c(0, 0.75, 0.5))
})

test_that("[IMPROVEMENT-5] corner cases", {
  expect_identical(scoth_changes_fast(rep(2.5, 10), "cusum_L2"), rep(0, 9))          # constant series
  # (with the linear cost, the original V(empty) of a constant series is lm rounding
  #  noise, not 0, so its scores are noise/noise; the fast version reproduces them)
  expect_identical(scoth_changes_fast(rep(2.5, 10), "reg_cusum_L2"),
                   scoth_changes(rep(2.5, 10), c_reg_cusum_L2, 1, 11, verbose = FALSE))
  expect_identical(scoth_changes_fast(c(1, 3), "cusum_L2"),
                   scoth_changes(c(1, 3), c_cusum_L2, 1, 3, verbose = FALSE))         # T = 2
  expect_identical(scoth_changes_fast(1:30, "cusum_L2", exact_max_len = Inf),       # integer input
                   scoth_changes(1:30, c_cusum_L2, 1, 31, verbose = FALSE))
  expect_false(fast_input_ok(c(1, NA, 3)))                                          # NA -> original algorithm
  expect_false(fast_input_ok(data.frame(a = c(1, 2), b = c("u", "v"))))             # text column -> original
  expect_true(fast_input_ok(data.frame(a = 1:5, b = c(0.5, 1, 2, 3, 1))))           # integer column is fine
  expect_true(fast_input_ok(data.frame(a = c(1, 2), b = c(3, 4))))
  expect_identical(fast_cost_name(c_reg_cusum_L2), "reg_cusum_L2")
  expect_null(fast_cost_name(c_box))
})

test_that("[IMPROVEMENT-5] O(1) prefix-sum costs agree with the exact costs on long series", {
  set.seed(9)
  T = 3000
  x = 50 + 0.01 * (1:T) + rep(c(0, 4, -2, 6), each = T / 4) + rnorm(T)  # trend + level shifts
  M = data.frame(a = x, b = rnorm(T, 100, 3))
  for (cn in c("cusum_L2", "reg_cusum_L2", "cusum_L2_multivariate", "reg_cusum_L2_multivariate")) {
    xx = if (grepl("multivariate", cn)) M else x
    oracle = fast_segment_cost(xx, cn, exact_max_len = 0)
    for (i in 1:40) {
      len = sample(3:T, 1); s = sample(1:(T - len + 1), 1); e = s + len
      a = oracle$cost_prefix(s, e); b = oracle$cost_exact(s, e)
      expect_lt(abs(a - b) / b, 1e-9)
    }
  }
})

test_that("[IMPROVEMENT-5] fast selection path gives identical partitions", {
  set.seed(21)
  x = c(rnorm(70), rnorm(50, 3), rnorm(80, -1))
  xm = data.frame(a = x, b = c(rnorm(100), rnorm(100, 2)))
  old = options(scoth.fast_selection = FALSE)
  on.exit(options(old))
  for (xx in list(x, xm)) {
    multi = is.data.frame(xx)
    T = if (multi) nrow(xx) else length(xx)
    for (cn in c("cusum_L2", "reg_cusum_L2")) {
      cn2 = if (multi) paste0(cn, "_multivariate") else cn
      cf = cost_fun_of(cn2)
      sc = scoth_scoring_func(xx, cn, verbose = FALSE)
      for (thr in c(0.03, 0.1)) {
        options(scoth.fast_selection = FALSE)
        a = scoth_recursive_partitions(thr, 10, sc$changepoint_gain, sc$cost_initial, xx, cf)
        options(scoth.fast_selection = TRUE)
        b = scoth_recursive_partitions(thr, 10, sc$changepoint_gain, sc$cost_initial, xx, cf)
        expect_identical(b, a)
      }
      # partition with a zero-cost segment and length-1 segments
      idx = c(5, 6, 30, 31, 32, 100)
      options(scoth.fast_selection = FALSE); pa = get_partition_from_idx(idx, xx, cf, 1, T + 1)
      options(scoth.fast_selection = TRUE);  pb = get_partition_from_idx(idx, xx, cf, 1, T + 1)
      expect_identical(pb, pa)
      pa$cost[2] = 0
      options(scoth.fast_selection = FALSE); ga = iterate_changepoint_gain_from_partition(pa, sc$changepoint_gain, sc$cost_initial)
      options(scoth.fast_selection = TRUE);  gb = iterate_changepoint_gain_from_partition(pa, sc$changepoint_gain, sc$cost_initial)
      expect_identical(gb, ga)
    }
  }
})
