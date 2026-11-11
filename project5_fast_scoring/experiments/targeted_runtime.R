# =============================================================================
# targeted_runtime.R   [IMPROVEMENT-5]  targeted experiment: run time of the scoring
#
# Question: how does the time of Algorithm 1 (scoring) grow with the length T of
# the series, for the original implementation (package scoth) and for the fast
# exact implementation (package scoth.fast)?
#
# Methods (same input series, same cost function):
#   original   : scoth::scoth_scoring_func                     (reference code, baseline)
#   fast       : scoth.fast::scoth_scoring_func(algorithm="fast"), exact_max_len = 64
#                (replica costs for segments <= 64 points, O(1) prefix sums otherwise)
#   fast_exact : same with exact_max_len = Inf (replica costs for every segment,
#                bit-identical to the original)
# Series: piecewise-constant mean (5 random change points, levels ~ N(0, 2^2)) + N(0,1) noise,
#         seeds 1, 2, 3 (one series per repetition), T in {250, ..., 100000}.
# The original is run only while one run is predicted to take <= 300 s (log-log slope
# of the last two measured lengths); larger T are EXTRAPOLATED with the slope fitted
# on its largest measured lengths (flag `extrapolated` in the CSV files).
# The cap depends on the speed of the machine: on a machine about 1.5-2 times slower
# than ours the original is already skipped at T = 32000. The environment variable
# SCOTH_MAX_SECONDS_ORIGINAL changes the cap (default 300; e.g. 600
# on a slower machine). Do not set it close to the times predicted for
# T = 100000, one per cost function (log lines "[<cost>] original skipped from
# T = 100000 on (predicted ... s ...)": 842-1461 s in our 5 runs). With a cap at or
# above the prediction of one cost function, the original also runs at T = 100000 for
# that cost, up to 3 times (one per series) of about 14-26 min each; with a cap at or
# above both predictions (the larger one: 1043-1461 s in our 5 runs), up to 6 times.
# For every run of the original, the fast scores are compared with the original ones
# (max |difference|, bit-identity) and the selections of Algorithm 2 are compared.
# Worst case of fast_exact (replica costs for every segment): two long exactly
# constant runs, where one segment grows by one point at each iteration.
# Last part: total time of the scoring stage on the 42 TCPD datasets (both costs).
#
# Usage (from the project folder):  Rscript experiments/targeted_runtime.R
# Outputs (results/): runtime_raw.csv, runtime_summary.csv, runtime_slopes.csv,
#                     runtime_original_extrapolation.csv, runtime_loglog.png,
#                     runtime_worst_case.csv, tcpd_scoring_time.csv
# Another output folder (relative to the project folder) can be given in the
# environment variable SCOTH_RESULTS_DIR, e.g. to keep a second run separate and
# compare the timings of the two runs.
# Run time: about 30-40 minutes (almost all of it is the original algorithm).
# =============================================================================
suppressMessages({ library(dplyr); library(ggplot2) })
stopifnot(requireNamespace("scoth", quietly = TRUE), requireNamespace("scoth.fast", quietly = TRUE))

OUT_DIR = gsub("\\\\", "/", Sys.getenv("SCOTH_RESULTS_DIR"))   # e.g. "results\run_2" in PowerShell
if (!nzchar(OUT_DIR)) OUT_DIR = "results"
RESULTS = file.path(here::here(), OUT_DIR)
T_GRID = c(250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 100000)
REPS = 3
MAX_SECONDS_ORIGINAL = as.numeric(Sys.getenv("SCOTH_MAX_SECONDS_ORIGINAL", "300"))  # cap of one run of the original
stopifnot(is.finite(MAX_SECONDS_ORIGINAL), MAX_SECONDS_ORIGINAL > 0)
COSTS = c("cusum_L2", "reg_cusum_L2")
METHODS = c("original", "fast", "fast_exact")
SEL_THRESHOLDS = c(0.03, 0.05, 0.1, 0.2)
TCPD_DATASETS = Inf                     # number of TCPD datasets used in the last part (all)
# quick check of the script (about 1-3 minutes): set the environment variable SCOTH_QUICK=1
if (nzchar(Sys.getenv("SCOTH_QUICK"))) {
  T_GRID = c(250, 500, 1000, 2000); REPS = 2; TCPD_DATASETS = 3
  RESULTS = file.path(tempdir(), "results_quick")
  message("QUICK mode")
}
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
message("outputs in ", RESULTS)
message("the original is skipped when one run is predicted to take more than ", MAX_SECONDS_ORIGINAL, " s")

make_series = function(T, seed) {
  set.seed(seed)
  K = 5
  cps = sort(sample(2:T, K))
  levels = rnorm(K + 1, sd = 2)
  rep(levels, diff(c(1, cps, T + 1))) + rnorm(T)
}

run_method = function(method, x, cost) {
  f = if (method == "original") {
    function() scoth::scoth_scoring_func(x, cost, verbose = FALSE)
  } else {
    options(scoth.fast_exact_max_len = if (method == "fast") 64 else Inf)
    function() scoth.fast::scoth_scoring_func(x, cost, verbose = FALSE, algorithm = "fast")
  }
  gc(verbose = FALSE)
  t = system.time(out <- f())[["elapsed"]]
  if (t < 0.2) {             # too short for the ~10 ms clock resolution: repeat and average
    n = ceiling(0.2 / max(t, 0.005))
    t = system.time(for (i in seq_len(n)) out <- f())[["elapsed"]] / n
  }
  list(seconds = t, score = out$changepoint_gain, cost_initial = out$cost_initial)
}

# selections of Algorithm 2 (levels 1..10 at a few thresholds), computed with the
# exact fast selection path; used only to compare two score vectors
selections = function(score, cost_initial, x, cost) {
  options(scoth.fast_selection = TRUE)
  cf = scoth.fast::all_cost_functions()[[cost]]
  lapply(SEL_THRESHOLDS[SEL_THRESHOLDS <= max(score)], function(thr) {
    p = scoth.fast::scoth_recursive_partitions(thr, 10, score, cost_initial, x, cf)$partitions
    lapply(p, function(d) d$s)
  })
}

raw = list()
for (cost in COSTS) {
  original_allowed = TRUE
  measured = data.frame(T = numeric(0), sec = numeric(0))
  for (T in T_GRID) {
    # ---- may the original run at this length? ----
    if (original_allowed && nrow(measured) >= 2) {
      n = nrow(measured)
      slope = diff(log(measured$sec[(n-1):n])) / diff(log(measured$T[(n-1):n]))
      predicted = measured$sec[n] * (T / measured$T[n])^max(slope, 1)
      if (predicted > MAX_SECONDS_ORIGINAL) {
        message(sprintf("[%s] original skipped from T = %d on (predicted %.0f s > %g s)", cost, T, predicted, MAX_SECONDS_ORIGINAL))
        original_allowed = FALSE
      }
    }
    t_orig = c()
    for (rep in seq_len(REPS)) {
      x = make_series(T, seed = rep)
      res = list()
      for (method in METHODS) {
        if (method == "original" && !original_allowed) next
        res[[method]] = run_method(method, x, cost)
      }
      for (method in names(res)) {
        r = res[[method]]
        row = data.frame(cost = cost, T = T, rep = rep, method = method, seconds = r$seconds,
                         max_abs_diff_vs_original = NA_real_, identical_to_original = NA,
                         same_selections_as_original = NA)
        if (method != "original" && !is.null(res$original)) {
          o = res$original
          row$max_abs_diff_vs_original = max(abs(r$score - o$score))
          row$identical_to_original = identical(r$score, o$score)
          row$same_selections_as_original = identical(selections(r$score, o$cost_initial, x, cost),
                                                      selections(o$score, o$cost_initial, x, cost))
        }
        raw[[length(raw) + 1]] = row
      }
      if (!is.null(res$original)) t_orig = c(t_orig, res$original$seconds)
      message(sprintf("[%s] T = %6d rep %d: %s", cost, T, rep,
                      paste(sprintf("%s %.2f s", names(res), sapply(res, `[[`, "seconds")), collapse = ", ")))
      # a single run of the original far above the limit: do not repeat it
      if (length(t_orig) > 0 && tail(t_orig, 1) > 1.5 * MAX_SECONDS_ORIGINAL) original_allowed = FALSE
    }
    if (length(t_orig) > 0) measured[nrow(measured) + 1, ] = c(T, median(t_orig))
  }
}
raw = bind_rows(raw)
write.csv(raw, file.path(RESULTS, "runtime_raw.csv"), row.names = FALSE)

# ---- medians, slopes, extrapolation ------------------------------------------
summ = raw %>% group_by(cost, method, T) %>%
  summarize(n_reps = n(), median_seconds = median(seconds), min_seconds = min(seconds), max_seconds = max(seconds),
            max_abs_diff_vs_original = suppressWarnings(max(max_abs_diff_vs_original, na.rm = TRUE)),
            all_identical_to_original = all(identical_to_original), all_same_selections = all(same_selections_as_original),
            .groups = "drop") %>%
  mutate(max_abs_diff_vs_original = ifelse(is.finite(max_abs_diff_vs_original), max_abs_diff_vs_original, NA),
         extrapolated = FALSE)

fit_slope = function(d) {
  m = lm(log10(median_seconds) ~ log10(T), data = d)
  data.frame(slope = unname(coef(m)[2]), intercept = unname(coef(m)[1]), r2 = summary(m)$r.squared,
             T_min = min(d$T), T_max = max(d$T), n_points = nrow(d))
}
# slopes on the large lengths (T >= 2000), where the asymptotic behaviour shows
slopes = summ %>% filter(T >= 2000) %>% group_by(cost, method) %>% group_modify(~ fit_slope(.x)) %>% ungroup()
# for the extrapolation of the original: slope of its 3 largest measured lengths
slopes_extra = summ %>% filter(method == "original") %>% group_by(cost) %>% slice_max(T, n = 3) %>%
  group_modify(~ fit_slope(.x)) %>% ungroup() %>% mutate(method = "original (3 largest T, used for extrapolation)")
slopes = bind_rows(slopes, slopes_extra)
write.csv(slopes, file.path(RESULTS, "runtime_slopes.csv"), row.names = FALSE)

extra = list()
for (cs in COSTS) {
  done = summ %>% filter(cost == cs, method == "original") %>% pull(T)
  sl = slopes_extra %>% filter(cost == cs)
  for (T in setdiff(T_GRID, done)) {
    extra[[length(extra) + 1]] = data.frame(cost = cs, method = "original", T = T, n_reps = 0L,
      median_seconds = 10^(sl$intercept + sl$slope * log10(T)), min_seconds = NA, max_seconds = NA,
      max_abs_diff_vs_original = NA, all_identical_to_original = NA, all_same_selections = NA, extrapolated = TRUE)
  }
}
summ = bind_rows(summ, bind_rows(extra)) %>% arrange(cost, method, T)
speed = summ %>% select(cost, method, T, median_seconds, extrapolated) %>%
  tidyr::pivot_wider(names_from = method, values_from = c(median_seconds, extrapolated)) %>%
  transmute(cost, T, original_seconds = median_seconds_original, original_extrapolated = extrapolated_original,
            fast_seconds = median_seconds_fast, fast_exact_seconds = median_seconds_fast_exact,
            speedup_fast = original_seconds / fast_seconds, speedup_fast_exact = original_seconds / fast_exact_seconds)
summ = summ %>% left_join(speed %>% select(cost, T, speedup_fast, speedup_fast_exact), by = c("cost", "T"))
write.csv(summ, file.path(RESULTS, "runtime_summary.csv"), row.names = FALSE)
print(as.data.frame(speed), digits = 4)
print(as.data.frame(slopes), digits = 3)

# ---- a second extrapolation of the original: time = a*T + b*T^2 -------------------
# The bookkeeping of the original is O(T^2) but each lm() call (linear cost) adds a
# term proportional to T; at T <= 32000 this linear term still lowers the log-log
# slope of the linear cost. This two-term least-squares fit (on T >= 2000) is a
# second, less conservative extrapolation; the figure and the summary use the
# log-log one, as planned.
two_term = bind_rows(lapply(COSTS, function(cs) {
  d = summ %>% filter(cost == cs, method == "original", !extrapolated, T >= 2000)
  loglog = summ %>% filter(cost == cs, method == "original", T == max(T_GRID))
  if (nrow(d) < 2) {   # only with a tiny SCOTH_MAX_SECONDS_ORIGINAL: too few lengths for this fit
    return(data.frame(cost = cs, T = max(T_GRID), a_seconds_per_point = NA_real_, b_seconds_per_point2 = NA_real_,
                      two_term_prediction_seconds = NA_real_,
                      loglog_prediction_seconds = loglog$median_seconds, loglog_is_extrapolated = loglog$extrapolated))
  }
  m = lm(median_seconds ~ 0 + T + I(T^2), data = d)
  data.frame(cost = cs, T = max(T_GRID), a_seconds_per_point = unname(coef(m)[1]), b_seconds_per_point2 = unname(coef(m)[2]),
             two_term_prediction_seconds = unname(predict(m, data.frame(T = max(T_GRID)))),
             loglog_prediction_seconds = loglog$median_seconds, loglog_is_extrapolated = loglog$extrapolated)
}))
write.csv(two_term, file.path(RESULTS, "runtime_original_extrapolation.csv"), row.names = FALSE)
print(two_term)

# ---- figure: log-log plot with fitted slopes ----------------------------------
METHOD_COLOURS = c(original = "#2a78d6", fast = "#eb6834", fast_exact = "#1baf7a")  # validated categorical slots 1-3
METHOD_LABELS = c(original = "original (scoth)", fast = "fast (exact_max_len = 64)", fast_exact = "fast_exact (exact_max_len = Inf)")
COST_LABELS = c(cusum_L2 = "quadratic cost (cusum_L2)", reg_cusum_L2 = "linear cost (reg_cusum_L2)")
fits = slopes %>% filter(!grepl("3 largest", method))
plot_df = summ %>% mutate(status = ifelse(extrapolated, "extrapolated (original only)", "measured (median of 3 runs)"))
lines_df = bind_rows(lapply(seq_len(nrow(fits)), function(i) {
  f = fits[i, ]
  data.frame(cost = f$cost, method = f$method, T = c(f$T_min, f$T_max),
             seconds = 10^(f$intercept + f$slope * log10(c(f$T_min, f$T_max))))
}))
extra_df = summ %>% filter(method == "original") %>% group_by(cost) %>%
  filter(T >= max(T[!extrapolated])) %>% ungroup()
# direct labels (text in grey ink, not in the series colour): slopes at the end of the
# fitted lines; fast and fast_exact overlap, so they share one label above their lines
sl3 = slopes %>% filter(grepl("3 largest", method)) %>% select(cost, slope3 = slope)
lab_orig = fits %>% filter(method == "original") %>% left_join(sl3, by = "cost") %>%
  transmute(cost, T = T_max, seconds = 2.6 * 10^(intercept + slope * log10(T_max)),
            label = sprintf("original: slope %.2f\n(%.2f on the 3 largest T)", slope, slope3))
lab_fast = fits %>% filter(method != "original") %>% group_by(cost) %>%
  summarize(T = max(T_max), seconds = 5 * 10^max(intercept + slope * log10(T_max)),
            label = paste(sprintf("%s: slope %.2f", method, slope), collapse = "\n"), .groups = "drop")
lab_df = bind_rows(lab_orig, lab_fast)
g = ggplot(plot_df, aes(T, median_seconds, colour = method)) +
  geom_line(data = lines_df, aes(T, seconds, colour = method), linewidth = 0.7) +
  geom_line(data = extra_df, aes(T, median_seconds, colour = method), linewidth = 0.7, linetype = "22") +
  geom_point(aes(shape = status), size = 2.6, stroke = 0.9) +
  geom_text(data = lab_df, aes(T, seconds, label = label), colour = "#52514e", size = 3.1, hjust = 1, inherit.aes = FALSE) +
  scale_colour_manual(values = METHOD_COLOURS, labels = METHOD_LABELS, breaks = names(METHOD_LABELS)) +
  scale_shape_manual(values = c("extrapolated (original only)" = 1, "measured (median of 3 runs)" = 16)) +
  scale_x_log10(breaks = T_GRID, labels = format(T_GRID, big.mark = ",", scientific = FALSE, trim = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), labels = c("0.01", "0.1", "1", "10", "100", "1,000")) +
  annotation_logticks(sides = "l", colour = "grey60") +
  facet_wrap(~ cost, labeller = labeller(cost = COST_LABELS)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "series length T (log scale)", y = "scoring time in seconds (log scale)", colour = NULL, shape = NULL,
       title = "Algorithm 1 (scoring): original implementation vs fast exact implementation",
       subtitle = "Lines: least-squares fit of log10(time) on log10(T) for T >= 2000 (slope 1 = linear, 2 = quadratic). Dashed: extrapolation of the original.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle = element_text(size = 8.5, colour = "#52514e"), strip.background = element_rect(fill = "grey95"))
ggsave(file.path(RESULTS, "runtime_loglog.png"), g, width = 10, height = 5.8, dpi = 150)


# ---- worst case of fast_exact: two long, exactly constant runs ----------------------
# x = 0,...,0, 1,...,1. Inside each run every gain is exactly 0, the ties are broken
# by the smallest index, so the cuts are removed from left to right and ONE segment
# grows by one point at each iteration. With replica costs for every segment
# (fast_exact) each new segment cost costs O(n) => O(T^2) in total, while the default
# fast version computes the long segments in O(1) (prefix sums + constant-run test).
WORST_T = if (nzchar(Sys.getenv("SCOTH_QUICK"))) c(1000, 2000) else c(2000, 8000, 32000)
worst = list()
for (cost in COSTS) for (T in WORST_T) {
  x = rep(c(0, 1), each = T / 2)
  res = lapply(c("fast", "fast_exact"), function(method) {
    r = lapply(seq_len(REPS), function(i) run_method(method, x, cost))
    list(method = method, seconds = median(sapply(r, `[[`, "seconds")), score = r[[1]]$score)
  })
  for (r in res) {
    worst[[length(worst) + 1]] = data.frame(cost = cost, T = T, method = r$method, median_seconds = r$seconds,
                                            identical_fast_vs_fast_exact = identical(res[[1]]$score, res[[2]]$score),
                                            max_abs_diff_fast_vs_fast_exact = max(abs(res[[1]]$score - res[[2]]$score)))
  }
  message(sprintf("[worst case %s] T = %d: fast %.2f s, fast_exact %.2f s", cost, T, res[[1]]$seconds, res[[2]]$seconds))
}
worst = bind_rows(worst) %>% group_by(cost, method) %>%
  mutate(loglog_slope = unname(coef(lm(log10(median_seconds) ~ log10(T)))[2])) %>% ungroup()
write.csv(worst, file.path(RESULTS, "runtime_worst_case.csv"), row.names = FALSE)
print(as.data.frame(worst))

# ---- total scoring-stage time on the 42 TCPD datasets --------------------------
kpis = scoth.fast::load_kpis(file.path(here::here(), "public-data"))
tcpd = list()
for (rep in seq_len(REPS)) {
  for (method in METHODS) {
    sec = 0                               # sum of the scoring times (gc() calls not included)
    for (i in seq_len(min(nrow(kpis), TCPD_DATASETS))) {
      k = kpis$kpis[[i]]
      x = if (ncol(k) == 1) k[, 1] else k
      for (cost in c("reg_cusum_L2", "cusum_L2")) sec = sec + run_method(method, x, cost)$seconds
    }
    tcpd[[length(tcpd) + 1]] = data.frame(rep = rep, method = method, seconds = sec)
    message(sprintf("TCPD scoring stage, rep %d, %s: %.1f s", rep, method, sec))
  }
}
tcpd = bind_rows(tcpd) %>% group_by(method) %>%
  summarize(median_seconds = median(seconds), min_seconds = min(seconds), max_seconds = max(seconds), .groups = "drop") %>%
  mutate(speedup_vs_original = median_seconds[method == "original"] / median_seconds)
write.csv(tcpd, file.path(RESULTS, "tcpd_scoring_time.csv"), row.names = FALSE)
print(as.data.frame(tcpd))
