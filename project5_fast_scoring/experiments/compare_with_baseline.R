# =============================================================================
# compare_with_baseline.R   [IMPROVEMENT-5]
#
# Puts the BASELINE (original scoth code, outputs of project0_reproduction, read
# only) and the IMPROVED pipeline (experiments/main_fast.R, outputs of this
# project) side by side on the 42 TCPD datasets:
#   1. scores of Algorithm 1 (max |difference|, bit-identity) and scoring times
#   2. selections of Algorithm 2 (all thresholds x levels) and their partitions
#   3. evaluation tables (F1, F1-biased, cover) row by row
#   4. the three protocols: ORACLE, BEST-SINGLE, FULLY-AUTOMATIC (in-sample + LOO-CV)
#      + per-dataset F1 differences + paired Wilcoxon signed-rank tests
#   5. wall-clock time of each stage of the pipeline (parsed from the two logs)
#
# Usage (from the project folder, after main_fast.R and extra_protocols.R):
#   Rscript experiments/compare_with_baseline.R
# Outputs (results/): tcpd_score_differences.csv, tcpd_selection_identity.csv,
#   tcpd_metric_differences.csv, tcpd_protocols_comparison.csv,
#   tcpd_per_dataset_F1_differences.csv, tcpd_wilcoxon.csv, tcpd_stage_timings.csv,
#   comparison_summary.txt                                  Run time: about 1 minute.
# =============================================================================
suppressMessages({ library(dplyr); library(tidyr) })
options(width = 200)    # print the tables of the report on one line

ROOT = here::here()
BASE_ROOT = file.path(dirname(ROOT), "project0_reproduction")       # baseline (read only)
BASE = file.path(BASE_ROOT, "output-data")
IMPR = file.path(ROOT, "output-data")
RESULTS = file.path(ROOT, "results")
dir.create(RESULTS, showWarnings = FALSE)
for (f in c(BASE, IMPR)) stopifnot(dir.exists(file.path(f, "summarizing_results_TCPD")))
report = c("Baseline (original scoth code)  vs  improved (scoth.fast, fast exact scoring)",
           strrep("=", 78), paste("baseline folder:", BASE), paste("improved folder:", IMPR), "")

## 1. Scores ------------------------------------------------------------------------
files = list.files(file.path(BASE, "algorithm_scoring_results"))
stopifnot(identical(files, list.files(file.path(IMPR, "algorithm_scoring_results"))))
score_cmp = bind_rows(lapply(files, function(f) {
  b = readRDS(file.path(BASE, "algorithm_scoring_results", f))
  i = readRDS(file.path(IMPR, "algorithm_scoring_results", f))
  stopifnot(identical(b$cost_method, i$cost_method))
  bind_rows(lapply(seq_len(nrow(b)), function(r) {
    sb = b$scoring[[r]]$score
    si = i$scoring[[r]]$score
    data.frame(dataset = b$dataset_name[r], cost = b$cost_method[r], T = length(sb) + 1,
               cost_initial_identical = identical(b$cost_initial[r], i$cost_initial[r]),
               bit_identical = identical(sb, si), max_abs_score_diff = max(abs(sb - si)),
               baseline_seconds = as.numeric(b$scoring_time_taken[r]),
               improved_seconds = as.numeric(i$scoring_time_taken[r]))
  }))
}))
write.csv(score_cmp, file.path(RESULTS, "tcpd_score_differences.csv"), row.names = FALSE)
w_time = wilcox.test(score_cmp$baseline_seconds, score_cmp$improved_seconds, paired = TRUE, exact = FALSE)
report = c(report, "1. Scores (Algorithm 1), 42 datasets x 2 costs = 84 runs",
  sprintf("   V(empty) identical in %d/84 runs; scores bit-identical in %d/84 runs", sum(score_cmp$cost_initial_identical), sum(score_cmp$bit_identical)),
  sprintf("   max |score difference| over all runs = %.3g", max(score_cmp$max_abs_score_diff)),
  sprintf("   scoring time (sum of the 84 runs): baseline %.1f s, improved %.2f s (speed-up x%.0f)",
          sum(score_cmp$baseline_seconds), sum(score_cmp$improved_seconds),
          sum(score_cmp$baseline_seconds) / sum(score_cmp$improved_seconds)),
  "   (times of two different pipeline runs; same-process comparison: results/tcpd_scoring_time.csv)",
  sprintf("   paired Wilcoxon signed-rank test on the 84 scoring times: V = %.0f, p = %.3g",
          unname(w_time$statistic), w_time$p.value), "")

## 2. Selections --------------------------------------------------------------------
same_partition_values = function(u, v) {    # s, e and cost of every level
  if (is.null(u$partitions) || is.null(v$partitions)) return(identical(u$partitions, v$partitions))
  length(u$partitions) == length(v$partitions) &&
    all(mapply(function(p, q) identical(p$s, q$s) && identical(p$e, q$e) && identical(unname(p$cost), unname(q$cost)),
               u$partitions, v$partitions))
}
sel_cmp = bind_rows(lapply(files, function(f) {
  b = readRDS(file.path(BASE, "algorithm_selection_results", f))
  i = readRDS(file.path(IMPR, "algorithm_selection_results", f))
  stopifnot(identical(b$cost_method, i$cost_method),
            identical(b$selection_threshold_included, i$selection_threshold_included),
            identical(b$selection_level, i$selection_level))
  data.frame(dataset = b$dataset_name, cost = b$cost_method,
             same_selection = mapply(identical, b$selection, i$selection),
             same_partitions = mapply(same_partition_values, b$selection_additional, i$selection_additional)) %>%
    group_by(dataset, cost) %>%
    summarize(n_rows = n(), n_identical_selections = sum(same_selection),
              n_identical_partitions = sum(same_partitions), .groups = "drop")
}))
write.csv(sel_cmp, file.path(RESULTS, "tcpd_selection_identity.csv"), row.names = FALSE)
report = c(report, "2. Selections (Algorithm 2): every (dataset, cost, threshold, level) row of the pipeline",
  sprintf("   identical selections (0/1 vector of change points): %d / %d rows",
          sum(sel_cmp$n_identical_selections), sum(sel_cmp$n_rows)),
  sprintf("   identical partitions (s, e, segment cost of every level): %d / %d rows",
          sum(sel_cmp$n_identical_partitions), sum(sel_cmp$n_rows)), "")

## 3. Evaluation tables ---------------------------------------------------------------
keys = c("dataset_name", "feature_name", "cost_method", "selection_threshold_included", "selection_level", "groundtruth_name")
metric_cmp = bind_rows(lapply(c("F1", "TCPD"), function(metric) {
  vals = if (metric == "F1") c("F1", "nb_gt", "nb_pred") else c("F1", "covering", "nb_gt", "nb_pred")
  b = readRDS(file.path(BASE, paste0("summarizing_results_", metric), "summary.RDS")) %>% select(all_of(c(keys, vals)))
  i = readRDS(file.path(IMPR, paste0("summarizing_results_", metric), "summary.RDS")) %>% select(all_of(c(keys, vals)))
  j = inner_join(b, i, by = keys, suffix = c("_base", "_impr"))
  bind_rows(lapply(vals, function(v) data.frame(
    evaluation = paste0("summarizing_results_", metric), column = v,
    n_rows_baseline = nrow(b), n_rows_improved = nrow(i), n_rows_matched = nrow(j),
    n_rows_identical = sum(mapply(identical, j[[paste0(v, "_base")]], j[[paste0(v, "_impr")]])),
    max_abs_diff = max(abs(j[[paste0(v, "_base")]] - j[[paste0(v, "_impr")]]), na.rm = TRUE))))
}))
write.csv(metric_cmp, file.path(RESULTS, "tcpd_metric_differences.csv"), row.names = FALSE)
report = c(report, "3. Evaluation tables (every dataset x cost x threshold x level x annotator row)",
  paste0("   ", capture.output(print(metric_cmp[, c("evaluation", "column", "n_rows_matched", "n_rows_identical", "max_abs_diff")], row.names = FALSE))), "")

## 4. The three protocols -------------------------------------------------------------
# same code as experiments/extra_protocols.R
protocols = function(folder) {
  res = readRDS(file.path(folder, "summarizing_results_F1", "summary.RDS")) %>%
    filter(groundtruth_name == "groundtruth") %>%
    rename(threshold = selection_threshold_included, level = selection_level) %>%
    select(dataset_name, cost_method, threshold, level, F1) %>%
    filter(is.finite(threshold), threshold > 0)
  oracle = res %>% group_by(dataset_name) %>% summarize(F1 = max(F1), .groups = "drop")
  per_ds_best_level = res %>% group_by(dataset_name, cost_method, threshold) %>% summarize(F1 = max(F1), .groups = "drop")
  bs_table = per_ds_best_level %>% group_by(cost_method, threshold) %>%
    summarize(mean_F1 = mean(F1), sd_F1 = sd(F1), .groups = "drop") %>% arrange(-mean_F1)
  bs = per_ds_best_level %>% filter(cost_method == bs_table$cost_method[1], threshold == bs_table$threshold[1])
  fa_table = res %>% group_by(cost_method, threshold, level) %>%
    summarize(mean_F1 = mean(F1), sd_F1 = sd(F1), .groups = "drop") %>% arrange(-mean_F1)
  fa = res %>% filter(cost_method == fa_table$cost_method[1], threshold == fa_table$threshold[1], level == fa_table$level[1])
  datasets = sort(unique(res$dataset_name))
  loo = sapply(datasets, function(d) {
    train = res %>% filter(dataset_name != d) %>% group_by(cost_method, threshold, level) %>%
      summarize(m = mean(F1), .groups = "drop") %>% arrange(-m, -threshold, level)
    ch = train[1, ]
    res %>% filter(dataset_name == d, cost_method == ch$cost_method, threshold == ch$threshold, level == ch$level) %>% pull(F1)
  })
  # oracle F1-biased and cover: the package function used for table_left.txt
  f1b = scoth.fast:::get_oracle("F1_biased", folder, full_results = TRUE) %>% select(dataset_name, F1_biased = eval)
  cov = scoth.fast:::get_oracle("cover", folder, full_results = TRUE) %>% select(dataset_name, cover = eval)
  per_dataset = oracle %>% rename(oracle_F1 = F1) %>%
    left_join(bs %>% select(dataset_name, best_single_F1 = F1), by = "dataset_name") %>%
    left_join(fa %>% select(dataset_name, fully_auto_F1 = F1), by = "dataset_name") %>%
    left_join(tibble(dataset_name = datasets, fully_auto_loo_F1 = as.numeric(loo)), by = "dataset_name") %>%
    left_join(f1b %>% rename(oracle_F1_biased = F1_biased), by = "dataset_name") %>%
    left_join(cov %>% rename(oracle_cover = cover), by = "dataset_name")
  settings = c(oracle_F1 = "per-dataset (cost, threshold, level)",
               best_single_F1 = sprintf("%s; threshold=%s; level per dataset", bs_table$cost_method[1], bs_table$threshold[1]),
               fully_auto_F1 = sprintf("%s; threshold=%s; level=%d", fa_table$cost_method[1], fa_table$threshold[1], fa_table$level[1]),
               fully_auto_loo_F1 = "chosen on the other 41 datasets",
               oracle_F1_biased = "per-dataset (table_left.txt)", oracle_cover = "per-dataset (table_left.txt)")
  list(per_dataset = per_dataset, settings = settings)
}
p_base = protocols(BASE)
p_impr = protocols(IMPR)
stopifnot(identical(p_base$per_dataset$dataset_name, p_impr$per_dataset$dataset_name))
cols = setdiff(names(p_base$per_dataset), "dataset_name")
per_ds = p_base$per_dataset["dataset_name"]
for (cl in cols) {
  per_ds[[paste0(cl, "_baseline")]] = p_base$per_dataset[[cl]]
  per_ds[[paste0(cl, "_improved")]] = p_impr$per_dataset[[cl]]
  per_ds[[paste0(cl, "_difference")]] = p_impr$per_dataset[[cl]] - p_base$per_dataset[[cl]]
}
write.csv(per_ds, file.path(RESULTS, "tcpd_per_dataset_F1_differences.csv"), row.names = FALSE)

fmt = function(v) sprintf("%.2f±%.2f", mean(v), sd(v))
prot = bind_rows(lapply(cols, function(cl) {
  b = p_base$per_dataset[[cl]]; i = p_impr$per_dataset[[cl]]; d = i - b
  wt = if (all(d == 0)) NULL else wilcox.test(i, b, paired = TRUE, exact = FALSE)
  data.frame(protocol = cl, baseline = fmt(b), improved = fmt(i),
             baseline_mean = mean(b), improved_mean = mean(i), mean_difference = mean(d),
             n_datasets_different = sum(d != 0),
             baseline_setting = p_base$settings[[cl]], improved_setting = p_impr$settings[[cl]],
             wilcoxon_p_value = if (is.null(wt)) NA_real_ else wt$p.value,
             wilcoxon_note = if (is.null(wt)) "all 42 differences are exactly 0: test not applicable (identical results)" else "paired, two-sided")
}))
write.csv(prot, file.path(RESULTS, "tcpd_protocols_comparison.csv"), row.names = FALSE)
write.csv(prot %>% select(protocol, n_datasets_different, wilcoxon_p_value, wilcoxon_note) %>%
            bind_rows(data.frame(protocol = "scoring time per run (84 runs)", n_datasets_different = sum(score_cmp$baseline_seconds != score_cmp$improved_seconds),
                                 wilcoxon_p_value = w_time$p.value, wilcoxon_note = "paired, two-sided, baseline vs improved seconds")),
          file.path(RESULTS, "tcpd_wilcoxon.csv"), row.names = FALSE)
report = c(report, "4. Protocols on the 42 datasets (mean±sd; F1 with margin 5, median annotator)",
  paste0("   ", capture.output(print(prot[, c("protocol", "baseline", "improved", "n_datasets_different", "improved_setting")], row.names = FALSE))),
  "", sprintf("   paired Wilcoxon signed-rank tests: %s", paste(unique(prot$wilcoxon_note), collapse = "; ")), "")

## 5. Stage timings (parsed from the logs of the two runs) -------------------------------
parse_stages = function(log_file) {
  if (!file.exists(log_file)) return(NULL)
  l = readLines(log_file, warn = FALSE)
  l = l[grepl("^INFO \\[", l)]
  ts = as.POSIXct(sub("^INFO \\[([^]]+)\\].*$", "\\1", l), format = "%Y-%m-%d %H:%M:%S")
  at = function(pattern) ts[grep(pattern, l)[1]]
  marks = c(scoring = at("Begin the scoring step"), selection = at("Begin the selection step"),
            evaluation_F1 = at("Begin the evaluation step for the F1 metric"),
            summarization_F1 = at("Begin the summarization step for the F1 metric"),
            evaluation_TCPD = at("Begin the evaluation step for the TCPD metric"),
            summarization_TCPD = at("Begin the summarization step for the TCPD metric"),
            tables_and_figure = at("Result for our methodology in the oracle setting"), end = tail(ts, 1))
  data.frame(stage = c(names(marks)[-length(marks)], "total"),
             seconds = c(as.numeric(diff(marks), units = "secs"), as.numeric(difftime(marks[["end"]], marks[["scoring"]], units = "secs"))))
}
tb = parse_stages(file.path(BASE_ROOT, "logs", "main_run.log"))
ti = parse_stages(file.path(ROOT, "logs", "main_fast_run.log"))
# the improved run also has output-data/stage_timings.csv (written by main_fast.R,
# 0.01 s resolution); the baseline (experiments/main.R of project0) only has its log
precise = function(file) {
  if (!file.exists(file)) return(data.frame(stage = character(0), improved_seconds_stage_timings = numeric(0)))
  st = read.csv(file)
  to_log_stage = c("1_scoring" = "scoring", "2_selection" = "selection", "3_evaluation_F1" = "evaluation_F1",
                   "4_summarization_F1" = "summarization_F1", "5_evaluation_TCPD" = "evaluation_TCPD",
                   "6_summarization_TCPD" = "summarization_TCPD", "7_tables" = "tables_and_figure",
                   "8_figure" = "tables_and_figure", "total" = "total")
  st %>% mutate(stage = unname(to_log_stage[stage])) %>% group_by(stage) %>%
    summarize(improved_seconds_stage_timings = round(sum(seconds), 2), .groups = "drop")
}
if (!is.null(tb) && !is.null(ti)) {
  timing = tb %>% rename(baseline_seconds = seconds) %>% left_join(ti %>% rename(improved_seconds = seconds), by = "stage") %>%
    left_join(precise(file.path(IMPR, "stage_timings.csv")), by = "stage") %>%
    mutate(baseline_minutes = round(baseline_seconds / 60, 2), improved_minutes = round(improved_seconds / 60, 2),
           # ratio of the two logs only when BOTH stages last >= 10 s: with a 1 s clock a
           # "1 s" stage may last anything from 0 to 2 s (the scoring stage of the improved run)
           speedup_from_logs = ifelse(baseline_seconds >= 10 & improved_seconds >= 10,
                                      round(baseline_seconds / improved_seconds, 1), NA),
           # baseline log time / improved stage_timings.csv time (baseline stage >= 10 s)
           speedup_with_stage_timings = ifelse(baseline_seconds >= 10,
                                               round(baseline_seconds / improved_seconds_stage_timings, 1), NA),
           changed_by_improvement5 = stage %in% c("scoring", "selection"))
  write.csv(timing, file.path(RESULTS, "tcpd_stage_timings.csv"), row.names = FALSE)
  report = c(report, "5. Wall-clock time of each stage of the two pipeline runs (baseline_seconds, improved_seconds: INFO lines",
             "   of the two logs, 1 s resolution; improved_seconds_stage_timings: output-data/stage_timings.csv, 0.01 s).",
             "   Only scoring and selection are changed, the other stages are the same code. Two different runs: the",
             "   machine load differs; the same-process comparison of the scoring stage is results/tcpd_scoring_time.csv.",
             paste0("   ", capture.output(print(timing, row.names = FALSE))), "")
}

report = c(report, "Baseline table_left.txt / table_right.txt:", readLines(file.path(BASE, "table_left.txt")), readLines(file.path(BASE, "table_right.txt")),
           "", "Improved table_left.txt / table_right.txt:", readLines(file.path(IMPR, "table_left.txt")), readLines(file.path(IMPR, "table_right.txt")))
writeLines(report, file.path(RESULTS, "comparison_summary.txt"))
cat(report, sep = "\n")
