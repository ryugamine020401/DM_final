# 算三種評估協定



suppressMessages({ library(dplyr); library(tidyr) })

OUTPUT_DATA_FOLDER = file.path(here::here(), "output-data")
summary_file = file.path(OUTPUT_DATA_FOLDER, "summarizing_results_F1", "summary.RDS")
res = readRDS(summary_file) %>%
  filter(groundtruth_name == "groundtruth") %>%            # median annotator only
  rename(threshold = selection_threshold_included, level = selection_level) %>%
  select(dataset_name, cost_method, threshold, level, F1)
# keep the 24 thresholds of the paper (drop the two degenerate ones added by the pipeline: 0 and +Inf)
res = res %>% filter(is.finite(threshold), threshold > 0)
n_ds = n_distinct(res$dataset_name)

fmt = function(v) sprintf("%.2f±%.2f", mean(v), sd(v))

## (1) Oracle -------------------------------------------------------------------
oracle = res %>% group_by(dataset_name) %>% summarize(F1 = max(F1), .groups = "drop")

## (2) Best-single (paper) ------------------------------------------------------
per_ds_best_level = res %>% group_by(dataset_name, cost_method, threshold) %>%
  summarize(F1 = max(F1), .groups = "drop")
bs_table = per_ds_best_level %>% group_by(cost_method, threshold) %>%
  summarize(mean_F1 = mean(F1), sd_F1 = sd(F1), .groups = "drop") %>% arrange(-mean_F1)
bs_best = bs_table[1, ]
best_single = per_ds_best_level %>% filter(cost_method == bs_best$cost_method, threshold == bs_best$threshold)

## (3) Fully-automatic: one (cost, threshold, level) for every dataset -----------
fa_table = res %>% group_by(cost_method, threshold, level) %>%
  summarize(mean_F1 = mean(F1), sd_F1 = sd(F1), .groups = "drop") %>% arrange(-mean_F1)
fa_best = fa_table[1, ]
fully_auto = res %>% filter(cost_method == fa_best$cost_method, threshold == fa_best$threshold, level == fa_best$level)

# best fixed level for the paper's recommended setting (linear cost, threshold 0.1)
paper_setting = fa_table %>% filter(cost_method == "reg_cusum_L2", abs(threshold - 0.1) < 1e-9) %>% arrange(level)

# leave-one-dataset-out CV of the fully-automatic global choice
datasets = sort(unique(res$dataset_name))
loo = sapply(datasets, function(d) {
  train = res %>% filter(dataset_name != d) %>% group_by(cost_method, threshold, level) %>%
    summarize(m = mean(F1), .groups = "drop") %>% arrange(-m, -threshold, level)
  ch = train[1, ]
  res %>% filter(dataset_name == d, cost_method == ch$cost_method, threshold == ch$threshold, level == ch$level) %>% pull(F1)
})

## Report -------------------------------------------------------------------------
lines = c(
  sprintf("Datasets: %d (F1 margin 5, median annotator)", n_ds),
  "",
  sprintf("(1) ORACLE (per-dataset cost, threshold, level)           F1 = %s", fmt(oracle$F1)),
  sprintf("(2) BEST-SINGLE (global cost+threshold, per-dataset level)  F1 = %s   [%s; threshold=%s]",
          fmt(best_single$F1), bs_best$cost_method, bs_best$threshold),
  sprintf("(3) FULLY-AUTOMATIC in-sample (global cost+threshold+level) F1 = %s   [%s; threshold=%s; level=%d]",
          fmt(fully_auto$F1), fa_best$cost_method, fa_best$threshold, fa_best$level),
  sprintf("(3') FULLY-AUTOMATIC leave-one-dataset-out CV                F1 = %s", fmt(loo)),
  "",
  "Paper setting (linear cost = reg_cusum_L2, threshold = 0.1): mean F1 for each FIXED level",
  paste0("  level ", paper_setting$level, ": ", sprintf("%.3f", paper_setting$mean_F1)),
  "",
  "Top-10 fully-automatic settings (cost, threshold, level):",
  capture.output(print(as.data.frame(head(fa_table, 10)), row.names = FALSE))
)
writeLines(lines, file.path(OUTPUT_DATA_FOLDER, "extra_protocols.txt"))
cat(paste(lines, collapse = "\n"), "\n")

per_dataset = oracle %>% rename(F1_oracle = F1) %>%
  left_join(best_single %>% select(dataset_name, F1_best_single = F1), by = "dataset_name") %>%
  left_join(fully_auto %>% select(dataset_name, F1_fully_auto = F1), by = "dataset_name") %>%
  left_join(tibble(dataset_name = datasets, F1_fully_auto_loo = as.numeric(loo)), by = "dataset_name")
write.csv(per_dataset, file.path(OUTPUT_DATA_FOLDER, "extra_protocols_per_dataset.csv"), row.names = FALSE)