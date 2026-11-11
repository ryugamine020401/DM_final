# =============================================================================
# main_fast.R   [IMPROVEMENT-5]  (adapted from experiments/main.R)
#
# Full TCPD pipeline of the paper (42 datasets, 2 costs, 24 thresholds, 10 levels)
# with the FAST EXACT scoring (Algorithm 1) and the fast exact selection path
# (Algorithm 2). Everything else (evaluation, tables, figure) is the original code.
#
# Usage (from the project folder project5_fast_scoring):
#   Rscript experiments/main_fast.R > logs/main_fast_run.log 2>&1
# Outputs: output-data/ (same structure as the baseline) and
#          output-data/stage_timings.csv (wall-clock time of each stage).
# Run time: about 9-15 minutes, depending on the load of the machine (most of it is
# the TCPD-metric evaluation, unchanged).
# =============================================================================
library(scoth.fast)          # [IMPROVEMENT-5] the improved package (NOT library(scoth))
library(dplyr)
library(tidyr)
library(tibble)
library(jsonlite)
library(logger)
# libraries used for the figure
library(ggplot2)
library(gridExtra)
library(grid)

## [IMPROVEMENT-5] switch on the fast exact algorithms -------------------------
options(scoth.scoring_algorithm = "fast",  # Algorithm 1 in O(T log T) (scoth_changes_fast)
        scoth.fast_exact_max_len = 64,     # segments <= 64 points: bit-identical replica costs
        scoth.fast_selection = TRUE)       # Algorithm 2: replica costs, no unused df2, vectorized zoom

## Parameters ----
PUBLIC_DATA_FOLDER = file.path(here::here(), "public-data")
OUTPUT_DATA_FOLDER = file.path(here::here(), "output-data")
cost_methods = c("reg_cusum_L2", "cusum_L2")
thresholds = c(seq(from = 0.03, to = 0.2, by = 0.01), 0.3, 0.4, 0.5, 0.6, 0.7, 1)
nb_levels = 10

# [IMPROVEMENT-5] the pipeline skips the files that already exist: start from an
# empty output folder so that the measured times are complete.
if (dir.exists(OUTPUT_DATA_FOLDER)) {
  message("Removing the previous ", OUTPUT_DATA_FOLDER)
  unlink(OUTPUT_DATA_FOLDER, recursive = TRUE)
}
dir.create(OUTPUT_DATA_FOLDER, showWarnings = FALSE, recursive = TRUE)

# [IMPROVEMENT-5] helper that runs one stage and records its wall-clock time
timings = data.frame(stage = character(0), seconds = numeric(0))
run_stage = function(stage, expr) {
  t0 = Sys.time()
  force(expr)
  sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings[nrow(timings) + 1, ] <<- list(stage, round(sec, 2))
  logger::log_info(sprintf("[IMPROVEMENT-5] stage '%s' took %.1f s", stage, sec))
  invisible(NULL)
}

## Loading benchmark data and ground truth ----
loaded_kpis = load_kpis(PUBLIC_DATA_FOLDER)
groundtruth = load_groundtruths(PUBLIC_DATA_FOLDER)

## Pipeline ----
run_stage("1_scoring", pipeline_scoring_func(loaded_kpis, cost_methods, OUTPUT_DATA_FOLDER))
run_stage("2_selection", pipeline_selection_func(thresholds, nb_levels, OUTPUT_DATA_FOLDER))

metric = "F1" # F1 metric
run_stage("3_evaluation_F1", pipeline_evaluation_func(metric, groundtruth, OUTPUT_DATA_FOLDER))
run_stage("4_summarization_F1", pipeline_summarization_func(metric, OUTPUT_DATA_FOLDER))

metric = "TCPD" # F1-biased and cover metrics from the TCPD paper
run_stage("5_evaluation_TCPD", pipeline_evaluation_func(metric, groundtruth, OUTPUT_DATA_FOLDER))
run_stage("6_summarization_TCPD", pipeline_summarization_func(metric, OUTPUT_DATA_FOLDER))

## Result tables ----
run_stage("7_tables", {
  get_oracle_table(OUTPUT_DATA_FOLDER) # table left
  get_single_best_table(OUTPUT_DATA_FOLDER) # table right
})

## Result figure ----
run_stage("8_figure", plot_figure_qualitative_results(OUTPUT_DATA_FOLDER, groundtruth)) # last figure

timings[nrow(timings) + 1, ] = list("total", round(sum(timings$seconds), 2))
write.csv(timings, file.path(OUTPUT_DATA_FOLDER, "stage_timings.csv"), row.names = FALSE)
print(timings)
cat(readLines(file.path(OUTPUT_DATA_FOLDER, "table_left.txt")), sep = "\n")
cat(readLines(file.path(OUTPUT_DATA_FOLDER, "table_right.txt")), sep = "\n")
