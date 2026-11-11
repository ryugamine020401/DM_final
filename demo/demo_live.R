# demo_live.R -- live demo for the presentation (about 5-10 seconds)
# Run it from the repository root:   Rscript demo\demo_live.R
# It needs both packages: scoth (original, from project0) and scoth.fast (improvement 5).

# 1) A toy series with two obvious changes: after x[4] and after x[8]
x = c(0, 0.1, -0.1, 0, 5, 5.1, 4.9, 5, 1, 1.1, 0.9, 1)
o = scoth::scoth_scoring_func(x, "cusum_L2", verbose = FALSE)
f = scoth.fast::scoth_scoring_func(x, "cusum_L2", verbose = FALSE, algorithm = "fast")
cat("Toy series, T = 12, quadratic cost\n")
cat("  scores:", sprintf("%.3f", o$changepoint_gain), "\n")
cat("  fast version gives identical scores:", identical(o$changepoint_gain, f$changepoint_gain), "\n")
cat("  score >= 0.1 (level 1): change after x[",
    paste(which(o$changepoint_gain >= 0.1), collapse = "] and after x["), "]\n", sep = "")

# 2) A real TCPD series: measles (T = 991), linear cost
k = scoth.fast::load_kpis(file.path("project5_fast_scoring", "public-data"))
y = k$kpis[[which(k$dataset_name == "measles")]][, 1]
t_orig = system.time(o2 <- scoth::scoth_scoring_func(y, "reg_cusum_L2", verbose = FALSE))[["elapsed"]]
t_fast = system.time(f2 <- scoth.fast::scoth_scoring_func(y, "reg_cusum_L2", verbose = FALSE, algorithm = "fast"))[["elapsed"]]
cat(sprintf("\nTCPD measles, T = %d, linear cost\n", length(y)))
cat(sprintf("  original %.2f s | fast %.3f s | speed-up x%.0f\n", t_orig, t_fast, t_orig / t_fast))
cat(sprintf("  max |score difference| = %.2g\n", max(abs(o2$changepoint_gain - f2$changepoint_gain))))
options(scoth.fast_exact_max_len = Inf)
f3 = scoth.fast::scoth_scoring_func(y, "reg_cusum_L2", verbose = FALSE, algorithm = "fast")
cat("  bit-identical with exact_max_len = Inf:", identical(o2$changepoint_gain, f3$changepoint_gain), "\n")
