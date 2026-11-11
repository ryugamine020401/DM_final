# 跑原作者的單元測試


suppressMessages(library(dplyr))
res = testthat::test_local(".", reporter = "summary", stop_on_failure = FALSE)
df = as.data.frame(res)
cat(sprintf("\ntests: %d   expectations: %d   failed: %d   errors: %d\n",
            nrow(df), sum(df$nb), sum(df$failed), sum(df$error)))
if (sum(df$failed) > 0 || any(df$error)) quit(status = 1)