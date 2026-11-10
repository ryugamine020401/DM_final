# =============================================================================
# install_packages.R
# Installs the CRAN packages needed by project0_reproduction and project5_fast_scoring.
# Usage (from the repository root):  Rscript setup\install_packages.R
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c("dplyr", "ggplot2", "gridExtra", "here", "jsonlite", "logger",
          "testthat", "tibble", "tidyr", "pkgload")

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, type = "binary", dependencies = c("Depends", "Imports", "LinkingTo"))
} else {
  message("Already installed: ", paste(pkgs, collapse = ", "))
}

# Final check
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = pkgs, installed = ok, row.names = NULL))
if (!all(ok)) stop("Some packages could not be installed: ", paste(pkgs[!ok], collapse = ", "))
cat("All requested packages are installed.\n")
