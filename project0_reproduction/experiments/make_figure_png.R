# 把論文 Fig. 4 另存成 PNG

suppressMessages({
  library(scoth)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

PUBLIC_DATA_FOLDER = file.path(here::here(), "public-data")
OUTPUT_DATA_FOLDER = file.path(here::here(), "output-data")   # read only
groundtruth = load_groundtruths(PUBLIC_DATA_FOLDER)

# The authors' function first draws the figure on the current graphics device (grid.arrange)
# and then saves it with ggsave(). Without an open device, Rscript would create an extra
# Rplots.pdf in the working folder; a null PDF device avoids that.
pdf(NULL)
plot_figure_qualitative_results(OUTPUT_DATA_FOLDER, groundtruth,
                                outfile = file.path(here::here(), "results", "figure_qualitative_results.png"))
invisible(dev.off())