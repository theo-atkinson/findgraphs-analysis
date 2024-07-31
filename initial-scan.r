## Initial Scan

# Setup
setwd("[PATH TO INITIAL SCAN OUTPUT DIRECTORY]")
f2_dir <- "[PATH TO F2 STATS DIRECTORY]"
populations <- c("[POPULATION NAMES]")
out_group <- "[OUTGROUP NAME]"
iterations <- [NUMBER OF ITERATIONS TO RUN FIND_GRAPHS]
max_admix <- [MAXIMUM NUMBER OF ADMIXTURE EVENTS TO RUN FIND_GRAPHS ON]

# Load required packages
library("foreach")
library("doParallel")
library("admixtools")
library("dplyr")
library("ggplot2")

# Obtain f-stats for required pops
f2_stats <- f2_from_precomp(f2_dir, pops = populations)

# Setting up parallel computing
core <- detectCores()
cl <- makeCluster(core)
registerDoParallel(cl)

cat(paste("\nInitial Scan Log -", Sys.time(), "\n"), file = "initial-scan-log.txt", append = TRUE)

# Run FindGraphs in parallel multiple times for various numbers of admixture events
cumulative_results <- foreach(i = 1:iterations, .combine = "c", .inorder = TRUE) %:%
    foreach(admixture_events = 1:max_admix, .combine = "c", .packages = c("admixtools", "dplyr"), .inorder = TRUE) %dopar% {
        # Include which iteration is being started to log file
        sink("initial-scan-log.txt", append = TRUE)
        cat(paste("----- STARTING ITERATION", i, "/", iterations, "FOR", admixture_events, "ADMIXTURE EVENTS -----\n"))

        # Run findGraphs and select best graph
        res <- find_graphs(f2_stats, numadmix = admixture_events, outpop = out_group, stop_gen = 10000, stop_gen2 = 30, numgraphs = 10, plusminus_generations = 10)
        min_scores <- res %>% slice_min(score)

        # Add the number of admixture_events to the row
        min_scores$admixture_events <- admixture_events

        # Run qpGraph to obtain worst residual score and also add to row
        worst_res <- qpgraph(f2_stats, min_scores$graph[[1]], return_fstats = TRUE, , numstart = 1000)
        min_scores$worst_residual <- worst_res$worst_residual

        cat(paste("----- COMPLETED ITERATION", i, "/", iterations, "FOR", admixture_events, "ADMIXTURE EVENTS -----\n"))

        return(list(min_scores))
    }
initial_scan_results <- bind_rows(cumulative_results)
initial_scan_results

# Save initial scan results tibble
saveRDS(initial_scan_results, file = "initial-scan-result.rds")

# Close parallel computing
stopCluster(cl)

# Plot score vs worst residual
gg <- ggplot(initial_scan_results, aes(x = as.factor(admixture_events), y = worst_residual)) +
    geom_point(position = position_jitter(w = 0.2, h = 0), alpha = 0.5) +
    stat_summary(fun = min, colour = "blue", geom = "line", aes(group = 1)) +
    scale_y_continuous(breaks = seq(0, max(initial_scan_results$worst_residual), by = 1)) +
    labs(
        x = "Number of Admixture Events",
        y = "Worst Residual"
    ) +
    theme_minimal()

ggsave("initial-scan-plot.jpg", plot = gg, width = 8, height = 6, dpi = 300)
