## Focused Scan

# Setup
setwd("[PATH TO FOCUSED SCAN OUTPUT DIRECTORY]")
f2_dir <- "[PATH TO F2 STATS DIRECTORY]"
populations <- c("[POPULATION NAMES]")
out_group <- "[OUTGROUP NAME]"
iterations <- [NUMBER OF ITERATIONS TO RUN FIND_GRAPHS]
num_admix <- [NUMBER OF ADMIXTURE EVENTS TO RUN FIND_GRAPHS ON]

library("foreach")
library("doParallel")
library("admixtools")
library("dplyr")

# Obtain f-stats for required pops
f2_stats <- f2_from_precomp(f2_dir, pops = populations)

# Setting up parallel computing
core <- detectCores()
cl <- makeCluster(core)
registerDoParallel(cl)

cat(paste("\nFocused Scan Log -", num_admix, "Admixture Events -", Sys.time(), "\n"), file = "focused-scan-log.txt", append = TRUE)

# Run FindGraphs in parallel for a set number of iterations with a specific number of admixture events
cumulative_results <- foreach(i = 1:iterations, .combine = "c", .packages = c("admixtools", "dplyr"), .inorder = TRUE) %dopar% {
    # Include which iteration is being started
    sink("focused-scan-log.txt", append = TRUE)
    cat(paste("----- STARTING ITERATION", i, "/", iterations, "-----\n"))

    # Run findGraphs
    res <- find_graphs(f2_stats, numadmix = num_admix, outpop = out_group, stop_gen = 10000, stop_gen2 = 30, numgraphs = 10, plusminus_generations = 10)
    min_scores <- res %>% slice_min(score)

    # Run qpGraph to obtain worst residual score and add to row
    worst_res <- qpgraph(f2_stats, min_scores$graph[[1]], return_fstats = TRUE, , numstart = 1000)
    min_scores$worst_residual <- worst_res$worst_residual

    cat(paste("----- COMPLETED ITERATION", i, "/", iterations, "-----\n"))
    sink()

    return(list(min_scores))
}

# Save focused scan results tibble
focused_scan_results <- bind_rows(cumulative_results)
focused_scan_results

saveRDS(focused_scan_results, file = "focused-scan-result.rds")

stopCluster(cl)
