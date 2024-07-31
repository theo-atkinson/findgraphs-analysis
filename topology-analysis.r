## Topology Analysis

library("foreach")
library("doParallel")
library("admixtools")
library("dplyr")
library("ggplot2")
library("data.table")

# Setup
setwd("[PATH TO TOPOLOGY ANALYSIS OUTPUT DIRECTORY]")
f2_dir <- "[PATH TO F2 STATS DIRECTORY]"
populations <- c("[POPULATION NAMES]")
focused_scan <- readRDS("[PATH TO FOCUSED SCAN .RDS FILE]")
max_tops <- [MAXIMUM NUMBER OF BEST TOPOLOGIES TO BE ANALYSED]
bootstraps <- [NUMBER OF BOOTSTRAP ITERATIONS]
admix_constraints <- tribble( # OPTIONAL, CAN LEAVE LIKE THIS
    ~pop, ~min, ~max,
)
event_constraints <- tribble( # OPTIONAL, CAN LEAVE LIKE THIS
    ~earlier1, ~earlier2, ~later1, ~later2,
)

# Obtain f-stats for required pops
f2_stats <- f2_from_precomp(f2_dir, pops = populations)

cat(paste("\nTopology Analysis -", Sys.time(), "\n"), file = "topology-analysis-log.txt", append = TRUE)

# Apply admixture constraints
if (!file.exists("focused-scan-constrained.RDS")) {
    constrained_list <- data.frame()
    for (i in seq_along(focused_scan$graph)) {
        if (satisfies_numadmix(focused_scan$graph[[i]], admix_constraints)) {
            constrained_list <- rbind(constrained_list, focused_scan[i, ])
        } else {
            if (nrow(event_constraints) > 1 && satisfies_eventorder(focused_scan$graph[[i]], event_constraints)) {
                constrained_list <- rbind(constrained_list, focused_scan[i, ])
            }
        }
    }
    focused_scan <- constrained_list
    focused_scan <- arrange(focused_scan, score) %>% mutate(row_id = row_number())

    saveRDS(focused_scan, file = "focused-scan-constrained.RDS")

    cat(paste("Constraints Applied (if applicable)\n"), file = "topology-analysis-log.txt", append = TRUE)
}

# 2. Plot score vs worst residual
if (!file.exists("score-WR-plot.jpg")) {
    focused_scan <- readRDS("focused-scan-constrained.RDS")
    gg <- ggplot(focused_scan, aes(x = score, y = worst_residual)) +
        geom_point(alpha = 0.2) +
        coord_cartesian(xlim = c(0, max(focused_scan$score + 2)), ylim = c(0, max(focused_scan$worst_residual + 0.2))) +
        theme_minimal()

    ggsave("score-WR-plot.jpg", plot = gg, width = 8, height = 6, dpi = 300)

    cat(paste("Score - Worst Residual Scatter Plotted\n"), file = "topology-analysis-log.txt", append = TRUE)
}

# 3. Plot histogram of result scores
if (!file.exists("scores-hist.jpg")) {
    focused_scan <- readRDS("focused-scan-constrained.RDS")
    gg <- ggplot(focused_scan, aes(x = score)) +
        geom_histogram(stat = "bin", binwidth = 0.1) +
        coord_cartesian(xlim = c(0, max(focused_scan$score + 2))) +
        theme_minimal()

    ggsave("scores-hist.jpg", plot = gg, width = 8, height = 6, dpi = 300)

    cat(paste("Score Histogram Plotted\n"), file = "topology-analysis-log.txt", append = TRUE)
}

# 4. Identifying graphs with the same topology
if (!file.exists("topology-list.csv")) {
    focused_scan <- readRDS("focused-scan-constrained.RDS")
    identical <- isomorphism_classes2(focused_scan$graph)
    focused_scan$topology <- identical
    top_list <- data.frame(focused_scan$row_id, focused_scan$topology, focused_scan$score, focused_scan$worst_residual)
    top_list <- arrange(top_list, focused_scan.topology)
    write.table(top_list, "topology-list.csv", sep = ",", quote = FALSE, row.names = FALSE)
    saveRDS(focused_scan, "focused-scan-constrained.RDS")
    cat(paste("Identical Graphs Identified and Topology List Saved\n"), file = "topology-analysis-log.txt", append = TRUE)
}

# 5. Function to perform topology identification and plotting
topology_identification <- function(top_num) {
    topology <- data.frame()
    for (i in seq_along(focused_scan$graph)) {
        if (focused_scan$topology[[i]] == top_num) {
            topology <- rbind(topology, focused_scan[i, ])
        }
    }

    dir.create(paste0("graph-topology-", top_num))
    identified_graphs <- data.frame(topology$row_id, topology$score, topology$worst_residual)
    write.table(identified_graphs, paste0("graph-topology-", top_num, "/graph-topology-", top_num, ".csv"), sep = ", ", quote = FALSE, row.names = FALSE)

    core <- detectCores()
    cl <- makeCluster(core)
    registerDoParallel(cl)
    foreach(i = seq_along(topology$graph), .packages = c("ggplot2", "admixtools"), .inorder = TRUE) %dopar% {
        graph <- plot_graph(topology$edges[[i]])
        ggsave(paste0("graph-topology-", top_num, "/graph-topology-", top_num, "(", topology$row_id[[i]], ").pdf"), plot = graph)
    }
    stopCluster(cl)

    cat(paste0(nrow(topology), " Graph(s) With Topology ", top_num, " Identified and Plotted \n"), file = "topology-analysis-log.txt", append = TRUE)
}

# Run loop of topology identification function for each topology
focused_scan <- readRDS("focused-scan-constrained.RDS")
for (i in unique(focused_scan$topology)) {
    if (!dir.exists(paste0("graph-topology-", i))) {
        topology_identification(i)
    }
    if (i == max_topologies) {
        break
    }
}

# 6. Function to run and save graph fits from bootstrap analysis
graph_fit <- function(top_num) {
    set.seed(16121998)
    graph_num <- focused_scan[which(focused_scan$topology == top_num, arr.ind = TRUE)[1], ]
    fit <- qpgraph_resample_multi(f2_stats, list(focused_scan$graph[[graph_num$row_id]]), nboot = bootstraps, f3basepop = "Jalla")
    saveRDS(fit, file = paste0("graph-fits/fit-top-", top_num, ".rds"))
}

# Run loop to calculate bootstrap fits for best graph of each topology
if (!file.exists("graph-fits")) {
    dir.create("graph-fits")
}
core <- detectCores()
cl <- makeCluster(core)
registerDoParallel(cl)
foreach(i = c(1, (focused_scan %>% distinct(topology) %>% filter(topology < max_topologies))$topology + 1), .packages = c("admixtools"), .inorder = TRUE) %dopar% {
    if (!file.exists(paste0("graph-fits/fit-top-", i, ".rds"))) {
        cat(paste("Calculating Bootstrap Fits For Topology ", i, "\n"), file = "topology-analysis-log.txt", append = TRUE)
        graph_fit(i)
        cat(paste("Bootstrap Fits For Topology ", i, " Calculated \n"), file = "topology-analysis-log.txt", append = TRUE)
    }
}
stopCluster(cl)

# 7. Compare and save bootstrap analysis fits
if (!file.exists("bootstrap-results")) {
    dir.create("bootstrap-results")
}
core <- detectCores()
cl <- makeCluster(core)
registerDoParallel(cl)
foreach(i = c(1, (focused_scan %>% distinct(topology) %>% filter(topology < max_topologies))$topology + 1), .packages = c("admixtools"), .inorder = TRUE) %dopar% {
    if (!file.exists(paste0("bootstrap-results/top-1_top-", i, ".csv"))) {
        fit_best <- readRDS(file = paste0("graph-fits/fit-top-1.rds"))
        fit_comp <- readRDS(file = paste0("graph-fits/fit-top-", i, ".rds"))
        test <- compare_fits(fit_best[[1]]$score, fit_comp[[1]]$score)
        write.table(test, paste0("bootstrap-results/top-1_top-", i, ".csv"), sep = ",", quote = FALSE, row.names = FALSE)
        cat(paste("Sinificance Testing For Topology ", i, " Saved \n"), file = "topology-analysis-log.txt", append = TRUE)
    }
}
stopCluster(cl)

cat(paste("Topology Analysis Complete -", Sys.time(), "\n"), file = "topology-analysis-log.txt", append = TRUE)
