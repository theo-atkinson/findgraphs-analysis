# Findgraphs Analysis
This repository includes a series of R scripts that can be used to carry out admixture analysis using the findgraphs tool from [ADMIXTOOLS 2](https://uqrmaie1.github.io/admixtools/index.html). 

These scripts follow the workflow for fitting admixture graphs as described in the [eLife paper](https://elifesciences.org/articles/85492) that introduces the tool.

Various parameters must be set at the top of the scripts before running.

See the NUS folder for an example script on how to create the input files, how to run scripts in the NUS HPC system, and for more info on workflow.

<br />

## Setup
The input files for computing f2 stats are created using _convertf_ from [eigensoft](https://github.com/argriffing/eigensoft)

### Computing f2 Stats:
See the [f-statistics page](https://uqrmaie1.github.io/admixtools/articles/fstats.html) on the admixtools site for more info.

Calculation of population level f2 stats example:
```
# Load package
library("admixtools")

# Set working directory
setwd("[DIRECTORY WITH INPUT FILES]")

# Calculate f2 stats
extract_f2("[INPUT FILES NAME]", "[OUTPUT DIRECTORY]", auto_only = FALSE, adjust_pseudohaploid = FALSE, overwrite = TRUE)
```
<br />

## Script 1: Initial Scan
The script "initial-scan.r" will run findgraphs in parallel a set number of times for a various number of admixture events. 
These can be set at the top of the script using the _iterations_ and _max_admix_ variables respectively.
#### Output
initial-scan-result.rds
* A tibble of the best fitting graphs from each iteration.

initial-scan-plot.jpg
* A plot for interpretting the results.
* Each point represents the best fitting graph for each findgraphs iteration.
* The blue line highlights the graph with the lowest worst residual for each number of admixture events.

initial-scan-log.txt
* A log file that will update while the script is running.

Specific [_find_graphs()_](https://uqrmaie1.github.io/admixtools/reference/find_graphs.html)
and [_qpgraph()_](https://uqrmaie1.github.io/admixtools/reference/qpgraph.html) parameters can be edited on lines 36 and 43 respectively.

<br />

## Script 2: Focused Scan
The script "focused-scan.r" works in a similar manner to the initial scan but runs findgraphs at a set number of admixture events using the _num_admix_ variable at the top of the script.

#### Output
focused-scan-result.rds
* A tibble of the best fitting graph from each iteration.

focused-scan-log.txt
* A log file that will update while the script is running.

Specific [_find_graphs()_](https://uqrmaie1.github.io/admixtools/reference/find_graphs.html)
and [_qpgraph()_](https://uqrmaie1.github.io/admixtools/reference/qpgraph.html) parameters can be edited on lines 33 and 37 respectively.

<br />

## Script 3: Topology Analysis
The script "topology-analysis.r" takes the results of the focused scan and outputs various utilites for interpreting the best fitting topologies.
#### Output
focused-scan-constrained.RDS
* A tibble of the best fitting graphs from the focused scan that follow the constraints defined as _admix_constraints_ and _event_constraints_ at the top of the script.
* See the documentation for [_satisfies_numadmix_](https://uqrmaie1.github.io/admixtools/reference/satisfies_numadmix.html) and
[_satisfies_eventorder_](https://uqrmaie1.github.io/admixtools/reference/satisfies_eventorder.html) on how to implement this.
* Note: This file will still be created and be identical to "focused-scan-result.rds" if no constraints are defined.

score-WR-plot.jpg
* A scatter plot of worst residual vs score.
* Each dot represents the best fitting graph for each focused scan iteration.

scores-hist.jpg
* A histogram of the scores from the best fitting graph of each focused scan iteration.

topology-list.csv
* A list of all of the best fitting graphs from each focused scan iteration, sorted by score and labelled with topology.

graph-topology-N
* A directory with a "graph-topology-N(X).pdf" plot of each graph and a "graph-topology-N.csv" file list of all graphs for topology N, sorted by score.

graph-fits
* A directory with "fit-top-N.rds" files containing bootstrap fits calculated for the best graph of each topology.

bootstrap-results
* A directory with "top-1_top-N.csv" files containing the significance test results between topology N against topology 1.

topology-analysis-log.txt
* A log file that will update while the script is running.

This script is designed to be run multiple times to analyse additional topologies if required. As such, if the output file/directory is already present in the 
parent directory, the script will not run this step again. 

For example, the topology analysis may first be run with the "max_tops" parameter at 5 in order to analyse the top 5 topologies, but after viewing the results more topologies must be analysed so the parameter is increased to 10. The script will not run the step to output "focused-scan-constrained.RDS", "score-WR-plot.jpg", "scores-hist.jpg", or "topology-list.csv", and will only output "graph-topology-N", "graph-fits", and "bootstrap-results" files for topologies 6-10.
  
The log file will indicate which steps have been performed.

<br />

## Increase Complexity

The authors suggest that once the best fitting topologie(s) for a given number of admixture events have been identified, this should be tested for robustness
by rerunning findgraphs at a greater complexity level.

This can be performed by rerunning the focused scan and topology analysis having set the number of admixture events to be one greater than what was previously set.
