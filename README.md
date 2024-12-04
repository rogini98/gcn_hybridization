# Effects of hybridization and gene flow on gene co-expression networks

This repository contains R scripts for simulating and analyzing the evolution of gene co-expression networks under two scenarios: immediate hybridization and continuous gene flow. The simulations track how network properties change across multiple generations and different migration rates.

## Overview

The codebase simulates two evolutionary scenarios and tracks three key network metrics:
- Mean Degree (network connectivity)
- Modularity (community structure)
- Transitivity (clustering)

The simulation framework includes:
- Hybridization simulation with parallel replicates
- Introgression tracking with variable migration rates
- Visualization of network evolution and allele frequencies
- Statistical analysis and parental comparisons
- Network metric tracking across generations

## Requirements

### R Packages
```R
required_packages <- c(
  "igraph",      # network analysis
  "ggplot2",     # plotting
  "dplyr",       # data manipulation
  "tidyr",       # data reshaping
  "patchwork",   # combining plots
  "parallel",    # parallel processing
  "doParallel",  # parallel backend
  "foreach",     # parallel iteration
  "pbapply",     # progress bars
  "Matrix",      # sparse matrices
  "pheatmap",    # heatmap visualization
  "RColorBrewer" # color palettes
)
```

## File Structure

```
├── 00-run_all_simulations.R          # Master script to run both simulations
├── 01-run_hybridization_simulation.R  # Runner script for hybridization
├── 02-run_introgression_simulation.R  # Runner script for introgression
│
├── sim_hybridization/               
│   ├── 01-network_functions_hybridizationSim.R   # Core functions
│   ├── 02-run_simulation.R                       # Simulation implementation
│   └── 03-generate_plots.R                       # Plot generation
│
├── sim_introgression/
│   ├── 01-network_functions_introgressionSim.R   # Core functions
│   ├── 02-visualization_functions.R              # Visualization functions
│   ├── 03-simulation_runners.R                   # Simulation implementation
│   └── 04-network_analysis.R                     # Network analysis
│
├── output-hybridization/
│   ├── plots/
│   │   ├── heatmaps/
│   │   └── networks/
│   ├── analysis_log.txt
│   └── simulation_results.RData
│
├── output-introgression/
│   ├── network_analysis/
│   │   ├── averaged_heatmaps/
│   │   ├── averaged_networks/
│   │   └── metrics/
│   └── tracking/
│
└── simulation_results/
    ├── complete_simulation_log.txt
    └── combined_summary.txt
```

## Scripts

### Main Scripts
**00-run_all_simulations.R**
- Master script to run both simulations sequentially
- Handles logging and summary generation
- Creates combined output reports

**01-run_hybridization_simulation.R** & **02-run_introgression_simulation.R**
- Individual runners for each simulation type
- Configure simulation parameters
- Manage output generation

### Hybridization Scripts (in sim_hybridization/)
**01-network_functions_hybridizationSim.R**
- Core functions for network analysis
- Parent population creation
- Genetic recombination
- Mating simulation
- Transcriptome generation

### Introgression Scripts (in sim_introgression/)
**01-network_functions_introgressionSim.R**
- Population genetics
- Migration handling
- Network evolution
- Allele frequency tracking

## Usage

1. Clone the repository:
```bash
git clone https://github.com/rogini98/gcn_hybridization.git
cd gcn_hybridization
```

2. Run all simulations:
```R
source("00-run_all_simulations.R")
```

Or run individual simulations:
```R
# For hybridization
source("01-run_hybridization_simulation.R")

# For introgression
source("02-run_introgression_simulation.R")
```

## Output

The analysis generates:
### Hybridization Results
- Network evolution plots
- Generation-specific heatmaps
- Network visualizations
- Summary statistics

### Introgression Results
- Migration rate effects
- Averaged network visualizations
- Network metric trajectories
- Combined analysis reports

## Parameters

Key simulation parameters:
```R
# Shared parameters
Chr.lengths <- c(10, 10, 5, 5, 10)  # Chromosome lengths
PopSize <- 500                      # Population size
n_replicates <- 10                  # Number of replicates
threshold <- 0.2                    # Network correlation threshold

# Introgression-specific
migration_rates_network <- exp(seq(-10, 0, by = 0.01))
```

## Citation

If you use this code in your research, please cite:

Runghen Rogini, Bolnick Daniel I. (2024). *Effects of hybridization and gene flow in GCN - Simulation Code* (Version v1.0.3). 
Zenodo. https://doi.org/10.5281/zenodo.14271749
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14271749.svg)](https://doi.org/10.5281/zenodo.14271749)

## Contact

Rogini Runghen: rogini.runghen@gmail.com
Daniel Bolnick: daniel.bolnick@uconn.edu

## Acknowledgments

Both Daniel Bolnick and Rogini Runghen contributed to this code (last updated 28 November 2024).
