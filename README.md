## Frequency Band Analysis for Multiple Stationary Time Series (FBAM)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13397030.svg)](https://doi.org/10.5281/zenodo.13397030)

This repository contains code that implements the genetic algorithm for 
optimization of the FBAM objective function and for replicating the simulation
studies and application to gait variability data described in the article and
supplementary material. 

This repository is structured as follows:

### FBAM optimization and helper functions

* The directory `R/` contains the code used to implement the optimization routine
described in the article based on a genetic algorithm. Included also are supplementary
code used in the simulation studies.
  * `R/fbam.R` contains all functions needed for optimization of the FBAM objective. The
  main function that should be used is the first function named `fbam()`.
  * `R/generate_sim_data.R` contains functions needed for generations of the simulation
  models described in the article. The main function used here is 
  `generate_simulation_data()`.

### Simulation studies
* The script `sim_study.R` can be used to replicate the simulation studies presented in 
the article. It accepts command line arguments described in the file header. It is recommended
to use this script by calling it from the command line with `Rscript`.
* The script `sim_report.R` generates the `.csv` file with the study results like those presented in the article. It accepts one command line argument: the directory where the output `.rda` files from the study run are stored.
* The directory `BMARD/` contains the code needed to replicate our comparison study
with the [BMARD](https://www.sciencedirect.com/science/article/abs/pii/S0167947321002437)
method. The files in that directory are as follows:
  * The directory `BMARD/` is a clone of the repository provided by the BMARD authors [here](https://github.com/Cuauhtemoctzin/BMARD).
  * `BMARD_comparison.R` is structured similarly to `sim_study.R` and accepts command line arguments. Again, it is recommended to run this script from the command line with `Rscript`.
  * `BMARD_report.R` generates the `.csv` file with the study results like those presented in the supplementary material. It accepts one command line argument: the directory where the output `.rda` files from the study run are stored.
  
### Gait analysis
* The script `gait_process.R` will download the gait variability study data
directly from [physionet.org](https://physionet.org/content/gaitndd/1.0.0/) and process it according to the method described in the article. Raw and processed data will be placed into
corresponding subdirectories in a newly created `data/` directory. Some figures,
including a reproduction of Figure 1 in the article, will be placed into the `figures/` directory (this will be created on the first run).
* The script `gait_analysis.R` will replicate all figures and results presented in
the application section of the article. Before running this file, ensure the
`gait_process.R` script has been successfully sourced.

*Note that the simulation studies were originally performed on high performance
computing clusters. If running these on a local machine with significantly fewer
resources, it is recommended to reduce the number of simulation replications.*

When running the scripts detailed above, the working directory **must**
be set to the root directory of this repository. 
The easiest way of doing this is by opening the `.Rproj` file
located within the root directory in [RStudio](https://posit.co/products/open-source/rstudio/).

## Dependencies

The `doParallel` package is required. It can be installed with

```
install.packages('doParallel')
```

## Processed gait data

Once the `gait.R` script has successfully executed, two files can be found in
`data/gait-processed/':

* `gait.rda` is a R data file which contains a list with the following components:
  * `x`: processed time series data in long format
  * `mtfreq`: Fourier frequencies used in estimation
  * `mtspec`: multitaper estimated relative power spectra in long format
  * `labels`: the group membership of each subject (e.g., "als", "control")
  * `patient_id`: the patient id found in each of the original files (e.g., "als1", "control1")
  * `noclust`: the output of `fbam` under setting (a)
  * `clust`: the output of `fbam` under setting (b)
* `gait.csv`: this is a data frame constructed from the above data and the 
subject level descriptions provided in the original data. It has the following
columns:
  * `patient_id`: same as above
  * `group`: same as `labels` above
  * `age`, `height`, `weight`: age (years), height (meters), and weight (kilograms)
  * `gender`: sex
  * `gait_speed`: average gait speed calculated over entire signal (m/s)
  * `duration_severity`: original duration/severity score (see [here](https://physionet.org/content/gaitndd/1.0.0/)) for more information
  * `duration_severity_std`: standardized duration/severity score to lie between `0` and `1`. A value of `0` means healthy control and `1` corresponds to worst case within the condition. See the `gait.R` script for more information.
  * `noclust_LF`, `noclust_HF`, `noclust_ratio`: low frequency average power, high frequency average power, and the ratio of low to high frequency power computed using the solution under setting (a).
  * `clust_LF`, `clust_HF`, `clust_ratio`: low frequency average power, high frequency average power, and the ratio of low to high frequency power computed using the solution under setting (b).
  * `clust_label`: cluster label under setting (b).
  * `noclust_endpoint`: the endpoint determined by FBAM under setting (a). This is the same for all subjects.
  * `clust_endpoint`: the endpoint determined by FBAM under setting (b). This is not the same for all subjects.
  * `ssp`: self-similarity parameter computed using `DFA`.

### Author Information

Connor K. Brubaker\
Ph.D. Candidate, Department of Statistics\
Texas A&M University\
College Station, Texas, United States\
Email: [brubaker@tamu.edu](mailto:brubaker@tamu.edu)

