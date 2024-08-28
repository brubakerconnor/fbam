# BMARD

The BMARD method stands for Bayesian Mixture Auto-Regressive Decomposition. BMARD is a stationary time series decomposition into latent processes such that each of the components is specified as a second-order autoregressive process AR(2). You can find the details and justifications of the method in the preprint https://arxiv.org/abs/2102.11971 where the main application is the analysis of brains signals if a rat during a non-spatial memory task. 

## Data

The method's performance can be tested by using this project; to that end, we share two types of data to try the BMARD method.

### SimulationDATABASES

The folder contains the databases of the periodograms computed for each of the simulations settings described in section 3. There are three scenarios:
A mixture of AR(2) processes
An autoregressive process of order 12
A moving average process of order 4

### Test_data 
At this point, we cannot share the original data, so we share PseudoData. In the experiment, rats learned to recognize five smells in a specific order; the pseudodata data is based on a trial where the smell C=Anise is presented in an incorrect order. There are two databases the signals before and after the odor delivery.

## Code
The code structure is described below; for details, consult the Reproduce_instructions.docx guide and BMARD_simuations_data_dictionary.txt file in the Guides folder.

### CPPcode
Main Rcpp code of the BMARD method the primary function to run the MCMC is called Spectral BDP. The folder contains the code to apply the Bernstein polynomial method of Choudhuri (2014). 

### Rcode_auxiliary
Auxiliary R code for MCMC post-processing is on the code ExtractionBMARDmaincomponentsmodes.R, which computes the mean and median curves and the mean of the parameters determining the optimal number of components of the multichain run. The folder also contains the code to reproduce the 1000 simulations of the scenarios mentioned before and the codes to summarize the multichain MCMC generated. 

## Instructions for use
### Guides

The Guides folder contains detailed instructions to reproduce Figure 2, Table 1, and Table 2 from section 3.

### to test a single example of the BMARD method

Run the following code to simulate a mixture of AR(2) processes and fit BMARD and compute a simple graphical summary of the results. This pipeline is replicated for an AR(12) process. 

BMARD_Simulation_Replication.R

To test BMARD on the pseudodata provided, run the following code.

BMARD_DATA_Replication.R
