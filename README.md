## Acoustic spatial capture-recapture study on vocalising bowhead whales

Contains all code used to do analyse the case study data, run the simulation study and visualise the results. 

The main scripts of interest are:

- simulate_1000_data_sets.R

	Simulates data that can be used for simulation studies. It uses simulate_data_Rcpp.R, and the data was used for the simulation study.
	
- ascrRcpp_build_script.R
	
	This script creates the package ascrRcpp, which is an custom R package with compiled Rcpp code to fit the ASCR models.
	
- varying_fits_to_real_data.R
	
	This script fits 35 models

- bootstrap_real_data.R

	This scripts creates 999 bootstrapped data sets, and fits the best model to all.
	
- create_grid_Albers.R
	
	Create a spatial grid that is used for the fitting, based on the Albers projection.
	
- fitting_to_simulated_data.R
	
	Fits models to simulated data for the simlation study.
	
- hidden_functions.R
	
	Contains some hidden functions required for other scripts.
	
- nll.R
	
	Derives the negative log likelihood using ascrRcpp.
	

### Other scripts

- main_Rcpp.R
	
	Fits an ASCR model to data using the ascrRcpp package.
	
- process_data_august 2020.R 
	
	Processes the raw data.
	

### Data files

- fits_1_35_nlminb_n=443.Rdata
	
	An Rdata file containing the results of 35 ASCR models with varying density specifications.

- detections_31-08-2010_successful_loc.csv

	A matrix of detections for every included call for every detector that was involved with the detection, where TRUE indicates a positive detection;

- received_levels_31-08-2010_successful_loc.csv

	A matrix of received sound levels for every included call for every detector that was involved with the detection, and NA otherwise;

- bearings_31-08-2010_successful_loc.csv

	A matrix of bearings for every included call for every detector that was involved with the detection, and -1 otherwise;

- DASARs.txt

	Some information on the DASARs.
