## Bowhead Whales

Contains all code used to do analyse the case study data, run the simulation study and visualise the results. 

The main scripts of interest are:

- simulate_main.R
	Simulates data that can be used for simulation studies. It uses simulate_data_Rcpp.R. 
- ascrRcpp_build_script.R
	This script creates the package ascrRcpp, which is an custom R package with compiled Rcpp code to fit the ASCR models.
- varying_fits_to_real_data.R
	This script fits 33 models
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