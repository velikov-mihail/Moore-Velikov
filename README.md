# Moore--Velikov
 Code used to create results in Moore and Velikov (WP, 2022), Oil price exposure and the cross section of stock returns

This repository contains code used to create the results in Moore and Velikov (WP, 2022), Oil price exposure and the cross section of stock returns. This code is to be used in conjunction with the MATLAB asset pricing package that accompanies Novy-Marx and Velikov (WP, 2022), Assaying Anomalies. 

The order of operations to replicate the results in Moore and Velikov (WP, 2022) is:

1. Download and follow the instructions for setting up the MATLAB Toolkit from https://github.com/velikov-mihail/AssayingAnomalies.git
	* The results in Moore and Velikov (WP, 2022) use the pre-release v0.3 of the MATLAB Toolkit.
3. Download the code in this repository.
4. Run mv.m. The script requires setting up the directories for the MATLAB asset pricing package repository and this repository. It starts a log file and calls multiple other scripts which perform the following functions:  
	* make_daily_betas.m creates and stores the daily betas to be used to calculate CARs for the construction of the weighted CAR3's
	* make_weighting_function.m creates and stores the weighted CAR3s
	* make_oil_prices.m download and stores the oil price data from FRED
	* make_oil_response_forecasts.m creates and stores the oil response forecast variables
	* make_ad_hoc_data_results.m creates and stores severeal auxiliary files
	* make_tables.m prints all tables in ~/Results/MooreVelikovTablesOutput.txt
	* make_figures.m creates and stores all figures in ~/Figures/
   
