# dynamin_model
code and inputs for models of dynamin clustering

#Data folder
Excel sheet containing example parameter sets found for optimal agreement with experiment.
These are a sample. Sorted by isoform and D_{sol} concentration.

#Python folder
Python code that solves the ODE version of the model. There is a separate file for each isoform due to the differences in the geometry and the experimental parameters. 
This code performs the genetic algorithm that finds sets of optimal parameter sets that minimize the chi-squared distance between the experimental growth in density in the clusters, and the model result.
These versions (GA_Opt_2Dens...ipynb) optimize all 8 parameters in Table S2 within specified ranges.

In GA_opt_loop_allIsoforms.ipynb, can switch between isoforms using 'isoform' keyword as 'AA','AB', or 'BB'.
This version uses fixed values of Dsol0 (loops over each value) and Dmem0.

*For Dyn1AB, if Dsol0 is small, and Dmem0 is high, following stimulation, copy numbers on membrane increase 7-fold, so if total copies pre-stimulation are too low, this generates negative copy numbers, so Dmem0 must be restricted to be lower.

#Matlab folder
Matlab code that also solves the ODE version of the model. Generates plots of the solutions.
Experimental growth curve in expdataBB_from4.dat

