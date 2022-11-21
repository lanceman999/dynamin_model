# dynamin_model
code and inputs for models of dynamin clustering

#Data folder
Excel sheet containing example parameter sets found for optimal agreement with experiment.
These are a sample. Sorted by isoform and D_{sol} concentration.

#Python folder
Python code that solves the ODE version of the model. There is a separate files for each isoform due to the differences in the geometry and the experimental parameters. 
This code performs the genetic algorithm that finds sets of optimal parameter sets that minimize the chi-squared distance between the experimental growth in density in the clusters, and the model result.

#Matlab folder
Matlab code that also solves the ODE version of the model. Generates plots of the solutions.

