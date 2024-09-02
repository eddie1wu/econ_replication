
### This script executes the whole process of replicating Figure 3

# One can safely ignore the warning "18 coefficients  not defined because the design matrix is rank deficient"
# Because it arises from C1_Data_preparation.R which removes perfectly collinear regressors

source("C1_Data_preparation.R")

source("C2_loop_bivariate.R")

source("gen_emptab.R")

source("PlotFigures.R")









