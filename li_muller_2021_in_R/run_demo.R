
### This script demonstrates how the function get_li_muller_CI() can be used

rm(list = ls())

source("get_li_muller_CI.R")


# Set working directory here
directory_home <- "/Users/eddiewu/Documents/Mon_travail/MY_PHD/Soonwoo/proj_many_controls/econ_replication/li_muller_2021_in_R"
setwd(directory_home)

# Specify the file paths 
data_path <- paste0(directory_home, "/data")
output_path <- paste0(directory_home, "/my_output")  # path to save output
cv_path <- paste0(directory_home, "/critical_values")  # path to the folder containing the critical values


# Load raw data
df <- read_dta(paste0(data_path, "/Pre-regression data.dta"))


# Obtain confidence interval
# Please check the function script for more details on the arguments etc
get_li_muller_CI(data_cleaning_fun = preprocess_data,
                 data = df,
                 kappa = 0,
                 alpha = 0.01,
                 output_path = output_path,
                 cv_path = cv_path)



### --- Generate fake data and test that the function works as expected

# Generate data
N <- 1000
NQ <- 20
NZ <- 20
X <- rnorm(N, 1, 1)
Q <- matrix(rnorm(NQ * N, 4, 1), nrow = N)
Z <- matrix(rnorm(NZ * N, 7, 1), nrow = N)
Y <- 2 * X + Q %*% rnorm(NQ, 10, 2) + Z %*% rnorm(NZ, 0, 1)
clustervar <- sample(1:20, N, replace = T)
df <- as.data.frame(cbind(Y, X, Q, Z, clustervar))

# Rename columns
Q_col <- paste0('Q', seq(1, NQ))
Z_col <- paste0('Z', seq(1, NZ))
col_names <- c('Y', 'X', Q_col, Z_col, 'clustervar')
names(df) <- col_names

# Test the function for various values of kappa and alpha
for (k in seq(0, 0.5, 0.02)) {
  for (a in c(0.01, 0.05, 0.1)) {
    
    CI <- get_li_muller_CI(data = df,
                           kappa = k,
                           alpha = a,
                           output_path = output_path,
                           cv_path = cv_path)
    print(CI)
  }
}






