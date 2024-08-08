# econ_replication

This folder contains the R code for replicating Figure 3 of Li and Muller 2021. The folder structure is as follows:

- **authors_output** contains the csv and xlsx outputs from the authors' original replication code. I use them for comparison to my own output to check for code correctness.

- **critical_values** contains tables of critical values which are loaded in R and used for calculations.

- **data** contains the pre-regression data in .dta format.

- **my_output** contains outputs from my own R code. The replicated Figure 3 is also in this folder.

Explanation of R scripts:
- [C1_Data_preparation.R](C1_Data_preparation.R): loads `Pre-regression data.dta`, runs the various regressions, and outputs the first few sheets of `Bivariatedata1.xlsx`.
- [C2_loop_bivariate.R](C2_loop_bivariate.R): loads `CV_vals.mat` and `Bivariatedata1.xlsx`, runs matrix computations, and outputs `Empirical.xlsx` and the last few sheets of `Bivariatedata1.xlsx`.
- [gen_emptab.R](gen_emptab.R): loads `Bivariatedata1.xlsx` and outputs `emptab.csv`.
- [PlotFigures.R](PlotFigures.R): loads the three .txt critical value tables and `emptab.csv`, run computations, and plots Figure 3 which is saved as `macchiavello.pdf`.
- [utils.R](utils.R): contains the utility functions for `C2_loop_bivariate.R`.

Running order:
1. Run `main.R` directly. 
OR
1. Run scripts in the order described in `main.R`.



