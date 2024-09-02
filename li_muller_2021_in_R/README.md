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
- [utils.R](utils.R): contains the utility functions for `C2_loop_bivariate.R` and `get_li_muller_CI.R`.
- [get_li_muller_CI.R](get_li_muller_CI.R): contains the function for computing the Li and Muller confidence interval, using raw data as input. Please refer to the notes in the script itself for more details.
- [main.R](main.R): run this script to replicate Li and Muller Figure 3, from raw data to the saved figure.
- [run_demo.R](run_demo.R): this script demonstrates how to use get_li_muller_CI() to obtain the Li and Muller confidence 
interval using raw data as input. It also contains some test cases that I wrote to ensure the function works correctly.

Running order for replicating Figure 3:
1. Run `main.R` directly. 
OR 
2. Run scripts individually in the order described in `main.R`.

To use get_li_muller_CI() in your project:
1. Copy my **critical_values** folder and [utils.R](utils.R) and paste them into your working directory.
2. Refer to the notes in [get_li_muller_CI.R](get_li_muller_CI.R) for more details on the function e.g. 
things that might require slight modification for your application.
3. Refer to [run_demo.R](run_demo.R) for a demonstration on using this function.



