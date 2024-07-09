
rm(list = ls())

library(dplyr)
library(haven)
library(openxlsx)
library(car)
library(estimatr)

##
directory_home <- "/Users/eddiewu/Documents/Mon_travail/MY_PHD/Soonwoo/proj_many_controls/econ_replication/li_muller_2021_in_R"
setwd(directory_home)


# Load data
df <- read_dta(paste0(directory_home, "/Pre-regression data.dta"))


# Data cleaning and transformations
df <- df %>% filter(!is.na(maxPA1past))

for (ii in 1:3) {
  df <- df %>% mutate(!!paste0("Iseason_", ii) := 1 * (season == ii))
}

lvl <- sort(unique(df$idrel))
for (ii in lvl) {
  df <- df %>% mutate(!!paste0("Iidrel_", ii) := 1 * (idrel == ii))
}

lvl <- sort(unique(df$firmid))
for (ii in lvl) {
  for (jj in 1:3) {
    df <- df %>% mutate(!!paste0("Ifirmid_", ii, "_x_season_", jj) := 1 * ((firmid == ii) & (season == jj)) )
  }
}

lvl <- sort(unique(df$buyerid))
for (ii in lvl) {
  for (jj in 1:3) {
    df <- df %>% mutate(!!paste0("Ibuyerid_", ii, "_x_season_", jj) := 1 * ((buyerid == ii) & (season == jj)) )
  }
}


# Variable definitions
Y <- "lnS"
X <- "lnpastinter"
Q <- c(names(df)[grep("Iidrel.*", names(df))],
       "maxPA1past",
       names(df)[grep("Iseason_.*", names(df))])
Z <- names(df)[grep("Ibuyerid_.*_x_.*", names(df))]
cl <- "firmid"


# Keep vars
df <- df[ , c(Y, X, Q, Z, cl)]

# Rename vars
df <- df %>%
  rename(Y = Y,
         X = X,
         clustervar = cl)

for (qq in 1:length(Q)) {
  names(df)[names(df) == Q[qq]] <- paste0("Q", qq)
}

for (zz in 1:length(Z)) {
  names(df)[names(df) == Z[zz]] <- paste0("Z", zz)
}


# Define vars
Y <- "Y"
X <- "X"
Q <- names(df)[grep("Q.*", names(df))]
Z <- names(df)[grep("Z.*", names(df))]
cl <- "clustervar"


# Drop collinear Q variables
formula <- as.formula(paste0(Y,
                             " ~ ",
                             X, " + ",
                             paste(Q, collapse = " + "),
                             " - 1"
                             ))
model <- lm_robust(formula,
                   se_type = "stata",
                   clusters = clustervar,
                   data = df)
out <- summary(model)
na_idx <- which(is.na(out$coefficients[, "Estimate"]))
df <- df[, !names(df) %in% names(na_idx)]
# Q <- setdiff(Q, names(na_idx))


# Drop collinear Z variables
formula <- as.formula(paste0(Y,
                             " ~ ",
                             X, " + ",
                             paste(Q, collapse = " + "), " + ",
                             paste(Z, collapse = " + "),
                             " - 1"
                             ))
model <- lm_robust(formula,
                   se_type = "stata",
                   clusters = clustervar,
                   data = df)
out <- summary(model)
na_idx <- which(is.na(out$coefficients[, "Estimate"]))
df <- df[, !names(df) %in% names(na_idx)]


# Re-define vars
Y <- "Y"
X <- "X"
Q <- names(df)[grep("Q.*", names(df))]
Z <- names(df)[grep("Z.*", names(df))]
cl <- "clustervar"


# Re-run the reg just to make sure
formula <- as.formula(paste0(Y,
                             " ~ ",
                             X, " + ",
                             paste(Q, collapse = " + "), " + ",
                             paste(Z, collapse = " + "),
                             " - 1"))
model <- lm_robust(formula,
                   se_type = "stata",
                   clusters = clustervar,
                   data = df)
summary(model)

# drop_list <- c("Z1", "Z16", "Z34", "Z37", "Z40", "Z46", "Z49", "Z64", "Z70", "Z73",
#                "Z85", "Z88", "Z100", "Z102", "Z106", "Z109", "Z112", "Z131")
# df <- df[, !names(df) %in% drop_list]


## Long regression
model <- lm_robust(formula,
                   se_type = "stata",
                   clusters = clustervar,
                   data = df)
coefL <- model$coefficients["X"]
tstatL <- model$coefficients["X"] / model$std.error["X"]
message(paste0("The long reg coef and stat of X are: ", coefL, " ", tstatL))


## Short regression
formula <- as.formula(paste0(Y,
                             " ~ ",
                             X, " + ",
                             paste(Q, collapse = " + "),
                             " - 1"))
model <- lm_robust(formula,
                   se_type = "stata",
                   clusters = clustervar,
                   data = df)
coefS <- model$coefficients["X"]
tstatS <- model$coefficients["X"] / model$std.error["X"]
message(paste0("The short reg coef and stat of X are: ", coefS, " ", tstatS))


## Variables for bivariate test
df <- df %>% filter(!is.na(model$fitted.values))
df <- df[ , c(Y, X, Q, Z, cl)]


## Statistics
# R square Y on X given Q 
formula <- as.formula(paste0(Y,
                             " ~ ",
                             paste(Q, collapse = " + "),
                             " - 1"))
model <- lm(formula, data = df)
MqY <- model$residuals

formula <- as.formula(paste0(X,
                             " ~ ",
                             paste(Q, collapse = " + "),
                             " - 1"))
model <- lm(formula, data = df)
MqX <- model$residuals

model <- lm(MqY ~ MqX - 1)
R2_XYgivenQ <- summary(model)$r.squared


# R square X on Z given Q
reg <- function(z) {
  formula <- as.formula(paste0(z, " ~ ",
                               paste(Q, collapse = " + "),
                               " - 1"))
  model <- lm(formula, data = df)
  out <- model$residuals
  return(out)
}

MqZ <- sapply(Z, reg)
model <- lm(MqX ~ MqZ - 1)
R2_ZXgivenQ <- summary(model)$r.squared



# Export results to Excel
NQ <- length( names(df)[grep("Q.*", names(df))] )
NZ <- length( names(df)[grep("Z.*", names(df))] )

# Get the number of obs
formula <- as.formula(paste0(Y,
                             " ~ ",
                             X, " + ",
                             paste(Q, collapse = " + "), " + ",
                             paste(Z, collapse = " + "),
                             " - 1"))
model <- lm_robust(formula,
                   se_type = "stata",
                   clusters = clustervar,
                   data = df)
N <- model$nobs

temp <- model$coefficients
Betlong <- data.frame(y1 = as.numeric(temp))
rownames(Betlong) <- names(temp)

Covlong <- model$vcov


# Save results
wb <- createWorkbook()

addWorksheet(wb, "Covlong")
writeData(wb, "Covlong", Covlong, rowNames = TRUE)

addWorksheet(wb, "Betlong")
writeData(wb, "Betlong", Betlong, rowNames = TRUE)

addWorksheet(wb, "Design")
writeData(wb, "Design", df)


# Specifications
specifications <- data.frame(
  Specification = c("Paper", "Baseline analysis", "Y variable", "X variable", "Q variable", "Z variable", "Other"),
  Value = c("Macchiavello, Morjaria (AER, September 2015)", "Table 5, Column 8", "log(Value)", 
            "log(# previous interactions)", "Season FEs, relationship FEs, max of past auction prices", 
            "Buyer x season FEs", "SEs clustered by seller; baseline includes Y,X,Q, and Z is added for bivariate test")
)
addWorksheet(wb, "Specification")
writeData(wb, "Specification", specifications)


# Statistics
statistics <- data.frame(
  Statistic = c("Number of observations",
                "Number of baseline controls Q",
                "Number of potential controls Z", 
                "Betahat long",
                "t-statistic long",
                "Betahat short",
                "t-statistic short (Matlab)", 
                "t-statistic short (Stata)",
                "R-squared of X on Y given Q",
                "R-squared of Z on X given Q",  
                "Correlation between long and short"),
  Value = c(N,
            NQ,
            NZ,
            coefL,
            NA,
            coefS,
            NA,
            tstatS,
            R2_XYgivenQ,
            R2_ZXgivenQ,
            NA))
addWorksheet(wb, "Statistic")
writeData(wb, "Statistic", statistics)


# QC
idx <- statistics$Statistic %in% c("Betahat long", "Betahat short", "t-statistic short (Stata)")
qc <- statistics[idx,]
addWorksheet(wb, "QC")
writeData(wb, "QC", qc)


# Save xlsx
saveWorkbook(wb, file = paste0(directory_home, "/Bivariatedata1.xlsx"), overwrite = TRUE)






