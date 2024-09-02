#
# Notes on get_li_muller_CI():
# 
# 1. The function preprocess_data() will need to be modified according to your dataset
#    because the data cleaning step is dataset specific. The function should load the raw data and
#    return a clean dataset that contains only [Y, X, Q1, Q2, ... Qm, Z1, Z2, ... Zn, clustervar] as its variables.
#    Alternatively, you can choose to clean your data in advance, let data_cleaning_fun be NULL which is the default,
#    and feed your clean dataset directly into get_li_muller_CI().
# 
# 2. You will need to set up a folder that contains the files with all the critical values,
#    i.e. the same thing as my critical_values folder.
# 
# 3. You will need utils.R in the same folder as this script because 
#    it sources the functions in utils.R.
#
# 4. Li and Muller originally use Mathematica to find the root of a function. 
#    I believe the Mathematica algorithm is a variant of the Newton method.
#    R does not have very good algorithms for finding roots. I believe the uniroot() function of R uses a 
#    bisecting algorithm. The sufficient conditions for the existence of a root is that the function must be continuous,
#    and the search interval [a ,b] must be such that f(a) and f(b) have opposite signs. 
#    In my code, I specify the interval as a neighborhood around b_init, which is
#    the starting point of the Newton method in Li and Muller. 
#    It works fine for Li and Muller's applications and my test cases. 
#    However, please feel free to adapt the uniroot() search interval to 
#    your application. To change the interval, ctrl/cmd + f and search for 
#    uniroot. There is only one place that requires change.
# 
# 
# Arguments:
# data_cleaning_fun - your function for cleaning raw data. It should return a 
#                     data frame with only columns [Y, X, Q1, Q2, ... Qm, Z1, Z2, ... Zn, clustervar].
#                     If your data is already cleaned, no need to specify this argument.
# 
# data - your input data.
# 
# kappa - the kappa value.
# 
# alpha - the significance level. Must be one of [0.01, 0.05, 0.1].
# 
# output_path - your path for saving temporary excel outputs. Can be any path you like.
# 
# cv_path - the path to your folder containing the critical value files.
# 

get_li_muller_CI <- function(data_cleaning_fun = NULL,
                             data,
                             kappa,
                             alpha,
                             output_path,
                             cv_path) {
  
  # Load the required packages and functions
  library(dplyr)
  library(haven)
  library(openxlsx)
  library(estimatr)
  library(R.matlab)
  library(stats)
  library(ggplot2)
  source("utils.R")
  
  # Clean data if you have written your data cleaning steps in preprocess_data()
  if (!is.null(data_cleaning_fun)) {
    df <- data_cleaning_fun(data)
  }
  
  # Define vars to be used in regression formula
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
                               " - 1"))
  model <- lm_robust(formula,
                     se_type = "stata",
                     clusters = clustervar,
                     data = df)
  out <- summary(model)
  na_idx <- which(is.na(out$coefficients[, "Estimate"]))
  df <- df[, !names(df) %in% names(na_idx)]
  
  # Drop collinear Z variables
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
  out <- summary(model)
  na_idx <- which(is.na(out$coefficients[, "Estimate"]))
  df <- df[, !names(df) %in% names(na_idx)]
  
  # Re-define vars since some collinear vars have been dropped
  Y <- "Y"
  X <- "X"
  Q <- names(df)[grep("Q.*", names(df))]
  Z <- names(df)[grep("Z.*", names(df))]
  cl <- "clustervar"
  
  # Re-define the formula for regression
  formula <- as.formula(paste0(Y,
                               " ~ ",
                               X, " + ",
                               paste(Q, collapse = " + "), " + ",
                               paste(Z, collapse = " + "),
                               " - 1"))
  
  # Run the long regression
  model <- lm_robust(formula,
                     se_type = "stata",
                     clusters = clustervar,
                     data = df)
  coefL <- model$coefficients["X"]
  tstatL <- model$coefficients["X"] / model$std.error["X"]
  
  # Run the short regression
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
  
  # Keep the observations that were used in the previous regression
  df <- df %>% filter(!is.na(model$fitted.values))
  df <- df[ , c(Y, X, Q, Z, cl)]
  
  
  ### Compute statistics
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
  
  # Get the number of variables and obs
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
  
  NQ <- length( names(df)[grep("Q.*", names(df))] )
  NZ <- length( names(df)[grep("Z.*", names(df))] )
  N <- model$nobs

  # Save the dataframe
  wb <- createWorkbook()
  addWorksheet(wb, "Design")
  writeData(wb, "Design", df)
  
  # Save the statistics
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
  
  # Save xlsx
  saveWorkbook(wb, file = paste0(output_path, "/Bivariatedata.xlsx"), overwrite = TRUE)
  
  
  ########################################
  
  location_CV <- file.path(cv_path, "CV_vals.mat")
  CV <- readMat(location_CV)[["CV"]]
  
  file_name <- "Bivariatedata.xlsx"
  data <- file.path(output_path, file_name)
  
  # Load the Design sheet
  ds <- read.xlsx(data, sheet = "Design")
  
  # Define variables
  varnames <- names(ds)
  Y <- as.numeric(ds$Y)
  X <- as.numeric(ds$X)
  
  Q_ind <- grepl("Q", varnames)
  Q <- as.matrix(ds[, Q_ind])
  KQ <- ncol(Q)
  Z_ind <- grepl("Z", varnames)
  Z <- as.matrix(ds[, Z_ind])
  XQ <- cbind(X, Q)
  
  lambda <- solve(t(XQ) %*% XQ) %*% (t(XQ) %*% Z)
  lambda <- lambda[1, ]
  
  # Do matrix computations
  MqY <- Y - Q %*% solve(t(Q) %*% Q) %*% (t(Q) %*% Y)
  MqX <- X - Q %*% solve(t(Q) %*% Q) %*% (t(Q) %*% X)
  MqZ <- Z - Q %*% solve(t(Q) %*% Q) %*% (t(Q) %*% Z)
  MxqY <- Y - XQ %*% solve(t(XQ) %*% XQ) %*% (t(XQ) %*% Y)
  N <- length(Y)
  
  
  ### Errors from short regression
  # Clustered standard errors
  cluster_groups <- as.numeric(ds$clustervar)
  cluster_groups_unique <- unique(cluster_groups)
  P <- length(cluster_groups_unique)
  np <- vector("list", P)
  for (ii in 1:P) {
    np[[ii]] <- which(cluster_groups == cluster_groups_unique[ii])
  }
  
  #
  ehat_short <- Y - XQ %*% solve(t(XQ) %*% XQ) %*% (t(XQ) %*% Y)
  sighat_e <- matrix(0, nrow = N, ncol = N)
  for (ii in 1:P) {
    tempind <- np[[ii]]
    temp <- ehat_short[tempind] %*% t(ehat_short[tempind])
    sighat_e[tempind, tempind] <- (temp + t(temp)) / 2
  }
  
  #
  QZ <- cbind(Q, Z)
  XQZ <- cbind(X, Q, Z)
  DLong <- solve(t(XQZ) %*% XQZ) %*% t(XQZ)
  DLong <- DLong[1, ]
  DShort <- solve(t(XQ) %*% XQ) %*% t(XQ)
  DShort <- DShort[1, ]
  D <- rbind(DLong, DShort)
  Dsig <- D
  VarLS <- Dsig %*% sighat_e %*% t(Dsig)
  sig_biv <- list(sigL2 = VarLS[1, 1], sigLS = VarLS[1, 2], sigS2 = VarLS[2, 2])
  
  # Coefficients
  coef_biv <- list(betaL = as.numeric(DLong %*% Y),
                   betaS = as.numeric(DShort %*% Y))
  
  
  ### Load the Statistic sheet
  Statistics <- read.xlsx(data, sheet = "Statistic")
  
  # Update statistics
  Statistics$Value[Statistics$Statistic == "t-statistic long"] <- coef_biv$betaL / sqrt(sig_biv$sigL2)
  Statistics$Value[Statistics$Statistic == "t-statistic short (Matlab)"] <- coef_biv$betaS / sqrt(sig_biv$sigS2)
  Statistics$Value[Statistics$Statistic == "Correlation between long and short"] <- sig_biv$sigLS / sqrt(sig_biv$sigL2 * sig_biv$sigS2)
  
  # Add new statistics
  R2ZXQ <- (t(MqX) %*% MqZ %*% solve(t(MqZ) %*% MqZ) %*% t(MqZ) %*% MqX) / (t(MqX) %*% MqX)
  R2ZYQ <- (t(MqY) %*% MqZ %*% solve(t(MqZ) %*% MqZ) %*% t(MqZ) %*% MqY) / (t(MqY) %*% MqY)
  
  MxqZ <- Z - XQ %*% solve(t(XQ) %*% XQ) %*% t(XQ) %*% Z
  normMqZ2 <- norm(MqZ %*% solve(t(MxqZ) %*% MxqZ) %*% t(MxqZ) %*% Y, "F")^2
  
  new_stats <- data.frame(
    Statistic = c(
      "(Matlab) Betahat long", "(Matlab) Betahat short", "(Matlab) norm(MqY)^2", "(Matlab) norm(MxqY)^2",
      "(Matlab) norm(MqX)^2", "(Matlab) Var(Betahat long)", "(Matlab) Cov(Betahat long, Betahat short)",
      "(Matlab) Var(Betahat short)", "(Matlab) R-squared of Z on X given Q", "(Matlab) R-squared of Z on Y given Q",
      "(Matlab) norm(MqZ*Long coefficient)^2"
    ),
    Value = c(
      coef_biv$betaL, coef_biv$betaS, norm(MqY, "F")^2, norm(MxqY, "F")^2, norm(MqX, "F")^2,
      sig_biv$sigL2, sig_biv$sigLS, sig_biv$sigS2,
      R2ZXQ,
      R2ZYQ,
      normMqZ2
    )
  )
  
  # Append new statistics
  Statistics <- rbind(Statistics, new_stats)
  
  # Write the updated statistics to Excel
  wb <- loadWorkbook(data)
  addWorksheet(wb, sheetName = "Matlab")
  writeData(wb, sheet = "Matlab", x = Statistics, colNames = TRUE, rowNames = FALSE)
  saveWorkbook(wb, file = data, overwrite = TRUE)
  
  
  ########################################
  
  file <- file.path(output_path, "Bivariatedata.xlsx")
  df <- read.xlsx(file, sheet = "Matlab")
  
  out <- c()
  
  out[1:3] <- df[1:3, "Value"]
  out[4] <- df[df$Statistic == "R-squared of Z on X given Q", "Value"]
  out[5] <- df[df$Statistic == "(Matlab) Betahat long", "Value"]
  out[6] <- df[df$Statistic == "(Matlab) Betahat short", "Value"]
  out[7] <- df[df$Statistic == "(Matlab) norm(MxqY)^2", "Value"]
  out[8] <- df[df$Statistic == "(Matlab) norm(MqY)^2", "Value"]
  out[9] <- df[df$Statistic == "(Matlab) norm(MqX)^2", "Value"]
  out[10] <- df[df$Statistic == "(Matlab) Var(Betahat long)", "Value"]
  out[11] <- df[df$Statistic == "(Matlab) Cov(Betahat long, Betahat short)", "Value"]
  out[12] <- df[df$Statistic == "(Matlab) Var(Betahat short)", "Value"]
  
  out <- data.frame(value = out)
  save_file <- file.path(output_path, "emptab.csv")
  write.table(out, file = save_file, sep = ",", row.names = FALSE, col.names = FALSE)
  
  
  ########################################
  
  # So that the optimizer can work when kappa = 0
  if (kappa == 0) {
    kappa <- 1e-8
  }
  
  ### Get critical values and define getcv() function
  cvtab <- list(read.table(file.path(cv_path, "cvmat_1.txt")),
                read.table(file.path(cv_path, "cvmat_5.txt")),
                read.table(file.path(cv_path, "cvmat_10.txt")))
  
  getcv <- function(j, kg, wb) {
    n <- 40
    x <- c(0.5 * (0.5 * log(kg^2 + wb^2) / log(200.0) + 1), 2 * atan(kg / wb) / pi)
    x <- pmin(pmax(x, 1/n), (n - 1)/n)
    i <- floor(n * x)
    w <- n * x - i
    wc <- 1 - w
    
    out <- wc[1] * wc[2] * cvtab[[j]][i[1], i[2]] +
      w[1] * wc[2] * cvtab[[j]][i[1] + 1, i[2]] +
      wc[1] * w[2] * cvtab[[j]][i[1], i[2] + 1] +
      w[1] * w[2] * cvtab[[j]][i[1] + 1, i[2] + 1]
    
    return(out)
  }
  
  getcv(2, 1.3, 3.0)
  
  
  ### Define getLRr() function
  minimize_x0 <- function(y1, y2, kg) {
    objective <- function(g) {
      y1^2 + (y2 - g)^2
    }
    
    out <- optimize(objective, c(-kg, kg))
    
    return(out$objective)
  }
  
  minimize_x1 <- function(y1, y2, kg, wb) {
    objective <- function(param) {
      g <- param[1]
      b <- param[2]
      (y1 - b)^2 + (y2 - wb * b - g)^2
    }
    
    res <- optim(c(0, 0), objective, method = "L-BFGS-B", lower = c(-kg, -Inf), upper = c(kg, Inf))
    
    return(res$value)
  }
  
  getLRr <- function(y1, y2, kg, wb) {
    x0 <- minimize_x0(y1, y2, kg)
    x1 <- minimize_x1(y1, y2, kg, wb)
    
    return(x0 - x1)
  }
  
  
  ### Define bhat
  bhat <- function(y1, y2, kg, wb) {
    if (
      ( y1 < 0 && (
        ( kg + wb * y1 >= y2 && ( (wb * y1 < y2 && y2 < 0) || y2 > 0 ) ) ||
        ( wb * y1 > y2 && (kg + y2 > wb * y1) && y2 < 0 )
      )
      ) ||
      ( y1 > 0 && (
        ( y2 > 0 && (
          ( kg + wb * y1 >= y2 && wb * y1 < y2 ) || ( wb * y1 > y2 && kg + y2 > wb * y1 )
        )
        ) ||
        (y2 < 0 && kg + y2 > wb * y1)
      )
      )
    ) {
      return(y1)
    }
    else if (kg + wb * y1 < y2) {
      return(((-kg * wb) + y1 + (wb * y2)) / (1 + wb^2))
    }
    else if (kg + y2 <= wb * y1) {
      return((y1 + wb * (kg + y2)) / (1 + wb^2))
    }
    else {
      return(NA)
    }
  }
  
  omfac <- 4
  
  ### Define variables using emptab
  file <- file.path(output_path, "emptab.csv")
  emptab <- read.csv(file, header = F)
  
  bs <- emptab[6, 1]
  bl <- emptab[5, 1]
  Om11 <- emptab[10, 1]
  Om12 <- emptab[11, 1]
  Om22 <- emptab[12, 1]
  compute_Kdel <- function(kappa) {
    return( kappa * sqrt(emptab[1, 1]) * sqrt(emptab[4, 1] / emptab[9, 1]) )
  }
  Omr <- Om12 / sqrt(Om11 * Om22)
  omsd <- sqrt(Om11 * Om22 - Om12^2)
  compute_kg <- function(kappa) {
    return(compute_Kdel(kappa) * sqrt(Om11) / omsd)
  }
  wb <- (Om11 - Om12)/omsd;
  
  y1 <- bl / sqrt(Om11)
  y2 <- (wb / sqrt(Om11) - sqrt(Om11) / omsd) * bl + sqrt(Om11) / omsd * bs
  ysign = 1
  Omi = solve(rbind(c(Om11, Om12), c(Om12, Om22)))
  
  # Determine which critical value table to use depending on alpha
  if (alpha == 0.01) {
    j <- 1
  } else if (alpha == 0.05) {
    j <- 2
  } else if (alpha == 0.10) {
    j <- 3
  } else {
    stop("alpha must be a float and one of 0.01, 0.05 or 0.1.")
  }
  
  kg <- compute_kg(kappa)
  cvx <- getcv(j, kg, wb)
  
  ### Compute confidence interval
  roots <- numeric()
  for (ix in seq(-1, 1, 2)) {
    equation <- function(b) {
      return( getLRr(y1 - b, y2 - wb * b, kg, wb) - cvx )
    }
    b_init <- bhat(y1, y2, kg, wb) + 2 * ix
    root <- uniroot(equation, c(b_init-2, b_init+2))$root
    roots <- c(roots, root)
  }
  
  roots <- sqrt(Om11) * roots

  return(roots)
}



preprocess_data <- function(df) {
  ### This function is specific to your raw data.
  ### It takes your raw data as input and returns a dataframe with only
  ### [Y, X, Q1, Q2, ... Qm, Z1, Z2, ... Zn, clustervar] as the column names.
  
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
  names(df)[names(df) == Y] <- "Y"
  names(df)[names(df) == X] <- "X"
  names(df)[names(df) == cl] <- "clustervar"
  
  for (qq in 1:length(Q)) {
    names(df)[names(df) == Q[qq]] <- paste0("Q", qq)
  }
  
  for (zz in 1:length(Z)) {
    names(df)[names(df) == Z[zz]] <- paste0("Z", zz)
  }
  
  return(df)
}




