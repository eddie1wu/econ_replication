
rm(list = ls())
library("R.matlab")
library("readxl")
library("xlsx")
library("dplyr")

source("utils.R")

##
directory_home <- "/Users/eddiewu/Documents/Mon_travail/MY_PHD/Soonwoo/proj_many_controls/replication_in_R"
setwd(directory_home)


##
location_CV <- "CV_vals.mat"
CV <- readMat(location_CV)[["CV"]]


examples <- "Bivariatedata1.xlsx"
data <- file.path(directory_home, examples)

alpha <- 0.05


##
ds <- read_excel(data, sheet = "Design")


##
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

##

MqY <- Y - Q %*% solve(t(Q) %*% Q) %*% (t(Q) %*% Y)
MqX <- X - Q %*% solve(t(Q) %*% Q) %*% (t(Q) %*% X)
MqZ <- Z - Q %*% solve(t(Q) %*% Q) %*% (t(Q) %*% Z)
MxqY <- Y - XQ %*% solve(t(XQ) %*% XQ) %*% (t(XQ) %*% Y)

N <- length(Y)



### Errors from short regression
## Clustered standard errors
cluster_groups <- as.numeric(ds$clustervar)
cluster_groups_unique <- unique(cluster_groups)
P <- length(cluster_groups_unique)
np <- vector("list", P)
for (ii in 1:P) {
  np[[ii]] <- which(cluster_groups == cluster_groups_unique[ii])
}


##
ehat_short <- Y - XQ %*% solve(t(XQ) %*% XQ) %*% (t(XQ) %*% Y)
sighat_e <- matrix(0, nrow = N, ncol = N)
for (ii in 1:P) {
  tempind <- np[[ii]]
  temp <- ehat_short[tempind] %*% t(ehat_short[tempind])
  sighat_e[tempind, tempind] <- (temp + t(temp)) / 2
}


##
QZ <- cbind(Q, Z)
XQZ <- cbind(X, Q, Z)
DLong <- solve(t(XQZ) %*% XQZ) %*% t(XQZ)
DLong <- DLong[1, ]
DShort <- solve(t(XQ) %*% XQ) %*% t(XQ)
DShort <- DShort[1, ]
D <- rbind(DLong, DShort)
rm(Dsig)
Dsig <- D
VarLS <- Dsig %*% sighat_e %*% t(Dsig)
sig_biv <- list(sigL2 = VarLS[1, 1], sigLS = VarLS[1, 2], sigS2 = VarLS[2, 2])



## Coefficients
coef_biv <- list(betaL = as.numeric(DLong %*% Y),
                 betaS = as.numeric(DShort %*% Y))

# QC
QC <- read_excel(data, sheet = "QC")
QC <- as.numeric(QC$Value)
stopifnot(abs((QC[1] - coef_biv$betaL) / QC[1]) < 10^(-3))
stopifnot(abs((QC[2] - coef_biv$betaS) / QC[2]) < 10^(-3))

# Long t-stat
tstatL <- coef_biv$betaL / sqrt(sig_biv$sigL2)



##
# 
R2_grid <- seq(0, 0.35, by = 0.00001)
Test <- numeric(length(R2_grid))


#
r0 <- norm(MqY, "F") * sqrt( t(lambda) %*% solve(t(MqZ) %*% MqZ) %*% lambda )
ii <- 1*(alpha == 0.01) + 2*(alpha == 0.05) + 3*(alpha == 0.1)
CV_struct <- CV[ , 1, ii]


##
# 
for (iTest in 1:length(R2_grid)) {
  R2_bound <- R2_grid[iTest]
  r_d <- r0 * sqrt(R2_bound)
  Test[iTest] <- logLik_GLR(coef_biv$betaL, coef_biv$betaS, sig_biv, r_d, CV_struct)
}


Test1 <- c(Test[-1], -999)
Testm1 <- c(-999, Test[-length(Test)])
Cutoffkappa <- data.frame(
  Test = Test[(Test != Test1) | (Test != Testm1)],
  R2_grid = R2_grid[(Test != Test1) | (Test != Testm1)]
)
Cutoffkappa$R2_grid <- Cutoffkappa$R2_grid * norm(MqY, "F")^2



# Write the results to Excel
write.xlsx(Cutoffkappa,
           file = data,
           sheetName = "Kappabound",
           col.names = FALSE,
           row.names = FALSE,
           append = TRUE)

idata <- 1
sheet_name <- as.character(idata)
results <- data.frame(
  Kappa2 = Cutoffkappa$R2_grid[2],
  sqrtKappa2_div_N = sqrt(Cutoffkappa$R2_grid[2] / N),
  Kappa2_div_normMxqY2 = Cutoffkappa$R2_grid[2] / (norm(MxqY, "F")^2)
)
write.xlsx(t(results),
           file = file.path(directory_home, 'Empirical.xlsx'),
           sheetName = sheet_name,
           col.names = FALSE,
           row.names = FALSE)




##
Statistics <- read_excel(data, sheet = "Statistic")

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
write.xlsx(Statistics,
           file = data,
           sheetName = "Matlab",
           append = TRUE)




