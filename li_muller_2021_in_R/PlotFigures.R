# Plot figure 3

rm(list = ls())

library(stats)
library(ggplot2)


## Set working directory here
directory_home <- "/Users/eddiewu/Documents/Mon_travail/MY_PHD/Soonwoo/proj_many_controls/econ_replication/li_muller_2021_in_R"
setwd(directory_home)

output_path <- paste0(directory_home, "/my_output")
cv_path <- paste0(directory_home, "/critical_values")


## Get critical values and define getcv() function
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



# Define bhat
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



### Macchiavello settings
iemp <- 1
kmax <- .35
omfac <- 4
fzind <- 4



### Define variables from emptab
file <- file.path(output_path, "emptab.csv")
emptab <- read.csv(file, header = F)

bs <- emptab[6, iemp]
bl <- emptab[5, iemp]
Om11 <- emptab[10, iemp]
Om12 <- emptab[11, iemp]
Om22 <- emptab[12, iemp]
compute_Kdel <- function(kappa) {
  return( kappa * sqrt(emptab[1, iemp]) * sqrt(emptab[4, iemp] / emptab[9, iemp]) )
}
# Kdel <- compute_Kdel(kappa)
Omr <- Om12 / sqrt(Om11 * Om22)
omsd <- sqrt(Om11 * Om22 - Om12^2)
compute_kg <- function(kappa) {
  return(compute_Kdel(kappa) * sqrt(Om11) / omsd)
}
# kg <- compute_kg(kappa)
wb <- (Om11 - Om12)/omsd;

y1 <- bl / sqrt(Om11)
y2 <- (wb / sqrt(Om11) - sqrt(Om11) / omsd) * bl + sqrt(Om11) / omsd * bs
ysign = 1
Omi = solve(rbind(c(Om11, Om12), c(Om12, Om22)))



### Compute data for plotting
pdata <- c()

for (kappa in seq(0, kmax, kmax/50)) {
  if (kappa == 0) {
    kappa <- 0.0000001
  }
  cvx <- numeric(3)
  kg <- compute_kg(kappa)

  for (j in 1:3) {
    cvx[j] <- getcv(j, kg, wb)
  }

  roots <- numeric()
  for (j in 1:3) {
    for (ix in seq(-1, 1, 2)) {
      equation <- function(b) {
        return( getLRr(y1 - b, y2 - wb * b, kg, wb) - cvx[j] )
      }
      b_init <- bhat(y1, y2, kg, wb) + 2 * ix
      print(b_init)
      root <- uniroot(equation, c(b_init-2, b_init+2))$root
      roots <- c(roots, root)
    }
  }

  pdata <- rbind(pdata, c(kappa, roots))
}



### Put into a list of matrices
xx <- list()

for (i in 2:ncol(pdata)) {
  temp <- matrix(nrow = nrow(pdata), ncol = 2)
  for (j in 1:nrow(pdata)) {
    temp[j, ] <- c(pdata[j, 1], sqrt(Om11) * pdata[j, i])
  }
  xx[[i - 1]] <- temp
}

temp <- matrix(nrow = nrow(pdata), ncol = 2)
for (j in 1:nrow(pdata)) {
  temp[j, ] <- c(pdata[j, 1], 0.5 * sqrt(Om11) * (pdata[j, 4] + pdata[j, 5]))
}
xx <- append(xx, list(temp), after = 0)



### Interpolation
temp <- xx[[fzind]]
x <- temp[, 1]
y <- temp[, 2]
f <- splinefun(x, y, method = "monoH.FC")

kstar <- uniroot(f, c(0, 0.5))$root

temp <- list(
  matrix( c(kstar, bl - 0.8 * omfac * sqrt(Om11), kstar, bl + 0.98 * omfac * sqrt(Om11)),
          nrow = 2,
          byrow = TRUE))
xx <- append(xx, temp)



### Final graph plot
plot <- ggplot() + 
  geom_line(aes(x = xx[[1]][,1], y = xx[[1]][,2]), linewidth = 0.7) +
  geom_line(aes(x = xx[[2]][,1], y = xx[[2]][,2]), linetype = 3) +
  geom_line(aes(x = xx[[3]][,1], y = xx[[3]][,2]), linetype = 3) +
  geom_line(aes(x = xx[[4]][,1], y = xx[[4]][,2])) +
  geom_line(aes(x = xx[[5]][,1], y = xx[[5]][,2])) +
  geom_line(aes(x = xx[[6]][,1], y = xx[[6]][,2]), linetype = 4) +
  geom_line(aes(x = xx[[7]][,1], y = xx[[7]][,2]), linetype = 4) +
  geom_line(aes(x = xx[[8]][,1], y = xx[[8]][,2]), linetype = 2) +
  geom_text(aes(x = xx[[8]][1,1], label="k*_{LR}", y = xx[[8]][1,2]), colour="black") +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  xlab(NULL) +
  ylab(NULL) + 
  theme_minimal()

plot

ggsave(file.path(output_path, "macchiavello.pdf"), plot)






