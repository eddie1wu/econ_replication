# Generate emptab.csv
rm(list = ls())

library(openxlsx)

## Set working directory here
directory_home <- "/Users/eddiewu/Documents/Mon_travail/MY_PHD/Soonwoo/proj_many_controls/econ_replication/li_muller_2021_in_R"
setwd(directory_home)

output_path <- paste0(directory_home, "/my_output")

file <- file.path(output_path, "Bivariatedata1.xlsx")
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





