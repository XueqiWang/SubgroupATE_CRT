##################################################
# Read data
# Size and Power
# Comparison with a back-of-the-envelope approach
##################################################

source("calcSubgroup.R")
source("powerSubgroup.R")

library(openxlsx)
library(mvtnorm)

############################
# Omnibus test
############################

result <- read.xlsx("conResults/empOmnibus.xlsx", colNames=TRUE)
result$power.diff <- result[, 7] - result[, 9]
colnames(result)[c(5:8)] <- c("size", "size.error_rate", "power", "power.error_rate")
result$scenario <- 1:nrow(result)

beta_2 <- 0.2
beta_4 <- 0.1

m <- c(rep(20,9), rep(50,9), rep(100,9))
rho_y <- rep(c(rep(0.02,3), rep(0.05,3), rep(0.1,3)), 3)
rho_s <- rep(c(0.1, 0.25, 0.5), 9)
table <- cbind(m, rho_y, rho_s)

nc <- numeric(27)
actual.power <- numeric(27)
for (i in 1:nrow(table)){
  m <- as.numeric(table[i,][1])
  rho_y <- as.numeric(table[i,][2])
  rho_s <- as.numeric(table[i,][3])
  
  np <- calc_omnibus(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=0, rho_s=0)
  n_c <- np[1]*(1+(m-1)*rho_y)
  n_c <- ceiling(n_c)
  if (n_c%%2 != 0){
    n_c <- n_c + 1
  }
  nc[i] <- n_c
  
  actual.power[i] <- as.numeric(power_omnibus(beta_2=beta_2, beta_4=beta_4, n=n_c, m=m, rho_y=rho_y, rho_s=rho_s))[3]
}

result$nc <- nc
result$actual.power <- actual.power
write.xlsx(result, file = "final_omnibus.xlsx", rowNames = FALSE)



############################
# Intersection-union test
############################

result <- read.xlsx("conResults/empIU.xlsx", colNames=TRUE)
result$power.diff <- result[, 7] - result[, 9]
colnames(result)[c(5:8)] <- c("size", "size.error_rate", "power", "power.error_rate")
result$scenario <- 1:nrow(result)

beta_2 <- 0.3
beta_4 <- 0.1

m <- c(rep(20,9), rep(50,9), rep(100,9))
rho_y <- rep(c(rep(0.02,3), rep(0.05,3), rep(0.1,3)), 3)
rho_s <- rep(c(0.1, 0.25, 0.5), 9)
table <- cbind(m, rho_y, rho_s)

nc <- numeric(27)
actual.power <- numeric(27)
for (i in 1:nrow(table)){
  m <- as.numeric(table[i,][1])
  rho_y <- as.numeric(table[i,][2])
  rho_s <- as.numeric(table[i,][3])
  
  np <- calc_IU(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=0, rho_s=0)
  n_c <- np[1]*(1+(m-1)*rho_y)
  n_c <- ceiling(n_c)
  if (n_c%%2 != 0){
    n_c <- n_c + 1
  }
  nc[i] <- n_c
  
  actual.power[i] <- as.numeric(power_IU(beta_2=beta_2, beta_4=beta_4, n=n_c, m=m, rho_y=rho_y, rho_s=rho_s))[3]
}

result$nc <- nc
result$actual.power <- actual.power
write.xlsx(result, file = "final_IU.xlsx", rowNames = FALSE)


