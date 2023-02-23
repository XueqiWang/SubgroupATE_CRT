####################################################
# Simulation program for the omnibus test
####################################################

source("gendata.R")
source("calcSubgroup.R")

library(nlme)
library(openxlsx)

beta_2 <- 0.2
beta_4 <- 0.1

m <- c(rep(20,9), rep(50,9), rep(100,9))
rho_y <- rep(c(rep(0.02,3), rep(0.05,3), rep(0.1,3)), 3)
rho_s <- rep(c(0.1, 0.25, 0.5), 9)
table <- cbind(m, rho_y, rho_s)

n <- numeric(27)
pred.power <- numeric(27)
for (i in 1:nrow(table)){
  m <- as.numeric(table[i,][1])
  rho_y <- as.numeric(table[i,][2])
  rho_s <- as.numeric(table[i,][3])
  
  np <- calc_omnibus(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=rho_y, rho_s=rho_s)
  n[i] <- np[1]
  pred.power[i] <- np[2]
}

table <- cbind(table, n)



# function to compute empirical power or empirical type I error
empirical_omnibus <- function(nullcase=F, parameter, nsims=2000){
  
  m <- as.numeric(parameter[1])
  rho_y <- as.numeric(parameter[2])
  rho_s <- as.numeric(parameter[3])
  n <- as.numeric(parameter[4])
  
  beta1 <- 0
  beta2 <- beta_2
  beta3 <- 0.15
  beta4 <- beta_4
  
  if (nullcase==T){
    beta2 <- 0
    beta4 <- 0
  }
  
  pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(29+2023*i)
    simdata <- gendata(beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
                       m=m, rho_y=rho_y, rho_s=rho_s, n=n)
    
    fit <- try(lme(Y ~ ctrt + s + ctrt:s, data=simdata, random= ~ 1 | cluster.id), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    R <- matrix(c(0,1,0,0,0,1,0,1), nrow=2, byrow=T)
    beta <- fit$coef$fixed

    test.stat <- as.numeric(t(R%*%beta) %*% solve(R%*%vcov(fit)%*%t(R)) %*% (R%*%beta))/2
    pvalue[i] <- 1-pf(test.stat, 2, n-2)
    
  }
  
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  return(c(empirical, error.rate))
}

# compute empirical power (and error rate)
empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power <- rbind(empirical.power, empirical_omnibus(parameter=table[i,]))
}

# compute empirical type I error (and error rate)
empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe <- rbind(empirical.tIe, empirical_omnibus(nullcase=T, parameter=table[i,]))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
names(result)[5:9] <- c("emp.tIe", "nconv.rate", "emp.power", "nconv.rate.p", "pred.power") 

write.xlsx(result, "conResults/empOmnibus.xlsx", rowNames=F)


