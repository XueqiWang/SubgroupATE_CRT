#######################################################################################
# Function to calculate the power for the omnibus test
#######################################################################################

power_omnibus <- function(alpha=0.05, sigma2_y=1, beta_2, beta_4, n, m, rho_y, rho_s, p_i=0.5, p_1=0.5){
  delta_0 <- beta_2
  delta_1 <- beta_2 + beta_4
  p_0 <- 1 - p_1
  
  sigma2_hte <- sigma2_y*(1-rho_y)*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m*(1+(m-2)*rho_y-(m-1)*rho_s*rho_y))
  sigma2_ate <- sigma2_y*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m)
  
  omega_d0 <- sigma2_ate + p_1/p_0*sigma2_hte
  omega_d1 <- sigma2_ate + p_0/p_1*sigma2_hte
  omega_d01 <- sigma2_ate - sigma2_hte
  
  delta_01 <- c(delta_0, delta_1)
  omega <- matrix(c(omega_d0, omega_d01, omega_d01, omega_d1), nrow=2, byrow=T)
  theta <- n*t(delta_01) %*% solve(omega) %*% delta_01
  pred.power <- pf(qf(1-alpha, 2, n-2), 2, n-2, ncp = theta, lower.tail = F)
  return(data.frame(rho_y=rho_y, rho_s=rho_s, t_test=pred.power))
}



#######################################################################################
# Function to calculate the power for the intersection-union test
#######################################################################################

power_IU <- function(alpha=0.05, sigma2_y=1, beta_2, beta_4, n, m, rho_y, rho_s, p_i=0.5, p_1=0.5){
  delta_0 <- beta_2
  delta_1 <- beta_2 + beta_4
  p_0 <- 1 - p_1
  
  sigma2_hte <- sigma2_y*(1-rho_y)*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m*(1+(m-2)*rho_y-(m-1)*rho_s*rho_y))
  sigma2_ate <- sigma2_y*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m)
  
  omega_d0 <- sigma2_ate + p_1/p_0*sigma2_hte
  omega_d1 <- sigma2_ate + p_0/p_1*sigma2_hte
  omega_d01 <- sigma2_ate - sigma2_hte
  
  mean_W <- c(delta_0/sqrt(omega_d0/n), delta_1/sqrt(omega_d1/n))
  cov_W <- (omega_d0)^(-1/2) * omega_d01 * (omega_d1)^(-1/2)
  sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)
  pred.power <- pmvt(df=n-2, lower=c(qt(1-alpha, n-2), qt(1-alpha, n-2)), upper=rep(Inf,2), delta=mean_W, sigma=sigma_W)
  return(data.frame(rho_y=rho_y, rho_s=rho_s, t_test=pred.power))
}


