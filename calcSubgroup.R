################################################################################################
# Function to calculate the required number of clusters for the omnibus test
################################################################################################

calc_omnibus <- function(alpha=0.05, beta=0.2, sigma2_y=1, beta_2, beta_4, m, rho_y, rho_s, p_i=0.5, p_1=0.5){
  delta_0 <- beta_2
  delta_1 <- beta_2 + beta_4
  p_0 <- 1 - p_1
  
  sigma2_hte <- sigma2_y*(1-rho_y)*(1+(m-1)*rho_y) / (p_i*(1-p_i)*p_1*p_0*m*(1+(m-2)*rho_y-(m-1)*rho_s*rho_y))
  sigma2_ate <- sigma2_y*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m)
  
  omega_d0 <- sigma2_ate + p_1^2*sigma2_hte
  omega_d1 <- sigma2_ate + p_0^2*sigma2_hte
  omega_d01 <- sigma2_ate - p_1*p_0*sigma2_hte
  
  pred.power <- 0
  n <- 4
  while (pred.power < 1-beta){
    n <- n + 2
    #non-centrality parameter (include the unknown of interest, n)
    delta_01 <- c(delta_0, delta_1)
    omega <- matrix(c(omega_d0, omega_d01, omega_d01, omega_d1), nrow=2, byrow=T)
    theta <- n*t(delta_01) %*% solve(omega) %*% delta_01
    pred.power <- pf(qf(1-alpha, 2, n-2), 2, n-2, ncp = theta, lower.tail = F)
  }
  return(c(n, pred.power))
}



################################################################################################
# Function to calculate the required number of clusters for the intersection-union test
################################################################################################

calc_IU <- function(alpha=0.05, beta=0.2, sigma2_y=1, beta_2, beta_4, m, rho_y, rho_s, p_i=0.5, p_1=0.5){
  delta_0 <- beta_2
  delta_1 <- beta_2 + beta_4
  p_0 <- 1 - p_1
  
  sigma2_hte <- sigma2_y*(1-rho_y)*(1+(m-1)*rho_y) / (p_i*(1-p_i)*p_1*p_0*m*(1+(m-2)*rho_y-(m-1)*rho_s*rho_y))
  sigma2_ate <- sigma2_y*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m)
  
  omega_d0 <- sigma2_ate + p_1^2*sigma2_hte
  omega_d1 <- sigma2_ate + p_0^2*sigma2_hte
  omega_d01 <- sigma2_ate - p_1*p_0*sigma2_hte
  
  pred.power <- 0
  n <- 4
  while (pred.power < 1-beta){
    n <- n + 2
    #non-centrality parameter (include the unknown of interest, n)
    mean_W <- c(delta_0/sqrt(omega_d0/n), delta_1/sqrt(omega_d1/n))
    cov_W <- (omega_d0)^(-1/2) * omega_d01 * (omega_d1)^(-1/2)
    sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)
    pred.power <- pmvt(df=n-2, lower=c(qt(1-alpha, n-2), qt(1-alpha, n-2)), upper=rep(Inf,2), delta=mean_W, sigma=sigma_W)
  }
  return(c(n, pred.power))
}



################################################################################################
# Function to calculate the required number of clusters for the interaction test
################################################################################################

calc_in <- function(alpha=0.05, beta=0.2, sigma2_y=1, beta_2, beta_4, m, rho_y, rho_s, p_i=0.5, p_1=0.5){
  sigma2_hte <- sigma2_y*(1-rho_y)*(1+(m-1)*rho_y) / (p_i*(1-p_i)*p_1*(1-p_1)*m*(1+(m-2)*rho_y-(m-1)*rho_s*rho_y))
  
  pred.power <- 0
  n <- 4
  while (pred.power < 1-beta){
    n <- n + 2
    #non-centrality parameter (include the unknown of interest, n)
    theta <- n * beta_4^2 / sigma2_hte
    pred.power <- pt(qt(alpha, df=n-2)+abs(beta_4)*sqrt(n/sigma2_hte), df=n-2)
  }
  return(c(n, pred.power))
}



################################################################################################
# Function to calculate the required number of clusters for testing the overall treatment effect
################################################################################################

calc_overall <- function(alpha=0.05, beta=0.2, sigma2_y=1, beta_2, beta_4, m, rho_y, rho_s, p_i=0.5, p_1=0.5){
  delta <- beta_2 + p_1*beta_4
  p_0 <- 1 - p_1
  
  sigma2_ate <- sigma2_y*(1+(m-1)*rho_y) / (p_i*(1-p_i)*m)
  
  pred.power <- 0
  n <- 4
  while (pred.power < 1-beta){
    n <- n + 2
    pred.power <- pt(qt(alpha, df=n-2)+abs(delta)*sqrt(n/sigma2_ate), df=n-2)
  }
  
  return(c(n, pred.power))
}


