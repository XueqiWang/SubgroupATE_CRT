#########################################################################
# Function to generate the simulated data sets for linear mixed modelling
#########################################################################

gendata <- function(beta1, beta2, beta3, beta4, rho_y, rho_s, sigma2y=1, n, m, p_1=0.5){
  
  #total number of individuals
  nm <- n*m
  
  #create individual identifications
  id <- seq(1, nm, 1)
  #assign cluster id
  cluster.id <- rep(1:n, each=m)
  data <- data.frame(cbind(id, cluster.id))
  
  #equally and randomly assign 0 or 1 to n clusters, representing the cluster-level treatment
  ctrt.assign <- sample(1:n, n/2)
  data$ctrt <- ifelse(data$cluster.id %in% ctrt.assign, 1, 0)
  
  #generate individual-level binary covariate
  q1 <- p_1*(1/rho_s-1)
  q2 <- (1-p_1)*(1/rho_s-1)
  p_c <- rbeta(n, q1, q2)
  s <- rbinom(n=nm, size=1, prob=rep(p_c, each=m))
  data$s <- s
  
  #randomly generate within-cluster variability for each individual
  data$epsilon <- rnorm(nm, 0, sqrt(1-rho_y))
  
  #randomly generate between-cluster variability for each cluster
  data$gamma <- rep(rnorm(n, 0, sqrt(rho_y)), each=m)
  
  #create outcome Y_ij's values
  data$Y <- beta1 + beta2*data$ctrt + beta3*data$s + beta4*data$s*data$ctrt + data$gamma + data$epsilon
  return(data)
}


