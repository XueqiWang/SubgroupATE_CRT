################################
# Application to the UMDEX study
################################

source("calcSubgroup.R")
source("powerSubgroup.R")

library(mvtnorm)
library(ggplot2)
library(directlabels)
library(cowplot)

# main assumptions
m <- 10
rho_y <- 0.04
rho_s <- 0.2
p_1 <- 0.36

delta_0 <- 0.7
delta_1 <- 0.5
beta_2 <- delta_0
beta_4 <- delta_1 - delta_0

# omnibus test
np_omnibus <- calc_omnibus(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1)
np_omnibus
n_omnibus <- np_omnibus[1] # 18
pred.power_omnibus <- np_omnibus[2]
# back-of-the-envelope approach
np_omnibus_b <- calc_omnibus(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=0, rho_s=0, p_1=p_1)
nc_omnibus_b <- np_omnibus_b[1]*(1+(m-1)*rho_y)
nc_omnibus_b <- ceiling(nc_omnibus_b)
if (nc_omnibus_b%%2 != 0){
  nc_omnibus_b <- nc_omnibus_b + 1
}
nc_omnibus_b # 20
pred.power_omnibus_b <- as.numeric(power_omnibus(beta_2=beta_2, beta_4=beta_4, n=nc_omnibus_b, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1))[3]
pred.power_omnibus_b

# intersection-union test
np_IU <- calc_IU(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1)
np_IU
n_IU <- np_IU[1] # 34
pred.power_IU <- np_IU[2]
# back-of-the-envelope approach
np_IU_b <- calc_IU(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=0, rho_s=0, p_1=p_1)
nc_IU_b <- np_IU_b[1]*(1+(m-1)*rho_y)
nc_IU_b <- ceiling(nc_IU_b)
if (nc_IU_b%%2 != 0){
  nc_IU_b <- nc_IU_b + 1
}
nc_IU_b # 42
pred.power_IU_b <- as.numeric(power_IU(beta_2=beta_2, beta_4=beta_4, n=nc_IU_b, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1))[3]
pred.power_IU_b

# interaction test
np_in <- calc_in(beta_2=beta_2, beta_4=beta_4, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1)
np_in
n_in <- np_in[1] # 284
pred.power_in <- np_in[2]


conPLOT <- function(m, rho_y_range, rho_s_range, test_type){
  con_plot <- NULL
  umdex <- NULL
  for (rho_y in seq(rho_y_range[1], rho_y_range[2], 0.001)){
    for (rho_s in seq(rho_s_range[1], rho_s_range[2], 0.05)){
      if (test_type == 1){
        umdex<-rbind(umdex, power_omnibus(beta_2=beta_2, beta_4=beta_4, n=n_omnibus, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1))
      } else {
        umdex<-rbind(umdex, power_IU(beta_2=beta_2, beta_4=beta_4, n=n_IU, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1))
      }
    }
  }
  fig <- ggplot()  +
    theme_bw() +
    ggtitle(bquote(m == ~.(m))) +
    xlab(expression(rho[y*'|'*s])) +
    ylab(expression(rho[s])) +
    ylim(c(rho_s_range[1], rho_s_range[2])) +
    scale_x_continuous(breaks=seq(rho_y_range[1], rho_y_range[2], 0.02), limits=c(0.01, 0.09)) +
    stat_contour(data = umdex, aes(x = rho_y, y = rho_s, z = t_test, colour = ..level..), 
                 breaks = round(quantile(umdex$t_test, seq(0, 1, 0.12)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(1, 0), legend.position=c(1, 0),
          text = element_text(size = 15))
  con_plot <- direct.label(fig, "bottom.pieces")
  return(con_plot)
}


# figure of the sensitivity analysis for the omnibus test
test_type <- 1

umdex_1_1<-conPLOT(10, c(0.01, 0.09), c(0, 1), 1)
umdex_1_2<-conPLOT(20, c(0.01, 0.09), c(0, 1), 1)
umdex_1_3<-conPLOT(30, c(0.01, 0.09), c(0, 1), 1)
umdex_1_4<-conPLOT(40, c(0.01, 0.09), c(0, 1), 1)
plot_grid(umdex_1_1, umdex_1_2, umdex_1_3, umdex_1_4, ncol=2)


# figure of the sensitivity analysis for the intersection-union test
test_type <- 2

umdex_2_1<-conPLOT(10, c(0.01, 0.09), c(0, 1), 2)
umdex_2_2<-conPLOT(20, c(0.01, 0.09), c(0, 1), 2)
umdex_2_3<-conPLOT(30, c(0.01, 0.09), c(0, 1), 2)
umdex_2_4<-conPLOT(40, c(0.01, 0.09), c(0, 1), 2)
plot_grid(umdex_2_1, umdex_2_2, umdex_2_3, umdex_2_4, ncol=2)


