################################
#        Numerical study
################################

source("powerSubgroup.R")

library(mvtnorm)
library(ggplot2)
library(gridExtra)

data_plot <- function(rho_y_range, rho_s_range, p_1_range, test_type){
  data_all <- NULL
  for (rho_s in rho_s_range){
    for (rho_y in rho_y_range){
      for (p_1 in p_1_range){
        if (test_type == 1){
          pred.power <- as.numeric(power_omnibus(beta_2=beta_2, beta_4=beta_4, n=n, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1))[3]
        } else{
          pred.power <- as.numeric(power_IU(beta_2=beta_2, beta_4=beta_4, n=n, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1))[3]
        }
        data_all <- rbind(data_all, c(n=n, m=m, rho_y=rho_y, rho_s=rho_s, p_1=p_1, pred.power=pred.power))
      }
    }
  }
  return(as.data.frame(data_all))
}

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#########################
# Omnibus test
#########################

n <- 30
m <- 100
beta_2 <- 0.3
beta_4 <- 0.1

p_1_range <- c(0.3, 0.5, 0.7)
test_type <- 1

# Vary rho_y
rho_y_range <- seq(0.01, 0.5, 0.01)
rho_s_range <- c(0.2, 0.5)
pdata_om_1 <- data_plot(rho_y_range, rho_s_range, p_1_range, test_type)

rho_s_name <- c("0.2" = "(A)~Vary~rho[y*'|'*s]~with~rho[s]==0.2", "0.5" = "(B)~Vary~rho[y*'|'*s]~with~rho[s]==0.5")
pplot_om_1 <- ggplot(data=pdata_om_1, aes(rho_y, pred.power, colour = factor(p_1), shape = factor(p_1))) +
  geom_line(cex = 1.1, aes(linetype=factor(p_1))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2"),
                     breaks = c("0.3", "0.5", "0.7")) +
  labs(title = "", #"Power of the omnibus test",
       y = "Power", x = expression(rho[y*'|'*s]),
       colour = expression(paste("Prevalence of the subgroup indicator")),
       linetype = " ") +
  scale_x_continuous(breaks=seq(0.1, 0.5, 0.1)) +
  scale_y_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0.28, 1)) +
  theme_bw() +
  facet_grid(. ~ rho_s, scale="free_x", space="free_x",
             labeller = labeller(rho_s = as_labeller(rho_s_name,  label_parsed))) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15)) +
  labs(color  = "Prevalence of the subgroup indicator",
       linetype = "Prevalence of the subgroup indicator")

# Vary rho_s
rho_y_range <- c(0.05, 0.2)
rho_s_range <- seq(0.01, 0.99, 0.01)
pdata_om_2 <- data_plot(rho_y_range, rho_s_range, p_1_range, test_type)

rho_y_name <- c("0.05" = "(C)~Vary~rho[s]~with~rho[y*'|'*s]==0.05", "0.2" = "(D)~Vary~rho[s]~with~rho[y*'|'*s]==0.2")
pplot_om_2 <- ggplot(data=pdata_om_2, aes(rho_s, pred.power, colour = factor(p_1), shape = factor(p_1))) +
  geom_line(cex = 1.1, aes(linetype=factor(p_1))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2"),
                     breaks = c("0.3", "0.5", "0.7")) +
  labs(title = "", #"Power of the omnibus test",
       y = "Power", x = expression(rho[s]),
       colour = expression(paste("Prevalence of the subgroup indicator")),
       linetype = " ") +
  scale_x_continuous(breaks=seq(0, 1, 0.2)) +
  scale_y_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0.28, 1)) +
  theme_bw() +
  facet_grid(. ~ rho_y, scale="free_x", space="free_x",
             labeller = labeller(rho_y = as_labeller(rho_y_name,  label_parsed))) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15)) +
  labs(color  = "Prevalence of the subgroup indicator",
       linetype = "Prevalence of the subgroup indicator")

mylegend <- g_legend(pplot_om_1)
pplot_om <- grid.arrange(arrangeGrob(pplot_om_1 + theme(legend.position="none"),
                                     pplot_om_2 + theme(legend.position="none"),
                                     nrow=2),
                         mylegend, nrow=2, heights=c(12, 1))



#########################
# Intersection-union test
#########################

n <- 30
m <- 100
beta_2 <- 0.3
beta_4 <- 0.1

p_1_range <- c(0.3, 0.5, 0.7)
test_type <- 2

# Vary rho_y
rho_y_range <- seq(0.01, 0.5, 0.01)
rho_s_range <- c(0.2, 0.5)
pdata_IU_1 <- data_plot(rho_y_range, rho_s_range, p_1_range, test_type)

rho_s_name <- c("0.2" = "(A)~Vary~rho[y*'|'*s]~with~rho[s]==0.2", "0.5" = "(B)~Vary~rho[y*'|'*s]~with~rho[s]==0.5")
pplot_IU_1 <- ggplot(data=pdata_IU_1, aes(rho_y, pred.power, colour = factor(p_1), shape = factor(p_1))) +
  geom_line(cex = 1.1, aes(linetype=factor(p_1))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2"),
                     breaks = c("0.3", "0.5", "0.7")) +
  labs(title = "", #"Power of the intersection-union test",
       y = "Power", x = expression(rho[y*'|'*s]),
       colour = expression(paste("Prevalence of the subgroup indicator")),
       linetype = " ") +
  scale_x_continuous(breaks=seq(0.1, 0.5, 0.1)) +
  scale_y_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0.19, 1)) +
  theme_bw() +
  facet_grid(. ~ rho_s, scale="free_x", space="free_x",
             labeller = labeller(rho_s = as_labeller(rho_s_name,  label_parsed))) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15)) +
  labs(color  = "Prevalence of the subgroup indicator",
       linetype = "Prevalence of the subgroup indicator")

# Vary rho_s
rho_y_range <- c(0.05, 0.2)
rho_s_range <- seq(0.01, 0.99, 0.01)
pdata_IU_2 <- data_plot(rho_y_range, rho_s_range, p_1_range, test_type)

rho_y_name <- c("0.05" = "(C)~Vary~rho[s]~with~rho[y*'|'*s]==0.05", "0.2" = "(D)~Vary~rho[s]~with~rho[y*'|'*s]==0.2")
pplot_IU_2 <- ggplot(data=pdata_IU_2, aes(rho_s, pred.power, colour = factor(p_1), shape = factor(p_1))) +
  geom_line(cex = 1.1, aes(linetype=factor(p_1))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2"),
                     breaks = c("0.3", "0.5", "0.7")) +
  labs(title = "", #"Power of the intersection-union test",
       y = "Power", x = expression(rho[s]),
       colour = expression(paste("Prevalence of the subgroup indicator")),
       linetype = " ") +
  scale_x_continuous(breaks=seq(0, 1, 0.2)) +
  scale_y_continuous(breaks=seq(0.2, 1, 0.2), limits=c(0.19, 1)) +
  theme_bw() +
  facet_grid(. ~ rho_y, scale="free_x", space="free_x",
             labeller = labeller(rho_y = as_labeller(rho_y_name,  label_parsed))) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15)) +
  labs(color  = "Prevalence of the subgroup indicator",
       linetype = "Prevalence of the subgroup indicator")

mylegend <- g_legend(pplot_IU_1)
pplot_IU <- grid.arrange(arrangeGrob(pplot_IU_1 + theme(legend.position="none"),
                                     pplot_IU_2 + theme(legend.position="none"),
                                     nrow=2),
                         mylegend, nrow=2, heights=c(12, 1))


