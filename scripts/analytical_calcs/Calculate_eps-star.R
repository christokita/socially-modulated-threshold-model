################################################################################
#
# Social interaction model: calculate epsilon* (analytical calculations)
#
################################################################################

rm(list = ls())

####################
# Source necessary scripts/libraries
####################
source("scripts/util/__Util__MASTER.R")
library(scales)

####################
# calculate epsilon* and beta* for n = 80 with infinite threshold limit
####################
# Let's assume all individuals are doing task 1 (or are inactive)
# Set parameters
tau <- 0.2
n <- 80
m <- 2
delta <- 0.8
beta <- seq(1, 1.25, 0.0001)
freq_activ <- delta
# freq_activ <- 1 / (1+tau)
n_1 <- freq_activ * n

beta_star_value <- (delta * n) / (delta * n - m)

# Calculate beta* and epsilon*
beta_star <- data.frame(epsilon = seq(0, 0.6, 0.01),
                        beta = rep(beta_star_value, length(seq(0, 0.6, 0.01))))
eps_star <- data.frame(beta = beta, 
                       epsilon_inac = rep(NA, length(beta)), 
                       epsilon_ac = rep(NA, length(beta)),
                       # epsilon_e11_e01_ratio = rep(NA, length(beta)),
                       epsilon_all = rep(NA, length(beta)))
for (i in 1:nrow(eps_star)) {
  E_01 <- n_1 / (n_1 + (n - n_1)) +  n_1 * ( (1) / (beta[i] * (n_1 - 1) + (n - n_1)) )  # expected number of interaction partners performing task 1 if i is inactive
  E_11 <- 2*beta[i]*(n_1 - 1) / (beta[i] * (n_1 - 1) + (n - n_1)) # expected number of interaction partners performing task 1 if i is task 1
  eps_01 <- delta / E_01
  eps_11 <- delta / E_11
  eps_star$epsilon_inac[i] <- eps_01
  eps_star$epsilon_ac[i] <- eps_11
  eps_star$epsilon_all[i] <- delta / (freq_activ * E_11 + (1 - freq_activ) * E_01)
}
eps_star$epsilon_all[eps_star$beta == 1.1]

# Calculate epsilon* for absolute loss of DOL (test?)
#for beta = 1.1, eps = 0.55, there seems to be points where there are as little as 10 individuals performing a task
n_1 <- delta/2 * n  + 1
n_2 <- delta/2 * n  - 1

eps_star_absolute <- data.frame(beta = beta,
                                E01 = rep(NA, length(beta)),
                                E02 = rep(NA, length(beta)),
                                E11 = rep(NA, length(beta)),
                                E12 = rep(NA, length(beta)),
                                E22 = rep(NA, length(beta)),
                                E21 = rep(NA, length(beta)),
                                epsilon_inac = rep(NA, length(beta)),
                                epsilon_ac = rep(NA, length(beta)),
                                epsilon_test = rep(NA, length(beta)),
                                # epsilon_e11_e01_ratio = rep(NA, length(beta)),
                                epsilon_all = rep(NA, length(beta)))
for (i in 1:nrow(eps_star_absolute)) {
  E_01 <- n_1 / (n_1 + (n - n_1)) +  n_1 * ( (1) / (beta[i] * (n_1 - 1) + (n - n_1)) )  # expected number of interaction partners performing task 1 if i is inactive
  E_02 <- n_2 / (n_2 + (n - n_2)) +  n_2 * ( (1) / (beta[i] * (n_2 - 1) + (n - n_2)) )  # expected number of interaction partners performing task 2 if i is inactive
  E_11 <- 2*beta[i]*(n_1 - 1) / (beta[i] * (n_1 - 1) + (n - n_1)) # expected number of interaction partners performing task 1 if i is task 1
  E_12 <- n_2 / (beta[i] * (n_1 - 1) + (n - n_1)) + n_2 * ( (1) / (beta[i] * (n_2 - 1) + (n - n_2)) )  # expected number of interaction partners performing task 1 if i is task 1
  E_22 <- 2*beta[i]*(n_2 - 1) / (beta[i] * (n_2 - 1) + (n - n_2)) # expected number of interaction partners performing task 2 if i is task 2
  E_21 <- n_1 / (beta[i] * (n_2 - 1) + (n - n_1)) + n_1 * ( (1) / (beta[i] * (n_1 - 1) + (n - n_1)) )  # expected number of interaction partners performing task 1 if i is task 2
  eps_star_absolute$E01[i] <- E_01 * (1 - delta)
  eps_star_absolute$E02[i] <- E_02 * (1 - delta)
  eps_star_absolute$E11[i] <- E_11 * delta
  eps_star_absolute$E12[i] <- E_12 * delta
  eps_star_absolute$E22[i] <- E_22 * delta
  eps_star_absolute$E21[i] <- E_21 * delta
  eps_delta_inac <- (E_01 - E_02)
  eps_delta_ac <-  (E_11 + E_21 - E_22 - E_12)
  eps_star_absolute$epsilon_inac[i] <- eps_delta_inac
  eps_star_absolute$epsilon_ac[i] <- eps_delta_ac
  eps_star_absolute$epsilon_test[i] <- (delta-(2*(n_2/n))) / (delta * (E_21 - E_22) + (1-delta) * (E_01 - E_02) )
  # eps_star_absolute$epsilon_test[i] <- delta*(E_11 + E_22) + (delta*(E_12 + E_21)) 
  # eps_star_absolute$epsilon_all[i] <- delta / (freq_activ * (E_11 - E_12) + (1 - freq_activ) * (E_01 - E_02))
  eps_star_absolute$epsilon_all[i] <- (delta-(2*(n_2/n))) / (delta*(eps_delta_ac) + (1-delta)*(eps_delta_inac))
}
eps_star_absolute$epsilon_test[eps_star_absolute$beta == 1.1] #should be ~0.54
eps_star_absolute$epsilon_test[eps_star_absolute$beta == 1.15] #should be ~0.6

eps_star_look <- eps_star_absolute %>% 
  select(-epsilon_inac, -epsilon_ac, -epsilon_all, -epsilon_test) %>% 
  melt(id.vars = "beta")
eps_star_look$value <- eps_star_look$value * 0.8

ggplot(data = eps_star_look, aes(x = beta, y = value, group = variable, col = variable)) +
  geom_hline(yintercept = (delta-(2*(n_2/n)))) +
  geom_line() +
  # scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "solid", "dotted")) +
  scale_colour_manual(values = c("#6a3d9a", "#cab2d6", "#1f78b4", "#fb9a99", "#a6cee3", "#e31a1c"))

ggplot(data = eps_star_absolute, aes(x = beta)) +
  # geom_hline(yintercept = (delta-(2*(n_2/n)))) +
  geom_line(aes(y = epsilon_test), color = "blue")
  

# Look at group size
Ns <- seq(0, n)
beta <- 1.2
epsilon <- 0.7
theta_shift <- data.frame(n2 = Ns,
                          d_delta2 = rep(NA, length(Ns)),
                          d_theta2 = rep(NA, length(Ns)),
                          d_thetadelt2 = rep(NA, length(Ns)),
                          d_delta1 = rep(NA, length(Ns)),
                          d_theta1 = rep(NA, length(Ns)),
                          d_thetadelt1 = rep(NA, length(Ns)),
                          rel_rate = rep(NA, length(Ns)))
for (i in 1:length(Ns)) {
  # Set variables
  n2 <- Ns[i]
  n1 <- delta*n - n2
  # Calculate values
  E01 <- n1/n + n1/(beta*(n1-1)+(n-n1))
  E02 <- n2/n + n2/(beta*(n2-1)+(n-n2)) 
  E21 <- n1/(beta*(n2-1)+(n-n2))  + n1/(beta*(n1-1)+(n-n1))
  E22 <- 2*beta*(n2-1)/(beta*(n2-1)+(n-n2))
  E12 <- n2/(beta*(n1-1)+(n-n1))  + n2/(beta*(n2-1)+(n-n2))
  E11 <- 2*beta*(n1-1)/(beta*(n1-1)+(n-n1))
  # Calculate expected change in theta
  d_theta2 <- delta*epsilon*(E21-E22) + (1-delta)*epsilon*(E01-E02)
  d_delta2 <- delta - m*n2/n
  d_theta1 <- delta*epsilon*(E12-E11) + (1-delta)*epsilon*(E02-E01)
  d_delta1 <- delta - m*n1/n
  # Return
  theta_shift$d_delta2[i] <- d_delta2
  theta_shift$d_theta2[i] <- d_theta2
  theta_shift$d_delta1[i] <- d_delta1
  theta_shift$d_theta1[i] <- -d_theta2
  # theta_shift$d_theta1[i] <- d_theta1
  theta_shift$d_thetadelt2[i] <- as.numeric(d_theta2 > d_delta2)
  theta_shift$d_thetadelt1[i] <- d_theta1 > d_delta1
  theta_shift$rel_rate[i] <-  theta_shift$d_thetadelt2[i] / theta_shift$d_thetadelt1[i]
}

ggplot(data = theta_shift, aes(x = n2)) +
  geom_vline(xintercept = delta/m*n, color = "grey30", linetype = "dashed") +
  geom_vline(xintercept = delta/m*n - 10, color = "grey30", linetype = "dotted") +
  geom_hline(yintercept = 0) +
  geom_line(aes(y = d_delta2), color = "lightblue") +
  geom_line(aes(y = d_theta2), color = "blue") +
  # geom_line(aes(y = d_delta1), color = "pink") +
  # geom_line(aes(y = d_theta1), color = "red") +
  # geom_line(aes(y = d_thetadelt2), color = "blue") +
  # geom_line(aes(y = d_thetadelt1), color = "red") +
  # geom_point(aes(y = d_thetadelt2), color = "red") +
  theme_classic()

####################
# Plot analytical "heatmap" equivalent figure
####################

# Corner points for epsilon* polygon
eps_corner1 <- data.frame(beta = 1, 
                            epsilon_inac = NA, 
                            epsilon_ac = NA,
                            epsilon_all = 0.60)
eps_corner2 <- data.frame(beta = 1.25, 
                            epsilon_inac = NA, 
                            epsilon_ac = NA,
                            epsilon_all = 0.60)
eps_poly <- eps_corner1 %>% 
  rbind(eps_star) %>% 
  rbind(eps_corner2)
eps_poly$DOL <- 0

# Polygon for zero DOL near zero epsilon
eps_low_poly <- data.frame(beta = c(1, 1, 1.25, 1.25), 
                           epsilon_inac = NA, 
                           epsilon_ac = NA,
                           epsilon_all = c(0, zero_corner, zero_corner, 0),
                           DOL = 0)

# Coner points for beta* points
beta_corner1 <- data.frame(epsilon = 0,
                           beta = 1)
beta_corner2 <- data.frame(epsilon = 0.60,
                           beta = 1)
beta_poly <- beta_corner1 %>% 
  rbind(beta_star) %>% 
  rbind(beta_corner2) %>% 
  mutate(Id = "BetaStar")
beta_poly$DOL <- 0

# DOL polygon
DOL_poly <- data.frame(beta = c(unique(beta_star$beta),
                                eps_star$beta[eps_star$beta >= unique(beta_star$beta)],
                                1.25),
                       epsilon = c(0, eps_star$epsilon_all[eps_star$beta >= unique(beta_star$beta)], 0))
DOL_poly$DOL <- 1

# Filter epsilon* and beta* lines 
intersection_point <- max(eps_star$epsilon_all[eps_star$beta <= beta_star_value])
eps_star_out <- eps_star %>% 
  filter(beta < beta_star_value)
eps_star_in <- eps_star %>% 
  filter(beta >= beta_star_value)
beta_star_out <- beta_star %>%
  filter(epsilon > intersection_point)
beta_star_in <- beta_star %>%
  filter(epsilon <= intersection_point)

# Plot
pal <- brewer_pal("seq", "GnBu")
pal <- pal(9)

gg_analytical_betaeps <- ggplot() +
  geom_polygon(data = eps_poly, aes(x = beta, y = epsilon_all, fill = DOL)) +
  geom_polygon(data = beta_poly, aes(x = beta, y = epsilon, fill = DOL)) +
  geom_polygon(data = DOL_poly, aes(x = beta, y = epsilon, fill = DOL)) +
  # geom_polygon(data = eps_low_poly, aes(x = beta, y = epsilon_all, fill = DOL)) +
  geom_line(data = beta_star, aes(x = beta, y = epsilon), size = 0.3, linetype = "dashed") +
  geom_line(data = eps_star, aes(x = beta, y = epsilon_all), size = 0.3) +
  geom_line(data = eps_star_absolute, aes(x = beta, y = epsilon_all), size = 0.3, linetype = "dotted") +
  # geom_hline(yintercept = zero_corner, size = 0.3) +
  # geom_line(data = beta_star_out, aes(x = beta, y = epsilon), size = 0.6, linetype = "dotted") +
  # geom_line(data = beta_star_in, aes(x = beta, y = epsilon), size = 0.6) +
  # geom_line(data = eps_star_out, aes(x = beta, y = epsilon_all), size = 0.6, linetype = "dotted") +
  # geom_line(data = eps_star_in, aes(x = beta, y = epsilon_all), size = 0.6) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1), expand = c(0, 0)) +
  scale_fill_gradientn(colours = pal, name = "Behavioral\nspecialization",
                       limits = c(0, 1)) +
  xlab(expression(paste("Interaction Bias (", italic(beta), ")"))) +
  ylab(expression(paste( "Social influence (", italic(epsilon), ")"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.3, color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        legend.position = "none")
gg_analytical_betaeps

ggsave(gg_analytical_betaeps, file ="output/AnalyticalResults/BetaEps_AnalyticalSpace.png", width = 45, height = 45, units = "mm", dpi = 300)
ggsave(gg_analytical_betaeps, file ="output/AnalyticalResults/BetaEps_AnalyticalSpace.svg", width = 45, height = 45, units = "mm")
