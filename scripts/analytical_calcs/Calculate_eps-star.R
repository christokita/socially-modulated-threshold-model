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

####################
# calculate epsilon* and beta* for n = 80 with infinite threshold limit
####################
# Let's assume all individuals are doing task 1 (or are inactive)
# Set parameters
tau <- 0.2
n <- 80
m <- 2
n_1 <- ( 1 / (1 + tau) ) * n
delta <- 0.8
beta <- seq(1, 1.25, 0.0001)

zero_corner <- 0.02

beta_star_value <- (delta * n) / (delta * n - m)

# Calculate beta* and epsilon*
beta_star <- data.frame(epsilon = seq(0, 0.6, 0.01),
                        beta = rep(beta_star_value, length(seq(0, 0.6, 0.01))))
eps_star <- data.frame(beta = beta, 
                       epsilon_inac = rep(NA, length(beta)), 
                       epsilon_ac = rep(NA, length(beta)),
                       epsilon_all = rep(NA, length(beta)))
for (i in 1:nrow(eps_star)) {
  E_01 <- n_1 * ( 1 / (beta[i] * (n_1 - 1) + (n - n_1)) ) # expected number of interaction partners performing task 1 if i is inactive
  E_11 <- 2*beta[i]*(n_1 - 1) / (beta[i] * (n_1 - 1) + (n - n_1)) # expected number of interaction partners performing task 1 if i is task 1
  eps_01 <- delta / E_01
  eps_11 <- delta / E_11
  eps_star$epsilon_inac[i] <- eps_01
  eps_star$epsilon_ac[i] <- eps_11
  eps_star$epsilon_all[i] <- (eps_11 * (1/(1+tau))) + (1-(1/(1+tau))) * eps_01
}

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
  # geom_hline(yintercept = zero_corner, size = 0.3) +
  # geom_line(data = beta_star_out, aes(x = beta, y = epsilon), size = 0.6, linetype = "dotted") +
  # geom_line(data = beta_star_in, aes(x = beta, y = epsilon), size = 0.6) +
  # geom_line(data = eps_star_out, aes(x = beta, y = epsilon_all), size = 0.6, linetype = "dotted") +
  # geom_line(data = eps_star_in, aes(x = beta, y = epsilon_all), size = 0.6) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1), expand = c(0, 0)) +
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
