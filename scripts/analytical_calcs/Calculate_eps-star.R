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
  eps_star$epsilon_all[i] <- delta / E_01
}
eps_star$epsilon_all[eps_star$beta == 1.02]

ggplot(data = eps_star, aes(x = beta)) +
  geom_line(aes(y = epsilon_all)) +
  geom_line(aes(y = epsilon_ac), color = "grey80") +
  theme_ctokita()

save(eps_star, beta_star, file = "output/AnalyticalResults/EpsStar-BetaStar_Calc.Rdata")


####################
# calculate epsilon* and beta* for n = 80 with infinite threshold limit (Method 2)
####################
# Calculate epsilon* for absolute loss of DOL (test?)

# Establish formulas
E01 <- function(n_1, n, beta) { n_1 / (n_1 + (n - n_1)) +  n_1 * ( (1) / (beta * (n_1 - 1) + (n - n_1)) ) }
E02 <- function(n_2, n, beta) { n_2 / (n_2 + (n - n_2)) +  n_2 * ( (1) / (beta * (n_2 - 1) + (n - n_2)) ) }
E11 <- function(n_1, n, beta) { 2 * beta * (n_1 - 1) / (beta * (n_1 - 1) + (n - n_1)) }
E22 <- function(n_2, n, beta) { 2 * beta * (n_2 - 1) / (beta * (n_2 - 1) + (n - n_2)) }
E12 <- function(n_1, n_2, n, beta) { n_2 / (beta * (n_1 - 1) + (n - n_1)) + n_2 * ( (1) / (beta * (n_2 - 1) + (n - n_2)) ) }
E21 <- function(n_1, n_2, n, beta) { n_1 / (beta * (n_2 - 1) + (n - n_1)) + n_1 * ( (1) / (beta * (n_1 - 1) + (n - n_1)) ) }

theta1 <- function(n_1, n_2, n, beta, delta) { delta * (E12(n_1, n_2, n, beta) - E11(n_1, n, beta)) + (1 - delta) * (E02(n_2, n, beta) - E01(n_1, n, beta)) }
theta2 <- function(n_1, n_2, n, beta, delta) { delta * (E21(n_1, n_2, n, beta) - E22(n_2, n, beta)) + (1 - delta) * (E01(n_1, n, beta) - E02(n_2, n, beta)) }
theta1_inac <- function(n_1, n_2, n, beta, delta) { (E02(n_2, n, beta) - E01(n_1, n, beta)) }

# Sweep by beta, keeping n fixed
Ns <- seq(30, 100, 1)
m <- 2
a <- m
betas <- seq(1, 1.25, 0.01)
epsilons <- seq(0, 0.8, 0.01)
delta <- 0.8

eps_star_sweep <- lapply(Ns, function(n) {
  eps_star_beta <- lapply(betas, function(beta) {
    # set n difference
    n_1_equal <- delta/m*n 
    n_2_equal <- n_1_equal
    n_1_diff <- n_1_equal + 1
    n_2_diff <- n_1_equal - 1
    # Calculate expected threshold change per individual
    # change0 <- theta1(n_1 = n_1_equal, n_2 = n_2_equal, n = n, beta = beta, delta = delta)
    # change1 <- theta1(n_1 = n_1_diff, n_2 = n_2_diff, n = n, beta = beta, delta = delta)
    change0 <- theta1_inac(n_1 = n_1_equal, n_2 = n_2_equal, n = n, beta = beta, delta = delta)
    change1 <- theta1_inac(n_1 = n_1_diff, n_2 = n_2_diff, n = n, beta = beta, delta = delta)
    theta_diff <- change1 - change0
    # Calculate expected delta change per individual
    delta_diff <- -a / n
    # Return
    to_return <- data.frame(beta = beta, 
                            theta_diff = theta_diff,
                            delta_diff,
                            eps_star = delta_diff / theta_diff)
    return(to_return)
  })
  eps_star_beta <- do.call('rbind', eps_star_beta) %>% 
    mutate(n = n)
})
eps_star_sweep <- do.call('rbind', eps_star_sweep)

ggplot(eps_star_sweep, aes(y = beta, x = n, fill = eps_star)) +
  geom_tile() +
  theme_ctokita()


####################
# calculate epsilon* and beta* for n = 80 with infinite threshold limit (Method 3)
####################
# Calculate epsilon* for absolute loss of DOL (test?)
# Establish formulas
E01 <- function(n_1, n, beta) { n_1 / (n_1 + (n - n_1)) +  n_1 * ( (1) / (beta * (n_1 - 1) + (n - n_1)) ) }
E02 <- function(n_2, n, beta) { n_2 / (n_2 + (n - n_2)) +  n_2 * ( (1) / (beta * (n_2 - 1) + (n - n_2)) ) }
E11 <- function(n_1, n, beta) { 2 * beta * (n_1 - 1) / (beta * (n_1 - 1) + (n - n_1)) }
E22 <- function(n_2, n, beta) { 2 * beta * (n_2 - 1) / (beta * (n_2 - 1) + (n - n_2)) }
E12 <- function(n_1, n_2, n, beta) { n_2 / (beta * (n_1 - 1) + (n - n_1)) + n_2 * ( (1) / (beta * (n_2 - 1) + (n - n_2)) ) }
E21 <- function(n_1, n_2, n, beta) { n_1 / (beta * (n_2 - 1) + (n - n_1)) + n_1 * ( (1) / (beta * (n_1 - 1) + (n - n_1)) ) }

theta1 <- function(n_1, n_2, n, beta, delta) { delta * (E12(n_1, n_2, n, beta) - E11(n_1, n, beta)) + (1 - delta) * (E02(n_2, n, beta) - E01(n_1, n, beta)) }
theta2 <- function(n_1, n_2, n, beta, delta) { delta * (E21(n_1, n_2, n, beta) - E22(n_2, n, beta)) + (1 - delta) * (E01(n_1, n, beta) - E02(n_2, n, beta)) }
theta1_inac <- function(n_1, n_2, n, beta, delta) { (E02(n_2, n, beta) - E01(n_1, n, beta)) }

# Sweep by beta, keeping n fixed
# Ns <- seq(30, 100, 1)
n <- 80
m <- 2
a <- m
betas <- seq(1, 1.25, 0.1)
epsilons <- seq(0, 0.8, 0.05)
delta <- 0.8

eps_star_sweep <- lapply(epsilons, function(epsilon) {
  eps_star_beta <- lapply(betas, function(beta) {
    # set n difference
    n_1 <- delta/m*n  + 1
    n_2 <- delta/m*n - 1
    # Calculate expected threshold change per individual
    d_theta1 <- theta1(n_1 = n_1, n_2 = n_2, n = n, beta = beta, delta = delta)
    d_theta2 <- theta1(n_1 = n_1, n_2 = n_2, n = n, beta = beta, delta = delta)
    d_inac <- theta1_inac(n_1 = n_1, n_2 = n_2, n = n, beta = beta, delta = delta)
    # Calculate average threshold change
    avg_d_theta <- delta/2*d_theta1 - delta/2*d_theta2 + (1-delta)*d_inac
    # Calculate delta change
    d_delta1 <- delta - m * n_1/n
    d_delta2 <- delta - m * n_2/n
    avg_d_delta <- (d_delta1 + d_delta2) / 2
    # Return
    to_return <- data.frame(beta = beta, 
                            d_theta1 = d_theta1 * epsilon,
                            d_theta2 = d_theta2 * epsilon,
                            d_inac = d_inac * epsilon,
                            d_delta1 =  d_delta1,
                            d_delta2 =  d_delta2,
                            avg_d_theta = avg_d_theta * epsilon,
                            avg_d_delta = avg_d_delta)
    return(to_return)
  })
  eps_star_beta <- do.call('rbind', eps_star_beta) %>% 
    mutate(epsilon = epsilon)
})





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
