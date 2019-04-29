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
# calculate expected threshold change for active and inactive individual
####################
# Let's assume all individuals are doing task 1 (or are inactive)
# Set parameters
tau <- 0.2
n <- 80
m <- 2
n_1 <- ( 1 / (1 + tau) ) * n
delta <- 0.8
beta <- seq(1, 1.25, 0.001)

beta_star_value <- (delta * n) / (delta * n - m)

# Calculate
beta_star <- data.frame(epsilon = seq(0, 0.65, 0.01),
                        beta_star = rep(beta_star_value, length(seq(0, 0.65, 0.01))))
eps_star <- data.frame(beta = beta, 
                       epsilon_inac = rep(NA, length(beta)), 
                       epsilon_ac = rep(NA, length(beta)),
                       epsilon_all = rep(NA, length(beta)))
for (i in 1:nrow(eps_star)) {
  E_01 <- n_1 * ( 1 / (beta[i] * (n_1 - 1) + (n - n_1)) )
  E_11 <- 2*beta[i]*(n_1 - 1) / (beta[i] * (n_1 - 1) + (n - n_1))
  eps_01 <- delta / E_01
  eps_11 <- delta / E_11
  eps_star$epsilon_inac[i] <- eps_01
  eps_star$epsilon_ac[i] <- eps_11
  eps_star$epsilon_all[i] <- (eps_11 * (1/(1+tau))) + (1-(1/(1+tau))) * eps_01
}

# Plot
ggplot(data = eps_star, aes(x = beta)) +
  # geom_line(aes(y = epsilon_inac), color = "blue") +
  # geom_line(aes(y = epsilon_ac), color = "red") +
  geom_line(aes(y = epsilon_all), linetype = "dashed") +
  geom_line(data = beta_star, aes(x = beta_star, y = epsilon)) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.65), breaks = seq(0, 0.6, 0.1), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

