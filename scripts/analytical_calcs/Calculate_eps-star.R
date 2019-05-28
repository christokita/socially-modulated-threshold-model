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
  eps_star$epsilon_all[i] <- delta / (tau/(1+tau) * E_01 + 1/(1+tau) * E_11)
}
eps_star$epsilon_all[eps_star$beta == 1.02]

ggplot(data = eps_star, aes(x = beta)) +
  geom_line(aes(y = epsilon_all)) +
  geom_line(aes(y = epsilon_inac), color = "lightblue") +
  geom_line(aes(y = epsilon_ac), color = "grey80") +
  theme_ctokita()

save(eps_star, beta_star, file = "output/AnalyticalResults/EpsStar-BetaStar_Calc.Rdata")

####################
# calculate epsilon* and beta* for n = 80 with infinite threshold limit (Method 2)
####################
n <- 80
delta <- 0.8
alpha <- 2
tau <- 0.2
beta <- 1.1
epsilon <- 0.55
fracs_n1 <- seq(0, 0.8, 0.01)

# Establish formulas
E01 <- function(n_1, n, beta) { n_1 / n +  n_1 * ( (1) / (beta * (n_1 - 1) + (n - n_1)) ) }
E02 <- function(n_2, n, beta) { n_2 / n +  n_2 * ( (1) / (beta * (n_2 - 1) + (n - n_2)) ) }

E11 <- function(n_1, n, beta) { 2 * beta * (n_1 - 1) / (beta * (n_1 - 1) + (n - n_1)) }
E22 <- function(n_2, n, beta) { 2 * beta * (n_2 - 1) / (beta * (n_2 - 1) + (n - n_2)) }

E12 <- function(n_1, n_2, n, beta) { n_2 / (beta * (n_1 - 1) + (n - n_1)) + n_2 * ( (1) / (beta * (n_2 - 1) + (n - n_2)) ) }
E21 <- function(n_1, n_2, n, beta) { n_1 / (beta * (n_2 - 1) + (n - n_2)) + n_1 * ( (1) / (beta * (n_1 - 1) + (n - n_1)) ) }

theta1 <- function(n_1, n_2, n, beta, tau) { (1/(1+tau)) * (E12(n_1, n_2, n, beta) - E11(n_1, n, beta)) + (tau/(1+tau)) * (E02(n_2, n, beta) - E01(n_1, n, beta)) }
theta2 <- function(n_1, n_2, n, beta, tau) { (1/(1+tau)) * (E21(n_1, n_2, n, beta) - E22(n_2, n, beta)) + (tau/(1+tau)) * (E01(n_1, n, beta) - E02(n_2, n, beta)) }
theta1_inac <- function(n_1, n_2, n, beta) { (E02(n_2, n, beta) - E01(n_1, n, beta)) }

# Calculate
delta_eps <- lapply(fracs_n1, function(fraction) {
  # set numbers
  n1 <- n * fraction
  n2 <- n * (0.8 - fraction)
  # Calcualte deltas
  delta1 <- delta - (alpha * fraction)
  delta2 <- delta - (alpha * (0.8 - fraction))
  # Calculate threhsold changes
  E_theta1 <- epsilon * theta1(n_1 = n1, n_2 = n2, n = n, beta = beta, tau = tau)
  E_theta2 <- epsilon * theta2(n_1 = n1, n_2 = n2, n = n, beta = beta, tau = tau)
  E_theta1inac <- epsilon * theta1_inac(n_1 = n1, n_2 = n2, n = n, beta = beta)
  # Return
  to_return <- data.frame(frac_n1 = fraction,
                          delta1 = delta1, 
                          delta2 = delta2,
                          E_theta1 = E_theta1,
                          E_theta2 = E_theta2,
                          E_theta1inac = E_theta1inac,
                          Diff1 = abs(E_theta1) - abs(delta1),
                          Diff2 = E_theta2 - delta2)
})
delta_eps <- do.call('rbind', delta_eps)

gg_deltaeps <- ggplot(data = delta_eps, aes(x = frac_n1)) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_vline(xintercept = 0.4, size = 0.3) +
  geom_line(aes(y = delta1), linetype = "dotted", color = "#1f78b4") +
  # geom_line(aes(y = delta2), linetype = "dotted", color = "#e31a1c") +
  geom_line(aes(y = E_theta1), color = "#1f78b4") +
  # geom_line(aes(y = E_theta2), color = "#e31a1c") +
  # geom_line(aes(y = E_theta1inac), color = "#a6cee3") +
  geom_line(aes(y = Diff1), linetype = "dashed", color = "#1f78b4") +
  # geom_line(aes(y = Diff2), linetype = "dashed", color = "#e31a1c") +
  theme_ctokita()

gg_deltaeps

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

theta1 <- function(n_1, n_2, n, beta, delta) { (1/(1+tau)) * (E12(n_1, n_2, n, beta) - E11(n_1, n, beta)) + (tau/(1+tau)) * (E02(n_2, n, beta) - E01(n_1, n, beta)) }
theta2 <- function(n_1, n_2, n, beta, delta) { (1/(1+tau)) * (E21(n_1, n_2, n, beta) - E22(n_2, n, beta)) + (tau/(1+tau)) * (E01(n_1, n, beta) - E02(n_2, n, beta)) }
theta1_inac <- function(n_1, n_2, n, beta, delta) { (E02(n_2, n, beta) - E01(n_1, n, beta)) }


