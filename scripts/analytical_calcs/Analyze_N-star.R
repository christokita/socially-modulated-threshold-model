################################################################################
#
# Analyze Simulated N star vs Analytical N star
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")


####################
# Load and process data
####################
# Load DOL values
load("output/AnalyticalResults/CalculateNstar-50_Sigma0-Epsilon0.1.Rdata")

# Summarise
sim_data <- parallel_simulations %>% 
  as.data.frame(.) %>% 
  group_by(beta, n, sim) %>% 
  mutate(DOL_change = Dind - lag(Dind, default = first(Dind), order_by = t)) %>% 
  filter(t == 6000) %>% 
  ungroup() %>% 
  group_by(beta, n) %>% 
  # select(beta, n, sim, DOL_change) %>% 
  summarise(DOL_change = mean(DOL_change))

# Calcualte N-star
betas <- unique(sim_data$beta)
sim_data$N_star_range <- FALSE
for (b in 1:length(betas)) { #find the rows within a beta range where DOL_change goes from negative to positive
  first_pos <- min(which(sim_data$DOL_change[sim_data$beta == betas[b]] > 0))
  set_of_rows_start <- (b-1) * 11 #account for which set of rows the indezxing is occuring in
  sim_data$N_star_range[set_of_rows_start + first_pos - 1] <- TRUE 
  sim_data$N_star_range[set_of_rows_start + first_pos] <- TRUE 
}
sim_data <- sim_data %>% 
  filter(N_star_range == TRUE)
Nstar_data <- sim_data %>% 
  select(beta) %>% 
  unique(.) %>% 
  mutate(N_star = NA)
for(b in betas) {
  lower_n <- min(sim_data$n[sim_data$beta == b])
  upper_n <- max(sim_data$n[sim_data$beta == b])
  lower_DOL <- min(sim_data$DOL_change[sim_data$beta == b])
  upper_DOL <- max(sim_data$DOL_change[sim_data$beta == b])
  slope <- (upper_DOL - lower_DOL) / (upper_n - lower_n)
  n_star <- -lower_DOL / slope #find x intercept
  n_star <- n_star + lower_n #make actual n*
  Nstar_data$N_star[Nstar_data$beta == b] <- n_star
}

# calculate
analytical_data <- data.frame(beta = seq(1.001, 1.255, 0.0001), n = rep(NA, length(seq(1.001, 1.255, 0.0001))))
for (i in 1:nrow(analytical_data)) {
  analytical_data[i, 2] <- (2 * analytical_data[i, 1]) / (0.8 * (analytical_data[i, 1] - 1))
}
analytical_data <- analytical_data %>% 
  filter(n < 102.5)

####################
# Plot
####################
gg_nstar <- ggplot() +
  geom_line(data = analytical_data,
            aes(x = beta, y = n),
            size = 0.2, 
            colour = "#fb9a99") +
  geom_point(data = Nstar_data,
             aes(x = beta, y = N_star),
             size = 0.4,
             colour = "#e31a1c") +
  xlab(expression(paste("Interaction bias ", italic(beta)))) +
  ylab(expression(paste("n"^"*"))) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  theme_ctokita() +
  theme(aspect.ratio = 1)
gg_nstar

ggsave("output/AnalyticalResults/NstarCalculation.png", width = 45, heigh = 45, units = "mm", dpi = 400)
ggsave("output/AnalyticalResults/NstarCalculation.svg", width = 45, heigh = 45, units = "mm")

