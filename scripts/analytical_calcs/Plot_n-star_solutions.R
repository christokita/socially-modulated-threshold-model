################################################################################
#
# Plot N star solutions for full solution vs simplified solution
#
################################################################################
rm(list = ls())
source("scripts/util/__Util__MASTER.R")

betas <- seq(1, 2, 0.001)
delta <- 0.8

nstar <- lapply(betas, function(beta) {
  full_solution_nom <- (4*beta-delta-4*beta*delta+5*beta^2*delta+sqrt(16*beta^2-8*beta*delta+16*beta^2*delta-8*beta^3*delta+delta^2+8*beta*delta^2-18*beta^2*delta^2+8*beta^3*delta^2+beta^4*delta^2))
  full_solution_denom <- 2*(-2*delta+2*beta*delta+delta^2-2*beta*delta^2+beta^2*delta^2)
  full_solution <- full_solution_nom / full_solution_denom
  simp_solution <- 2*beta / (delta * (beta-1))
  to_return <- data.frame(beta = beta, nstar_simp = simp_solution, nstar_full = full_solution)
})
nstar <- do.call('rbind', nstar)
nstar_melt <- nstar %>% 
  melt(id.vars = "beta")


# Plot
gg_nstar_comp <- ggplot(nstar_melt, aes(x = beta, group = variable, color = variable)) +
  geom_line(aes(y = value), size = 0.3) +
  ylab("n*") +
  xlab(expression(paste("Interaction bias (", beta, ")"))) +
  scale_y_continuous(limits = c(0, 30)) +
  scale_color_manual(name = "Solution", values = c("#0d75ff", "#d60036"), labels = c("Simplified", "Full")) +
  theme_ctokita() +
  theme(aspect.ratio = 1,
        legend.position = "right")
gg_nstar_comp

ggsave(gg_nstar_comp, file = "output/AnalyticalResults/FullVsSimplifiedSolution.png", dpi = 400, width = 95, height = 65, units = "mm")
