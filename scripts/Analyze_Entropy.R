################################################################################
#
# Comparing various specialization plots
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)

####################
# Load data
####################
load("output/Rdata/_ProcessedData/Entropy/Sigma0-Epsilon0.1-Beta1.1.Rdata")
compiled_data$Model <- "Social"
entropy_data <- compiled_data
rm(compiled_data)

load("output/Rdata/_ProcessedData/Entropy/Sigma0.05-Epsilon0-Beta1.1.Rdata")
compiled_data$Model <- "Fixed_Sigma0.05"
entropy_data <- rbind(entropy_data, compiled_data)
rm(compiled_data)


####################
# Summarise data
####################
# Calculate mean and SE
entropy <- entropy_data %>% 
  group_by(Model, n) %>% 
  summarise(Mean = mean(Dind),
            SE = sd(Dind) / sqrt(length(Dind)))


####################
# Plot
####################
gg_entropy <- ggplot(data = entropy, aes(x = n, colour = Model)) +
  geom_line(aes(y = Mean),
            size = 0.4) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                width = 0) +
  geom_point(aes(y = Mean),
             size = 2) +
  theme_classic() +
  ylab("Division of Labor") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme(axis.text = element_text(colour = "black", size = 6),
        axis.title = element_text(size = 7),
        aspect.ratio = 1)
gg_entropy
