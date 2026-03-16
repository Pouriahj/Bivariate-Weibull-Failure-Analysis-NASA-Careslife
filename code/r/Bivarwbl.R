# ============================================================
# Posterior visualization for the bivariate Weibull model
# ============================================================
# This script reproduces data for figures:
#   1. Autocorrelation plots for HMC chains
#   2. Density/Box Plot comparison between ABC-MH and ABC-HMC
#   3. 2D binned scatter plots for pairwise parameters
#
#
# Author: Pouria Hajizadeh
#
# ============================================================
#
# ---------------------------
# 1. Required libraries
# ---------------------------
library(rstan)
library(readxl)
library(bayesplot)
library(ggplot2)
library(R.matlab)
library(coda)
library(tidyr)
library(dplyr)
library(gridExtra)

# ---------------------------
# 2. Input files
# ---------------------------
# Experimental paired dataset (used when fitting Stan)
data_file <- "data/raw/data.xlsx"

# Stan model file
stan_file <- "stan/bivariate_weibull_hmc.stan"

# MATLAB file containing ABC-MH samples
# Assumption:
#   the .mat file contains a matrix named 'Freqq'
mh_mat_file <- "data/processed/ooutFreq.mat"

# ---------------------------
# 3. Load experimental data
# ---------------------------
dat <- read_excel(data_file)

# ---------------------------
# 4. Run / load HMC posterior from Stan
# ---------------------------
# Adjust warmup/iter/LENGTH as needed for your actual workflow.
tp1 <- 30

fit <- stan(
  file   = stan_file,
  data   = list(xy = dat, LENGTH = tp1),
  warmup = 100,
  iter   = 200,
  chains = 4
)

print(fit)

# Extract HMC posterior samples as matrix
hmc_samples <- as.matrix(fit)
colnames(hmc_samples) <- c("theta1", "theta2", "beta1", "beta2", "delta")

# ---------------------------
# 5. Load ABC-MH posterior samples from MATLAB
# ---------------------------
mat_data <- readMat(mh_mat_file)

# Assumes MATLAB variable name is 'Freqq'
mh_samples <- as.matrix(mat_data$Freqq)
colnames(mh_samples) <- c("theta1", "theta2", "beta1", "beta2", "delta")

# Convert MH samples to mcmc object when useful
mh_mcmc <- as.mcmc(mh_samples)

# ============================================================
# FIGURE 1: Autocorrelation plots for HMC chains
# ============================================================

bayesplot::color_scheme_set("brightblue")

acf_plot <- mcmc_acf(
  fit,
  pars = c("theta1", "theta2", "theta3", "theta4", "theta5"),
  lags = 20
) +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "serif", size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  ) +
  labs(title = "", x = "", y = "Autocorrelation")

print(acf_plot)

# ============================================================
# FIGURE 2: Density comparison between ABC-MH and ABC-HMC
# ============================================================

# Convert both posterior matrices to long format
hmc_df <- as.data.frame(hmc_samples) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(method = "ABC HMC")

mh_df <- as.data.frame(mh_samples) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(method = "ABC MH")

density_df <- bind_rows(hmc_df, mh_df)

# Nice labels for plotting
density_df$parameter <- factor(
  density_df$parameter,
  levels = c("theta1", "beta1", "theta2", "beta2", "delta"),
  labels = c(expression(theta[1]), expression(beta[1]), expression(theta[2]), expression(beta[2]), expression(delta))
)

density_plot <- ggplot(density_df, aes(x = value, color = method)) +
  geom_density(linewidth = 1.1) +
  facet_wrap(~ parameter, scales = "free", ncol = 3, labeller = label_parsed) +
  scale_color_manual(values = c("ABC HMC" = "deeppink", "ABC MH" = "blue")) +
  labs(x = "", y = "Probability Density", color = "") +
  theme_classic(base_family = "serif", base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

print(density_plot)

# ============================================================
# FIGURE 3: 2D binned scatter / posterior pair heatmaps
# ============================================================

# Parameter pairs shown in the paper-style plot
pair_list <- list(
  c("theta1", "theta2"),
  c("theta1", "beta1"),
  c("theta1", "beta2"),
  c("theta1", "delta"),
  c("theta2", "beta1"),
  c("theta2", "beta2"),
  c("theta2", "delta"),
  c("beta1",  "beta2"),
  c("beta1",  "delta"),
  c("beta2",  "delta")
)
