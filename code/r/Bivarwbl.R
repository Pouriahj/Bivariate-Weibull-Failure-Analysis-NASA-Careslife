library(rstan);
library(readxl);
library(plotly);
library(bayesplot);
library(ggplot2);
library(rstanarm);
library(R.matlab);
library(ggmcmc);
library(coda);


dat = read_excel("data/raw/data.xlsx")

# Create the bubble chart
ggplot(dat, aes(x = Strength, y = TTF, size = Time_to_failure)) +
  geom_point(shape = 21, fill = "pink") +
  scale_size_area(max_size = 15) +
  ggtitle("Bubble Chart of Paired Data") +
  xlab("Strength") +
  ylab("Time to Failure")


d1 = dat[ , 1];
d2 = dat[ , 2];

tp1 = 30;

fit = stan(file = "bivarweibullt1.stan",
           data = list(xy = dat, LENGTH = tp1),
           warmup = 100,
           iter = 200,
           chains = 4);

print(fit)

mcmc_chain = as.matrix(fit)

plot(density(mcmc_chain[,'theta5']))

write.csv(mcmc_chain, "matrix_from_fit.csv")

# Defining new parameter names for fit object

thetan1 <- paste("\u03B8", "\u2081", sep = "");
thetan2 <- paste("\u03B8", "\u2082", sep = "");
thetan3 <- paste("\u03B2", "\u2081", sep = "");
thetan4 <- paste("\u03B2", "\u2082", sep = "");
thetan5 <- "\u03B4"

bayesplot::color_scheme_set("brightblue")

new_names <- c("thetan1" = "theta1", "thetan2" = "theta2", "thetan3" = "theta3", "thetan4" = "theta4", "thetan5" = "theta5")
fit_renamed <- rename_vars(fit, new_names)


fit_matrix <- as.matrix(fit)
colnames(fit_matrix) <- c(thetan1, thetan2, thetan3, thetan4, thetan5)
mcmc_areas(fit_matrix)

# Rename columns to desired labels
colnames(params) <- c("thetan1", "thetan2", "thetan3", "thetan4", "thetan5")

# Plot the renamed parameters

# Plot credible intervals for parameters
mcmc_areas(fit, prob = 0.9, alpha = 0.3, pars = c("theta1", "theta3")) 
mcmc_areas(fit, prob = 0.9, alpha = 0.3, pars = c("theta2", "theta4"))
mcmc_areas(fit, prob = 0.9, alpha = 0.3, pars = c("theta5"))

# Read the .mat file
mat_data <- readMat("ooutFreq.mat")

# Extract the desired rows from the matrix
data <- as.matrix(mat_data$Freqq)
mat_mcmc = as.mcmc(data);

colnames(data) <- c("theta1", "theta2", "theta3", "theta4", "theta5")

# Plot the data using mcmc_areas
bayesplot::color_scheme_set("pink")
mcmc_areas(data)

# This one wasn't good
color_scheme_set("darkgray")
mcmc_parcoord(fit)

# This one doesn't work out
bayesplot::color_scheme_set("red")
ppc_dens_overlay(y = mcmc_chain[ , "theta1"] , yrep = posterior_predict(fit, draws = 50))

colnames(data) <- c("theta1", "theta2","theta3", "theta4", "theta5")
df <- as.data.frame(data)
bayesplot::color_scheme_set("pink")
mcmc_pairs(df, pars = colnames(df), off_diag_args = list(size = 1.5))



bayesplot::color_scheme_set("pink")
mcmc_pairs(fit, pars = c("theta1", "theta2","theta3", "theta4", "theta5"),
           off_diag_args = list(size = 1.5))

bayesplot::color_scheme_set("pink")
mcmc_pairs(mat_mcmc, pars = c("theta1", "theta2","theta3", "theta4", "theta5"),
           off_diag_args = list(size = 1.5))
# bayesplot_theme_update(panel.background = element_rect(fill = "white"),
                       #panel.border = element_blank(),
                       #panel.grid.minor = element_blank(),
                       #panel.grid.major = element_blank(),
                       #axis.line = element_line(colour = "white"));

# This visualization is good enough | Run it to see
color_scheme_set("brightblue")
mcmc_combo(fit, pars = c("theta1", "theta2", "theta3", "theta4", "theta5"))
bayesplot_theme_update(panel.border = element_rect(size=0.5, colour="black"),
                       panel.grid.minor = element_blank(),
                       panel.grid.major = element_line(colour = "gray"),
                       axis.line = element_line(size=0.5, colour = "black"));


bayesplot::color_scheme_set("viridis")
mcmc_trace_highlight(fit, pars = "theta1", highlight = 10)



x <- example_mcmc_draws()

color_scheme_set("mix-blue-red")
mcmc_combo(
  fit,
  combo = c("dens_overlay", "trace"),
  pars = c("theta1", "theta3"),
  transformations = list(theta1 = "log"),
  gg_theme = legend_none()
)


color_scheme_set("red")
mcmc_scatter(
  as.matrix(fit),
  pars = c("theta1", "theta3"), 
  np = nuts_params(fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)



fit2 = stan_glm(mpg ~ ., data = dat)
posterior = as.matrix(fit2)

plot_title = ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("Strength", "TTF"),
           prob = 0.8) + plot_title




# Visualizing 3D surface of CDF function

theta_1 = 42.29
theta_2 = 7.16
beta_1 = 9.93
beta_2 = 2.37
delta = 0.7

F <- function(x, y) {
  1 - exp(-((x/theta_1)^(beta_1/delta) + (y/theta_2)^(beta_2/delta))^delta)
}

x <- seq(0, 70, length.out = 100)
y <- seq(0, 10, length.out = 100)
z <- outer(x, y, F)

library(RColorBrewer)
library(viridis)
pink_palette <- brewer.pal(9, "PuBu")

persp(x, y, z, theta = 30, phi = 30, col = pink_palette, shade = 0.75, ltheta = -120,
      ticktype = "detailed", type = 's', d = 3, lphi = 0, expand = 0.5, sub = "", xlab = "Strength", ylab = "Time to Failure", zlab = "Probability of Failure")



# Load ggplot2 library
library(ggplot2)

theme_set(theme_classic())

# Use mcmc_chain as the data for the plot
df = data.frame(parameter = rep(c("theta1", "theta2", "theta3", "theta4", "theta5")
                                , each = nrow(mcmc_chain)), 
                value = c(mcmc_chain[, "theta1"], mcmc_chain[, "theta2"]
                          , mcmc_chain[, "theta3"], mcmc_chain[, "theta4"]
                          , mcmc_chain[, "theta5"]))

df = data.frame(parameter = c(rep("theta1", each = nrow(mcmc_chain)),
                              rep("theta1", nrow(data))),
                value = c(mcmc_chain[, "theta1"], data[, "theta1"]))


ggplot(data = mcmc_chain, aes(x = parameter, y = value, fill = parameter))+
  geom_boxplot(alpha=1, color="black")+
  scale_fill_manual(values=c("#f2cd0d", "#0d32f2", "#f60997", "#09F668", "#ba7045"))+
  scale_alpha_manual(values=c(1, 0.8, 0.6, 0.4, 0.2))+
  xlab("Parameters")+
  ylab("Values")+  
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(plot.background = element_rect(fill = "white"))

ggplot(data = df, aes(x = parameter, y = value, fill = parameter))+
  geom_boxplot(alpha=1, color="black")+
  scale_fill_manual(values=c("#f2cd0d", "#0d32f2", "#f60997", "#09F668", "#ba7045"))+
  scale_alpha_manual(values=c(1, 2.5, 0.5, 0.3, 0.1))+
  xlab("Parameters")+
  ylab("Values")+
  theme(plot.background = element_rect(fill = "white"))

  legend(c("thetan1", "thetan1", "thetan1", "thetan1", "thetan5"))





boxplot(mcmc_chain[, c(1, 3)], col = col1)

# Second subplot
col2 <- c(rgb(0, 1, 0), rgb(0, 0.5, 0))
boxplot(mcmc_chain[, c(2, 4)], col = col2)






















par(mfrow = c(1, 2))
boxplot(mcmc_chain[, 1], col = rgb(242,205,13, max = 255, alpha = 114))
boxplot(data[, 1], col = rgb(242,205,13, max = 255))


par(mfrow = c(1, 2))
boxplot(mcmc_chain[, 1], col = rgb(13,50,242, max = 255, alpha = 114))
boxplot(data[, 1], col = rgb(13,50,242, max = 255))


par(mfrow = c(1, 2))
boxplot(mcmc_chain[, 1], col = rgb(246,9,151, max = 255, alpha = 114))
boxplot(data[, 1], col = rgb(246,9,151, max = 255))


par(mfrow = c(1, 2))
boxplot(mcmc_chain[, 1], col = rgb(9,246,104, max = 255, alpha = 114))
boxplot(data[, 1], col = rgb(9,246,104, max = 255))


par(mfrow = c(1, 2))
boxplot(mcmc_chain[, 1], col = rgb(186,112,69, max = 255, alpha = 114))
boxplot(data[, 1], col = rgb(186,112,69, max = 255))










df = data.frame(parameter = rep(c("theta1", "theta1")
                                , each = nrow(mcmc_chain)), 
                value = c(mcmc_chain[, "theta1"], mcmc_chain[, "theta2"]
                          , mcmc_chain[, "theta3"], mcmc_chain[, "theta4"]
                          , mcmc_chain[, "theta5"]))
ggplot(data = df, aes(x = parameter, y = value, fill = parameter))+
  geom_boxplot(alpha=0.5, color="black")+
  scale_fill_manual(values=c("red", "red", "red", "red", "red"))+
  scale_alpha_manual(values=c(1, 2.5, 0.5, 0.3, 0.1))+
  xlab("Parameters")+
  ylab("Values")









library(plotly)

theta_1 <- 42.29
theta_2 <- 7.16
beta_1 <- 9.93
beta_2 <- 2.37
delta <- 0.7

F <- function(x, y) {
  1 - exp(-((x/theta_1)^(beta_1/delta) + (y/theta_2)^(beta_2/delta))^delta)
}

# Generate a grid of x and y values to plot the function
x <- seq(0, 100, length.out = 1000)
y <- seq(0, 30, length.out = 1000)
xy_grid <- expand.grid(x = x, y = y)
z <- with(xy_grid, F(x, y))
xyz <- cbind(xy_grid, z = z)

# Plot the 3D surface
p <- plot_ly(xyz, x = x, y = y, z = z, type = "surface")
p <- layout(p, scene = list(xaxis = list(title = "X"), 
                            yaxis = list(title = "Y"), 
                            zaxis = list(title = "F(X, Y)")))
p




library(ggplot2)

# Create sample data
set.seed(123)
n <- 1000
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)
y1 <- x1 + x2 + rnorm(n, mean = 0, sd = 0.5)
y2 <- 2*x1 + x2 + rnorm(n, mean = 0, sd = 0.5)
df <- data.frame(x1 = c(x1, x1), x2 = c(x2, x2), y = c(y1, y2), 
                 group = rep(c("group1", "group2"), each = n))

# Create scatter plot with two regression lines
ggplot(df, aes(x1, x2, color = group)) + 
  geom_point(aes(size = 2)) +
  scale_color_manual(values = c("blue", "red")) +
  geom_smooth(data = subset(df, group == "group1"), method = "lm", 
              se = FALSE, color = "blue", formula = y ~ x1 + x2) +
  geom_smooth(data = subset(df, group == "group2"), method = "lm", 
              se = FALSE, color = "red", formula = y ~ x1 + x2) +
  labs(title = "Bivariate Distribution Parameter Estimates",
       x = "Parameter 1", y = "Parameter 2") +
  theme_bw(base_size = 14)




mcmc_trace(fit, pars = c("theta1", "theta2", "theta3"), inc_warmup = TRUE)
mcmc_areas(fit, pars = c("theta1", "theta2", "theta3"), prob = 0.8, inc_warmup = TRUE)
mcmc_acf(fit, pars = c("theta1", "theta2", "theta3"), lag.max = 20, inc_warmup = TRUE)
mcmc_n_eff(fit, pars = c("theta1", "theta2", "theta3"), inc_warmup = TRUE)


bayesplot::color_scheme_set("brightblue")
mcmc_acf(fit, pars = c("theta1", "theta2", "theta3", "theta4", "theta5"), lag.max = 20, inc_warmup = TRUE) +
  theme(panel.grid = element_blank(),
        text = element_text(family = "times", size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")) +
  labs(title = "", x = "", y = "") +
  scale_color_manual(values = c("blue", "red", "green")) +
  theme(plot.title = element_text(size = 12, color = "black"))

