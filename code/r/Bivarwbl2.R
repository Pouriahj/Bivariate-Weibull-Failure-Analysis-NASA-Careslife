library(rstan);
library(readxl);
library(plotly);
library(bayesplot);
library(ggplot2);

dat_mod_cp = stan_model("bivarweibullt1.stan")

dat = read_excel("C:/Users/PouriaHj/Documents/newdata.xlsx")

d1 = dat[ , 1];
d2 = dat[ , 2];

tp1 = 30;

fit_cp = sampling(dat_mod_cp, data = list(xy = dat, LENGTH = tp1), 
                  control = list(adapt_delta = 0.9))

posterior_cp = as.array(fit_cp)

lp_cp = log_posterior(fit_cp)
head(lp_cp)

np_cp = nuts_params(fit_cp)
head(np_cp)

# for the second model
lp_ncp = log_posterior(fit_cp)
np_ncp = nuts_params(fit_cp)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp)

color_scheme_set("pink")
mcmc_pairs(posterior_cp, np = np_cp, pars = c("theta1", "theta2", "theta3", "theta4","theta5"),
           off_diag_args = list(size = 0.75))

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_cp, pars = "theta1", np = np_cp) + 
  xlab("Post-warmup iteration")

mcmc_nuts_divergence(np_cp, lp_cp, chain = 4)

color_scheme_set("darkgray")
pp_check(fit_cp, type = "overlaid")

bayesplot::color_scheme_set("brightblue")
# bayesplot::mcmc_areas(fit, prob = 0.5, prob_outer = 0.9)
mcmc_areas(x = fit_cp, 
           pars = c("theta1", "theta3"),
           point_est = "mean"
)+ theme(
  legend.key = element_rect(linewidth = 1, colour = "red")
)
