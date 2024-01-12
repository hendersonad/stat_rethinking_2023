library(tidyverse)
library(GGally)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggmcmc)
library(patchwork)
library(here)
here::set_here("/Users/lsh1510922/Documents/github/stat_rethinking_2023/brms")
fs::dir_create("fits")

data(rugged, package = "rethinking")

d <- rugged

head(d)
glimpse(d)


d <- 
	d %>%
	mutate(log_gdp = log(rgdppc_2000))

dd <-
	d %>%
	drop_na(rgdppc_2000) %>% 
	mutate(log_gdp_std = log_gdp / mean(log_gdp),
				 rugged_std  = rugged / max(rugged),
				 cid         = ifelse(cont_africa == 1, "1", "2")) %>% 
	mutate(rugged_std_c = rugged_std - mean(rugged_std))


dat_slim <-
	dd %>%
	mutate(cid = as.integer(cid)) %>% 
	select(log_gdp_std, rugged_std, cid, rugged_std_c) %>% 
	list()

glimpse(dat_slim)

hist(dat_slim[[1]][, "log_gdp_std"])

plot_slim <- data.frame(dat_slim[[1]])
ggplot(plot_slim, aes((x = log_gdp_std))) + 
	geom_histogram(bins = 10) +
	facet_wrap(~cid)

ggplot(plot_slim, aes((x = rugged_std))) + 
	geom_histogram(bins = 10) +
	facet_wrap(~cid) +
	geom_vline(xintercept = 0.215, lty = 2, col = 4)

b9.1 <- 
	brm(data = dd, 
			family = gaussian,
			bf(log_gdp_std ~ 0 + a + b * (rugged_std - 0.215), 
				 a ~ 0 + cid, 
				 b ~ 0 + cid,
				 nl = TRUE),
			prior = c(prior(normal(1, 0.1), class = b, coef = cid1, nlpar = a),
								prior(normal(1, 0.1), class = b, coef = cid2, nlpar = a),
								prior(normal(0, 0.3), class = b, coef = cid1, nlpar = b),
								prior(normal(0, 0.3), class = b, coef = cid2, nlpar = b),
								prior(exponential(1), class = sigma)),
			chains = 1, cores = 1,
			seed = 9,
			file = "fits/b09.01")

print(b9.1) # note: Tail_ESS reflects quality of mixing in the tails of the posterior
						# whereas bulk_ESS is more relevant to the centre of the posterior

as_draws_df(b9.1) %>% 
	pivot_longer(b_a_cid1:sigma) %>% 
	group_by(name) %>% 
	tidybayes::mean_hdi(value, width = .89)


## multiple chains
parallel::detectCores() #4

b9.1b <- 
	update(b9.1, 
				 chains = 4, cores = 4,
				 seed = 9,
				 file = "fits/b09.01b")
print(b9.1b)

pairs(b9.1b, 
			off_diag_args = list(size = 1/5, alpha = 1/5))
post <- as_draws_df(b9.1b) 

post %>% 
	select(b_a_cid1:sigma) %>% 
	GGally::ggpairs()

plot(b9.1b)
plot(conditional_effects(b9.1b), points = TRUE)
pp_check(b9.1b)
loo(b9.1b, b9.1)

## Checking the chain 
plot(b9.1b, widths = c(1, 2))
bayesplot::mcmc_trace(post, 
											pars = vars(b7_a_cid1:sigma),
											facet_args = list(ncol = 3)) +
	scale_color_viridis_d() + 
	labs(title("My plot"))
# with warmups
pwarmup <- ggmcmc::ggs(b9.1b) %>%
	janitor::clean_names() %>% 
	mutate(chain = factor(chain)) %>% 
	# 
	ggplot(aes(x = iteration, y = value, colour = chain)) +
	geom_line() +
	annotate(geom = "rect", 
					 xmin = 0, xmax = 1000, ymin = -Inf, ymax = Inf,
					 fill = "gray80", alpha = 1/6, lty = 0) +
	scale_color_viridis_d() +
	facet_wrap(~parameter, scales = "free_y") +
	theme_ggdist()
pwarmup
# focus on warmpups
pwarmup +
	xlim(c(0, 50))

# ACF plots
bayesplot::mcmc_acf(post, 
										pars = vars(b_a_cid1:sigma),
										lags = 7) + 
	theme_ggdist()

# TRANK 
as_draws_df(b9.1b) %>%  
	bayesplot::mcmc_rank_overlay(pars = vars(b_a_cid1:sigma)) +
	scale_color_viridis_d() +
	theme_ggdist()


## Taming Bad chains 
b9.2 <-
	brm(data = list(y = c(-1, 1)), 
			family = gaussian,
			y ~ 1,
			prior = c(prior(normal(0, 1000), class = Intercept),
								prior(exponential(0.0001), class = sigma)),
			iter = 2000, warmup = 1000, chains = 3,
			seed = 9,
			file = "fits/b09.02")
print(b9.2)	 # not good! 

p1 <-
	as_draws_df(b9.2) %>% 
	mcmc_trace(pars = vars(b_Intercept:sigma))

p2 <-
	as_draws_df(b9.2) %>% 
	mcmc_rank_overlay(pars = vars(b_Intercept:sigma))


(
	(p1 / p2) &
		scale_color_viridis_d() &
		theme(legend.position = "none")
) +
	plot_annotation(subtitle = "These chains are not healthy") 

## so we regularise them slightly with decent priors that are less vague and stupid
b9.3 <-
	brm(data = list(y = c(-1, 1)), 
			family = gaussian,
			y ~ 1,
			prior = c(prior(normal(1, 10), class = Intercept),
								prior(exponential(1), class = sigma)),
			iter = 2000, warmup = 1000, chains = 3,
			seed = 9,
			file = "fits/b09.03")

print(b9.3)	 # better! 

post <- as_draws_df(b9.3)

# left
p1 <-
	post %>%
	select(b_Intercept) %>%
	
	ggplot(aes(x = b_Intercept)) +
	geom_density(trim = T) +
	geom_line(data = tibble(x = seq(from = -15, to = 15, length.out = 50)),
						aes(x = x, y = dnorm(x = x, mean = 0, sd = 10)),
						color = "dodgerblue", linetype = 2) +
	geom_line(data = tibble(x = seq(from = -15, to = 15, length.out = 50)),
						aes(x = x, y = dnorm(x = x, mean = 0, sd = 1000)),
						color = "red", linetype = 2) +
	xlab(expression(alpha))

# right
p2 <-
	post %>%
	select(sigma) %>%
	
	ggplot(aes(x = sigma)) +
	geom_density(trim = T) +
	geom_line(data = tibble(x = seq(from = 0, to = 10, length.out = 50)),
						aes(x = x, y = dexp(x = x, rate = 1)),
						color = "dodgerblue", linetype = 2) +
	geom_line(data = tibble(x = seq(from = 0, to = 10, length.out = 50)),
						aes(x = x, y = dexp(x = x, rate = 0.0001)),
						color = "red2", linetype = 2) +
	labs(x = expression(sigma),
			 y = NULL) +
	coord_cartesian(xlim = c(0, 10),
									ylim = c(0, 0.7))

# combine
(
	p1 + p2 &
		theme_ggdist()
) +
	plot_annotation(subtitle = "Prior (dashed) and posterior (solid) distributions for the\nmodel with weakly-informative priors, b9.3")

