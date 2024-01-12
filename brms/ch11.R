library(tidyverse)
library(ggthemes)
library(GGally)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggmcmc)
library(patchwork)
library(here)

library(wesanderson)
library(flextable)

here::set_here("/Users/lsh1510922/Documents/github/stat_rethinking_2023/brms")
fs::dir_create("fits")

## some theme nonsense
theme_set(
	theme_default() + 
		theme_tufte() 
)

# Binomial ----------------------------------------------------------------
data(chimpanzees, package = "rethinking")
d <- chimpanzees
rm(chimpanzees)

d %>% 
	distinct(prosoc_left, condition) %>% 
	mutate(description = str_c("Two food items on ", c("right and no partner",
																										 "left and no partner",
																										 "right and partner present",
																										 "left and partner present"))) %>% 
	flextable() %>% 
	width(width = c(1, 1, 4))


d <- d %>% 
	mutate(treatment = factor(1 + prosoc_left + 2 * condition)) %>% 
	# this will come in handy, later
	mutate(labels = factor(treatment,
												 levels = 1:4,
												 labels = c("r/n", "l/n", "r/p", "l/p")))

d %>% 
	count(condition, treatment, prosoc_left)

## BIG sd on the prior (log odds scale so bad idea)
b11.1 <-
	brm(data = d, 
			family = binomial,
			pulled_left | trials(1) ~ 1,
			prior(normal(0, 10), class = Intercept),
			seed = 11,
			sample_prior = T,
			file = "fits/b11.01")

b11.1b <-
	brm(data = d, 
			family = binomial,
			pulled_left | trials(1) ~ 1,
			prior(normal(0, 1.5), class = Intercept),
			seed = 11,
			sample_prior = T,
			file = "fits/b11.01b")

# wrangle
bind_rows(prior_draws(b11.1),
					prior_draws(b11.1b)) %>% 
	mutate(p = inv_logit_scaled(Intercept),
				 w = factor(rep(c(10, 1.5), each = n() / 2),
				 					 levels = c(10, 1.5))) %>% 
	# plot
	ggplot(aes(x = p, fill = w)) +
	geom_density(linewidth = 0, alpha = 3/4, adjust = 0.1) +
	scale_fill_manual(expression(italic(w)), values = c("dodgerblue", "orangered")) +
	scale_y_continuous(NULL, breaks = NULL) +
	labs(title = expression(alpha%~%Normal(0*", "*italic(w))),
			 x = "prior prob pull left")

## Doing it Index style (not with dummy vars) - ????
# gives us an overall intercept (a) and individual intercepts for each level of 
# treatment 
# w = 10
b11.2 <- 
	brm(data = d, 
			family = binomial,
			bf(pulled_left | trials(1) ~ a + b,
				 a ~ 1, 
				 b ~ 0 + treatment,
				 nl = TRUE),
			prior = c(prior(normal(0, 1.5), nlpar = a),
								prior(normal(0, 10), nlpar = b, coef = treatment1),
								prior(normal(0, 10), nlpar = b, coef = treatment2),
								prior(normal(0, 10), nlpar = b, coef = treatment3),
								prior(normal(0, 10), nlpar = b, coef = treatment4)),
			iter = 2000, warmup = 1000, chains = 4, cores = 4,
			seed = 11,
			sample_prior = T,
			file = "fits/b11.02")

# w = 0.5
b11.3 <- 
	brm(data = d, 
			family = binomial,
			bf(pulled_left | trials(1) ~ a + b,
				 a ~ 1, 
				 b ~ 0 + treatment,
				 nl = TRUE),
			prior = c(prior(normal(0, 1.5), nlpar = a),
								prior(normal(0, 0.5), nlpar = b, coef = treatment1),
								prior(normal(0, 0.5), nlpar = b, coef = treatment2),
								prior(normal(0, 0.5), nlpar = b, coef = treatment3),
								prior(normal(0, 0.5), nlpar = b, coef = treatment4)),
			iter = 2000, warmup = 1000, chains = 4, cores = 4,
			seed = 11,
			sample_prior = T,
			file = "fits/b11.03")


# Dummy approach - the four "treatments" correspond to 00, 01, 10, 11 combo of food on left/partner
# alternative is to include these binaries and their interaction in a model 
b10.2 <-
	brm(data = d, 
			family = binomial,
			formula = pulled_left | trials(1) ~ 1 + prosoc_left + condition + condition:prosoc_left,
			prior = c(prior(normal(0, 10), class = Intercept),
								prior(normal(0, 10), class = b)),
			seed = 10,
			file = "fits/b10.02")

## are these equivalent?? Should be, but how to work it out? 
# Get the fitted values from the interaction model
# create all possibloe values of condition and prosoc_left
fitdf <- expand.grid(
	prosoc_left = 0:1,
	condition = 0:1
) 

# get fitted values for the 4 possible combos
dummy_estimates <- fitted(b10.2, 
							newdata = fitdf) %>%
	as_tibble() %>% 
	bind_cols(fitdf) %>% 
	distinct(Estimate, Q2.5, Q97.5, condition, prosoc_left) %>% 
	mutate(x_axis = str_c(prosoc_left, condition, sep = "/")) %>%
	mutate(x_axis = factor(x_axis, levels = c("0/0", "1/0", "0/1", "1/1"))) %>% 
	rename(pulled_left = Estimate) %>% 
	arrange(condition, prosoc_left)
dummy_estimates

# for the index approach. Get all the draws
# Mean is the intercept + beta for that treatment. and inv_logit_scale it
index_estimates <- as_draws_df(b11.3) %>% 
	pivot_longer(b_b_treatment1:b_b_treatment4) %>% 
	mutate(treatment = str_remove(name, "b_b_treatment"),
				 mean      = inv_logit_scaled(b_a_Intercept + value)) %>%
	group_by(treatment) %>% 
	mean_qi(mean)

index_estimates %>%
	mutate(condition = c(0,0,1,1),
				 prosoc_left = c(0,1,0,1),
				 x_axis = c("0/0", "1/0", "0/1", "1/1")) %>%
	dplyr::select(pulled_left = mean,
								Q2.5 = .lower,
								Q97.5 = .upper,
								prosoc_left,
								condition,
								x_axis) %>% 
	bind_rows(dummy_estimates) %>% 
	janitor::clean_names() %>% 
	mutate(model = c(rep("index", 4), rep("dummy", 4))) %>% 
	#
	ggplot(aes(y = x_axis, x = pulled_left, xmin = q2_5, xmax = q97_5, colour = model)) +
	geom_point(position = position_dodge(width = 0.25)) +
	geom_linerange(position = position_dodge(width = 0.25)) +
	theme_ggdist()
	
plot(conditional_effects(b11.3), points = TRUE)

# wrangle
prior <-
 	bind_rows(prior_draws(b11.2),
						prior_draws(b11.3)) %>% 
	mutate(w  = factor(rep(c(10, 0.5), each = n() / 2),
										 levels = c(10, 0.5)),
				 p1 = inv_logit_scaled(b_a + b_b_treatment1),
				 p2 = inv_logit_scaled(b_a + b_b_treatment2)) %>% 
	mutate(diff = abs(p1 - p2)) 

# plot
prior %>% 
	ggplot(aes(x = diff, fill = w)) +
	geom_density(linewidth = 0, alpha = 3/4, adjust = 0.1) +
	scale_fill_manual(expression(italic(w)), values = c("dodgerblue", "orangered")) +
	scale_y_continuous(NULL, breaks = NULL) +
	labs(title = expression(alpha%~%Normal(0*", "*italic(w))),
			 x = "prior diff between treatments")



# Aggregated binomial and postcheck ---------------------------------------
data(UCBadmit, package = "rethinking")
d <- UCBadmit
rm(UCBadmit)

d <- 
	d %>%  
	mutate(gid  = factor(applicant.gender, levels = c("male", "female")),
				 case = factor(1:n()))
b11.7 <-
	brm(data = d, 
			family = binomial,
			admit | trials(applications) ~ 0 + gid,
			prior(normal(0, 1.5), class = b),
			iter = 2000, warmup = 1000, cores = 4, chains = 4,
			seed = 11,
			file = "fits/b11.07") 


as_draws_df(b11.7) %>% 
	mutate(diff_a = b_gidmale - b_gidfemale,
				 diff_p = inv_logit_scaled(b_gidmale) - inv_logit_scaled(b_gidfemale)) %>% 
	pivot_longer(contains("diff")) %>% 
	group_by(name) %>% 
	mean_qi(value, .width = .89)


# poisson regression  -----------------------------------------------------
data(Kline, package = "rethinking")
d <- Kline
rm(Kline)

d <-
	d %>%
	mutate(log_pop_std = (log(population) - mean(log(population))) / sd(log(population)),
				 cid         = contact)

# but what to do about priors? 
# a normal prior for the rate, on the outcome scale has a lognormal prior: 

# sd = 10
curve(dlnorm(x,0,10),from=0,to=100,n=200)
# big spike at 0 and a huge tail. How big? 
rnorm(1e4, 0, 10) %>% exp() %>% mean() # BIG 

curve(dlnorm(x,0,10),from=0,to=100,n=200)
curve(dlnorm(x,0,2),from=0,to=100,n=200, col = 2, add = T)
curve(dlnorm(x,0, 0.2),from=0,to=100,n=200, col = 3, add = T)
curve(dlnorm(x,0,0.5),from=0,to=100,n=200, col = 4, add = T)
curve(dlnorm(x,3,0.5),from=0,to=100,n=200, col = 5, add = T)
curve(dlnorm(x,4,0.5),from=0,to=100,n=200, col = 6, add = T)

# so flat priors are a terror! and are best avoided
curve(dlnorm(x,3,0.5),from=0,to=100,n=200)
#with mean: 
exp(3 + ((0.5^2) / 2))

# cool. Now we need a prior for the beta
N <- 100 
a <- rnorm( N, 3, 0.5 ) ## fix the intercept prior at our nice sensible value
b <- rnorm( N, 0, 10 ) ## assume a 'conventional' flat prior
plot( NULL, xlim = c(-2,2), ylim = c(0,100))
for(i in 1:N ){
	curve( exp( a[i] + b[i]*x), add = TRUE, col = rethinking::grau() )
}
# pivots around zero because standardised pop size, 
# so either catstrophic decline just after mean 
# of massive take off just before

# in other words: BAD

b <- rnorm( N, 0, 0.2)
plot( NULL, xlim = c(-2,2), ylim = c(0,100))
for(i in 1:N ){
	curve( exp( a[i] + b[i]*x), add = TRUE, col = rethinking::grau() )
}
## flatter relationships are now possible

## but let's look back on the log(pop) and raw pop scale (because standrdised are confusing)
set.seed(11)
n <- 100
prior <-
	tibble(i = 1:n,
				 a = rnorm(n, mean = 3, sd = 0.5),
				 b = rnorm(n, mean = 0, sd = 0.2)) %>% 
	expand_grid(x = seq(from = log(100), to = log(200000), length.out = 100))

# left
p1 <-
	prior %>% 
	ggplot(aes(x = x, y = exp(a + b * x), group = i)) +
	geom_line(linewidth = 1/4, alpha = 2/3,
						color = wes_palette("Moonrise2")[4]) +
	labs(subtitle = expression(beta%~%Normal(0*', '*0.2)),
			 x = "log population",
			 y = "total tools") +
	coord_cartesian(xlim = c(log(100), log(200000)),
									ylim = c(0, 500))
# right
p2 <-
	prior %>% 
	ggplot(aes(x = exp(x), y = exp(a + b * x), group = i)) +
	geom_line(linewidth = 1/4, alpha = 2/3,
						color = wes_palette("Moonrise2")[4]) +
	labs(subtitle = expression(beta%~%Normal(0*', '*0.2)),
			 x = "population",
			 y = "total tools") +
	coord_cartesian(xlim = c(100, 200000),
									ylim = c(0, 500))

# combine
p1 | p2


## Fit some models with our new sexy priors
# intercept only
b11.9 <-
	brm(data = d, 
			family = poisson,
			total_tools ~ 1,
			prior(normal(3, 0.5), class = Intercept),
			iter = 2000, warmup = 1000, chains = 4, cores = 4,
			seed = 11,
			file = "fits/b11.09") 

# interaction model
b11.10 <-
	brm(data = d, 
			family = poisson,
			bf(total_tools ~ a + b * log_pop_std,
				 a + b ~ 0 + cid,
				 nl = TRUE),
			prior = c(prior(normal(3, 0.5), nlpar = a),
								prior(normal(0, 0.2), nlpar = b)),
			iter = 2000, warmup = 1000, chains = 4, cores = 4,
			seed = 11,
			file = "fits/b11.10") 
