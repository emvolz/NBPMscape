# Calculating a latent period distribution using the following distributions as
# computed by
# He, X., Lau, E.H.Y., Wu, P. et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat Med 26, 672–675 (2020). https://doi.org/10.1038/s41591-020-0869-5
# The code for computing the following distributions is available at
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0869-5/MediaObjects/41591_2020_869_MOESM4_ESM.r
# 1) serial interval (symptom onset to symptom onset period)
# 2) infectiousness
# 3) infection to symptom onset

# After running He et al's code
# ~/mSCAPE/1_epidemic_modelling/parameters/latent_period/He_2020_Nat_Med/41591_2020_869_MOESM4_ESM_KDedited.r
# we have

# Create x-range using the min and max of the three distributions generated in He et al's code
x <- seq(-10, 20, 0.1)

# Serial interval distribution parameters from He et al's (2020) code
min.serial <- -7
ser.par <- c(8.123758,0.6361684)
# Function from He et al's code
# dgamma.shift <- function(x,min.serial,gpar1,gpar2) dgamma(x-min.serial,gpar1,gpar2)
serial_interval_shifted_gamma <- dgamma.shift(x+0.5, min.serial, ser.par[1], ser.par[2])
plot(x = x, y = serial_interval_shifted_gamma)
lines(x = x, y = serial_interval_shifted_gamma, col="red") # overlays plot in He et al's code

# Inferred distribution of infectiousness - parameters from He et al's (2020) code
inf.par <- c(20.51651, 1.592124, 12.27248)
infectiousness_gamma <- dgamma(x+inf.par[3], shape = inf.par[1], rate = inf.par[2])
plot(x = x, y = infectiousness_gamma)
lines(x = x, y = infectiousness_gamma, col = "red") # overlays plot in He et al's code

# Incubation period - parameters from He et al's (2020) code and
# from Li et al, N Engl J Med 2020;382:1199-207. DOI: 10.1056/NEJMoa2001316
# lognormal mean = 5.2; 95% CI = c(4.1, 7.0)
ln.par1 = 1.434065
ln.par2 = 0.6612
incubation_period_lnorm <- dlnorm(x, ln.par1, ln.par2)
plot(x = x, y = incubation_period_lnorm)
lines(x = x, y = incubation_period_lnorm, col="red") # overlays plot in He et al's code

# Load required libraries
library(fitdistrplus)
library(MASS)

# Function to generate sample data
generate_sample_data <- function(n) {
  # Generate sample data for each distribution
  #serial_interval <- rnorm(n, mean = 4, sd = 2.9)
  serial_interval <- rgamma(n, shape = ser.par[1], rate = ser.par[2]) -0.5 +min.serial
  #infectiousness <- rnorm(n, mean = -2, sd = 1.5)  # Days relative to symptom onset
  infectiousness <- rgamma(n, shape = inf.par[1], rate = inf.par[2]) -inf.par[3]
  #incubation <- rlnorm(n, meanlog = 1.621, sdlog = 0.418)  # Log-normal distribution
  incubation <- rlnorm(n, ln.par1, ln.par2)

  list(infectiousness = infectiousness,
       incubation = incubation,
       serial_interval = serial_interval)
}

# Function to fit a shifted log-normal distribution
fit_shifted_lognormal <- function(data) {
  # Find the minimum value to determine the shift
  shift <- min(data) - 0.01  # Subtract a small value to ensure all data is positive after shifting

  # Shift the data
  shifted_data <- data - shift

  # Fit log-normal to shifted data
  fit <- fitdist(shifted_data, "lnorm")

  # Return fit and shift
  list(fit = fit, shift = shift)
}


# Generate sample data
n <- 1000000

synthetic_data <- generate_sample_data(n)

par(mfrow=c(4,1), mar=c(5,5,1,1))
hist( synthetic_data$incubation, probability = TRUE, breaks = 1000, xlim = c(-10,20))
hist( synthetic_data$infectiousness, probability = TRUE, breaks = 1000, xlim = c(-10,20))
hist( synthetic_data$serial_interval, probability = TRUE, breaks = 1000, xlim = c(-10,20))

# Estimate start of infectiousness (t_start)
t_start <- quantile( synthetic_data$infectiousness, 0.025)  # Using 2.5th percentile as start of infectiousness

# Calculate latent period for each incubation period
latent_period <- synthetic_data$incubation - abs(t_start)
hist(latent_period, probability = TRUE, breaks = 1000)

##############
# KD
# Try fitting to positive values of latent period
latent_period_positive <- subset( latent_period, latent_period >= 0 )
mean( latent_period_positive ) # [1] 3.591764
median( latent_period_positive ) # [1] 2.352222

fit_latent_positive_gamma <- fitdist(latent_period_positive, "gamma")
fit_latent_positive_lnorm <- fitdist(latent_period_positive, "lnorm")
fit_latent_positive_nbinom <- fitdist(latent_period_positive, "nbinom")
fit_latent_positive_weibull <- fitdist(latent_period_positive, "weibull")

plot( fit_latent_positive_gamma )
plot( fit_latent_positive_lnorm )
plot( fit_latent_positive_nbinom )
plot( fit_latent_positive_weibull )

summary( fit_latent_positive_gamma )
summary( fit_latent_positive_lnorm )
summary( fit_latent_positive_nbinom )
summary( fit_latent_positive_weibull )


latent_positive_samples_gamma <- rgamma(n, shape = fit_latent_positive_gamma$estimate["shape"],
                                          rate = fit_latent_positive_gamma$estimate["rate"])
mean( latent_positive_samples_gamma )
median( latent_positive_samples_gamma )

# Convert rate to scale
fit_latent_positive_gamma_scale <- mean( latent_positive_samples_gamma ) * unname( fit_latent_positive_gamma$estimate["rate"] )

# Compute 95% CI
# Calculate the 95% confidence interval
#latent_gamma_ci_lower <- qgamma(0.025, shape = n * sample_mean^2 / var(data), scale = var(data) / sample_mean)
latent_gamma_ci_lower <- qgamma(0.025, shape = unname( fit_latent_positive_gamma$estimate["shape"])
                                     , rate = unname( fit_latent_positive_gamma$estimate["rate"]) )
latent_gamma_ci_upper <- qgamma(0.975, shape = unname( fit_latent_positive_gamma$estimate["shape"])
                                     , rate = unname( fit_latent_positive_gamma$estimate["rate"]) )
hist(latent_positive_samples_gamma)

# Print the results
cat("95% Confidence Interval:", ci_lower, "-", ci_upper)

latent_positive_samples_weibull <- rweibull(n, shape = fit_latent_positive_weibull$estimate["shape"],
                                               scale = fit_latent_positive_weibull$estimate["scale"])
mean( latent_positive_samples_weibull )
median( latent_positive_samples_weibull )
#Compute 95% CI
# Calculate the 95% confidence interval
#latent_gamma_ci_lower <- qgamma(0.025, shape = n * sample_mean^2 / var(data), scale = var(data) / sample_mean)
latent_weibull_ci_lower <- qweibull(0.025, shape = unname( fit_latent_positive_weibull$estimate["shape"])
                                       , scale = unname( fit_latent_positive_weibull$estimate["scale"] ) )
latent_weibull_ci_lower <- qweibull(0.975, shape = unname( fit_latent_positive_weibull$estimate["shape"])
                                     , scale = unname( fit_latent_positive_weibull$estimate["scale"] ) )

# Plot histogram of latent period samples
hist(latent_positive_samples, breaks = 30, probability = TRUE,
     main = "Estimated Latent Period Distribution",
     xlab = "Days")

# Add fitted curve to the histogram
curve(dgamma(x, shape = fit_latent_positive_gamma$estimate["shape"],
             scale = fit_latent_positive_gamma$estimate["rate"]),
              add = TRUE, col = "red", lwd = 2)

# Calculate and print mean and median of the latent period
mean_latent <- mean(latent_positive_samples)
median_latent <- median(latent_positive_samples)
cat("Mean latent period:", mean_latent, "days\n")
cat("Median latent period:", median_latent, "days\n")

#######################################################


# Fit a distribution to the latent period data
fit_latent <- fitdist(latent_period, "gamma")

# Plot the fitted distribution
plot(fit_latent)

# Print summary statistics
summary(fit_latent)

# Generate samples from the fitted latent period distribution
latent_samples <- rgamma(n, shape = fit_latent$estimate["shape"],
                         rate = fit_latent$estimate["rate"])

# Plot histogram of latent period samples
hist(latent_samples, breaks = 30, probability = TRUE,
     main = "Estimated Latent Period Distribution",
     xlab = "Days")

# Add fitted curve to the histogram
curve(dgamma(x, shape = fit_latent$estimate["shape"],
             rate = fit_latent$estimate["rate"]),
             add = TRUE, col = "red", lwd = 2)

# Calculate and print mean and median of the latent period
mean_latent <- mean(latent_samples)
median_latent <- median(latent_samples)
cat("Mean latent period:", mean_latent, "days\n")
cat("Median latent period:", median_latent, "days\n")

# Try fitting with a shifted lognormal distribution because there are some negative latent periods
fit_result <- fit_shifted_lognormal( latent_period )

# Plot the fitted distribution
x_range <- seq( min( latent_period ), max( latent_period ), length.out = 100 )
y_values <- dlnorm(x_range - fit_result$shift,
                   meanlog = fit_result$fit$estimate["meanlog"],
                   sdlog = fit_result$fit$estimate["sdlog"])

hist(latent_period, breaks = 30, probability = TRUE,
     main = "Estimated Latent Period Distribution",
     xlab = "Days")
lines(x_range, y_values, col = "red", lwd = 2)

# Print summary statistics
summary(fit_result$fit)
cat("Shift:", fit_result$shift, "\n")

# Generate samples from the fitted latent period distribution
latent_samples <- rlnorm(n, meanlog = fit_result$fit$estimate["meanlog"],
                         sdlog = fit_result$fit$estimate["sdlog"]) + fit_result$shift

# Calculate and print mean and median of the latent period
mean_latent <- mean(latent_samples)
median_latent <- median(latent_samples)
cat("Mean latent period:", mean_latent, "days\n")
cat("Median latent period:", median_latent, "days\n")

# What is the mean and median if only include positive latent periods
mean( subset( latent_samples, latent_samples >= 0 ) )
median( subset( latent_samples, latent_samples >= 0 ) )

# Calculate the proportion of negative latent periods
prop_negative <- mean(latent_samples < 0)
cat("Proportion of negative latent periods:", prop_negative, "\n")

################################################################################

# Using code and data from
# Shi Zhao, Biao Tang, Salihu S Musa, Shujuan Ma, Jiayue Zhang, Minyan Zeng, Qingping Yun, Wei Guo, Yixiang Zheng, Zuyao Yang, Zhihang Peng, Marc KC Chong, Mohammad Javanbakht, Daihai He, Maggie H. Wang,
# Estimating the generation interval and inferring the latent period of COVID-19 from the contact tracing data,
# Epidemics, Volume 36, 2021, 100482, ISSN 1755-4365, https://doi.org/10.1016/j.epidem.2021.100482.
# (https://www.sciencedirect.com/science/article/pii/S1755436521000359)

# From ~/mSCAPE/1_epidemic_modelling/parameters/latent_period/Zhao_2021_Epidemics/supp_code.R

#load(file = 'supp_dataset.RData')
load(file = "~/mSCAPE/1_epidemic_modelling/parameters/latent_period/Zhao_2021_Epidemics/supp_dataset.RData")

library(MASS)
#source('supp_functions_for_dist.R')
source('~/mSCAPE/1_epidemic_modelling/parameters/latent_period/Zhao_2021_Epidemics/supp_functions_for_dist.R')


n.sample = 999999
# GT.mean: (mean = 6.7, min = 5.4, max = 7.6)
# GT.sd: (mean = 1.8, min = 0.3, max = 3.8)
# incub.mean: (mean = 6.8, min = 6.2, max = 7.5)
# incub.sd: (mean = 4.1, min = 3.7, max = 4.8)
GT.mean = 6.7    # 95%CI: 5.4, 7.6
GT.sd = 1.8    # 95%CI: 0.3, 3.8
incub.mean = 6.8    # 95%CI: 6.2, 7.5
incub.sd = 4.1    # 95%CI: 3.7, 4.8


GT.sample = get.gamma.sample(X.mean = GT.mean, X.sd = GT.sd, n.sample = n.sample)
incub.sample = get.gamma.sample(X.mean = incub.mean, X.sd = incub.sd, n.sample = n.sample)
#
comb.random.samples = get.comb.gamma.sample(GT.mean = GT.mean
                                            , GT.sd = GT.sd
                                            , incub.mean = incub.mean
                                            , incub.sd = incub.sd
                                            , n.sample = n.sample)
##############
# KD
hist( GT.sample , xlim=c(-20,40) , breaks = 1000 )
hist( comb.random.samples , xlim=c(-20,40) , breaks = 1000)
##############


cdf.approx.fun = get.cdf.from.data(comb.random.samples)
#
# plot(c(-10:30), get.pd.array(X.array = c(-10:30), x.win = 1, approx.cdf = cdf.approx.fun), type = 'l')
#

pair.cov = cov(comb.mat$incub, comb.mat$SI)
SI.mean = mean(comb.random.samples)
SI.sd = sd(comb.random.samples)


# Simulate bivariate normal data
mu.array <- c(incub.mean, SI.mean)                                  # Mean
sigma.mat <- matrix(c(incub.sd^2, pair.cov, pair.cov, SI.sd^2), 2)  # Covariance matrix
# eigen(sigma.mat)

# Generate sample from N(mu, Sigma)
num.rv = 99999
bi.rv <- mvrnorm(num.rv, mu = mu.array, Sigma = sigma.mat)       # from Mass package
#head(bi.rv)


X.mean = incub.mean
X.sd = incub.sd
gamma.shape = (X.mean / X.sd) ^2
gamma.rate = X.mean / (X.sd ^2)
incub.rv = qgamma(p = pnorm(q = bi.rv[,1], mean = incub.mean, sd = incub.sd)
                  , shape = gamma.shape
                  , rate = gamma.rate)

#

X.array = seq(-99,99, length.out = 9999)
cdf.array = cdf.approx.fun(X.array)
SI.rv = apply(X = as.matrix(pnorm(q = bi.rv[,2], mean = SI.mean, sd = SI.sd))
              , MARGIN = 1, function(x){X.array[which.min(abs(cdf.array -x))[1]]})


sum(SI.rv < incub.rv) / num.rv
plot(SI.rv)
hist(SI.rv)
hist(incub.rv)
hist(SI.rv < incub.rv)


################################################################################
# Fitting different statistical distributions to the Zhao et al (2021)
# latent period: mean = 3.3 days (95% CI 0.2, 7.9)

# Load required packages
library(fitdistrplus)
library(ggplot2)

# Generate synthetic latent period data based on mean and CI
set.seed(123)
latent_period_mean <- 3.3
CI_95_upper <- 7.9
CI_95_lower <- 0.2
# Because 95% CI ~ [ mean - 1.96sd, mean + 1.96sd]
latent_period_sd <- ( CI_95_upper - CI_95_lower ) / (2*1.96)

latent_period_shape <- latent_period_mean^2 / latent_period_sd^2

latent_period_scale <- latent_period_sd^2 / latent_period_mean

# latent_data <- rgamma(1000, shape = 3.3, scale = (7.9 - 0.2) / 3.3)
latent_data <- rgamma(1000000, shape = latent_period_shape, scale = latent_period_scale )

# Visualize the data
ggplot(data.frame(latent_data), aes(x = latent_data)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Synthetic Latent Period Data", x = "Latent Period (days)", y = "Frequency")

# Fit gamma distribution
fit_gamma <- fitdist(latent_data, "gamma")
summary(fit_gamma)
fit_gamma_scale <- mean( latent_data ) * unname( fit_gamma$estimate["rate"] )
# 95% CI
latent_gamma_ci_lower <- qgamma(0.025, shape = unname( fit_gamma$estimate["shape"] )
                                     , scale = fit_gamma_scale )
latent_gamma_ci_upper <- qgamma(0.975, shape = unname( fit_latent_positive_gamma$estimate["shape"])
                                , scale = fit_latent_positive_gamma_scale )

# Fit Weibull distribution
fit_weibull <- fitdist(latent_data, "weibull")
summary(fit_weibull)

# Fit log-normal distribution
fit_lognormal <- fitdist(latent_data, "lnorm")
summary(fit_lognormal)

# Compare fitted distributions visually
par(mfrow = c(2, 2))
plot(fit_gamma)
plot(fit_weibull)
plot(fit_lognormal)

# Compare goodness-of-fit metrics
gof <- gofstat(list(fit_gamma, fit_weibull, fit_lognormal))
print(gof)


########################################

# Infer gamma dist from mean and 95% CI

fit_gamma <- function(mean, q025, q975) {
  # Objective function to minimize
  objective <- function(params) {
    shape <- params[1]
    scale <- params[2]

    # Calculate the difference between target and estimated values
    mean_diff <- abs(mean - shape*scale)
    q025_diff <- abs(q025 - qgamma(0.025, shape, scale))
    q975_diff <- abs(q975 - qgamma(0.975, shape, scale))

    return(mean_diff^2 + q025_diff^2 + q975_diff^2)
  }

  # Initial guess for shape and scale
  initial_guess <- c(1,1)#1.9765, 0.7002)

  # Optimize using optim
  result <- optim(par = initial_guess, fn = objective, method = "L-BFGS-B", lower = c(0.000001, 0.000001))

  # Return the fitted shape and rate parameters
  return(list(shape = result$par[1], scale = result$par[2]))
}

# Example usage
mean_value <- 3.3
q025_value <- 0.2
q975_value <- 7.9

fitted_params <- fit_gamma( mean_value, q025_value, q975_value)
print(paste("Fitted shape:", round(fitted_params$shape, 4)))
print(paste("Fitted scale:", round(fitted_params$scale, 4)))

sample <- rgamma(1000000, shape=fitted_params$shape, scale = fitted_params$scale )
mean(sample)
quantile(sample, probs=0.025)
quantile(sample, probs=0.975)
mean(sample)

# Plot
par(mfrow=c(1,1), mar=c(5,5,1,1))
hist(sample)
x_range <- seq(0,20,0.1)
plot( x_range, dgamma( x_range, shape=fitted_params$shape, scale = fitted_params$scale )
      , typ = "l", lwd = 2
      , xlab = "Days between infection and becoming infectious"
      , ylab = "Density"
      , cex.lab = 1.3
      , cex.axis = 1.3
      )
abline(h=0)
