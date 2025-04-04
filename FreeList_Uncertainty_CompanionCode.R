### Companion code for 'Modelling uncertainty around free-list cultural salience scores'
### Created by Dan Major-Smith & Ben Purzycki
### R Version 4.4.1


####################################################################
### Clear workspace, set working directory and install/load packages
rm(list=ls())
Sys.setenv(LANG = "en")

#setwd("C:/Users/au776065/Dropbox/Major-Smith-Purzycki Shared/Projects/FreeListUncertainty")
setwd("")

#install.packages("devtools")
#library("devtools")
#install_github('alastair-JL/AnthroTools')
library(AnthroTools)

## Note that, for brms to work, 'stan' will also need to be installed on your computer - see https://mc-stan.org/install/
#install.packages("brms")
library(brms)

#install.packages("ordbetareg")
library(ordbetareg)

#install.packages("tidyverse")
library(tidyverse)

#install.packages("ggdist")
library(ggdist)


###################################################################
### Read in Tyvan Virtues/Morality data
dat <- read.delim("FL10_Virtues.txt", comment.char="#")

# This data set contains information on participant ID ("Subj"), the order in which traits were listed ("Order"), listed traits of a 'good Tyvan' ("Good"), and coding of these traits into themes/categories ("GC"). The other variables are not relevant here, so will remove them. Also drop cases with missing data
head(dat)
summary(dat)
str(dat)

dat <- dat[, c("GC1", "GC2", "Subj", "Order", "GC")]
dat <- dat[complete.cases(dat), ]
head(dat)

# Calculate item salience
FL.sal <- CalculateSalience(dat, 
                            Subj = "Subj",
                            Order = "Order", 
                            CODE = "GC",
                            Salience = "GC.S")

# Calculate Smith's S for each code
(S <- SalienceByCode(FL.sal, 
                     Subj = "Subj",
                     CODE = "GC", 
                     Salience = "GC.S",
                     dealWithDoubles = "MAX"))

# In the methods below, we want to include 0s in the data for people who did not list said code, and hence have an item salience of 0. Can do this using the 'FreeListTable' command
FL.sal0 <- FreeListTable(FL.sal,
                         Subj = "Subj",
                         Order = "Order",
                         CODE = "GC",
                         Salience = "GC.S",
                         tableType = "MAX_SALIENCE")



### To begin with, we'll use 'hard-working' as an example trait to demonstrate these methods

# For ease, remove dash from 'hard-working' variable name
colnames(FL.sal0)[colnames(FL.sal0) == "hard-working"] <- "hardworking"

## Descriptive statistics and histogram of hardworking

# Recall that, when including item salience's of 0, the mean item salience equals Smith's S
summary(FL.sal0$hardworking)
S[S$CODE == "hard-working", ]

# Histogram
hist(FL.sal0$hardworking, main = "", xlab = "Item salience for 'hard-working'",
     cex.axis = 1.2, cex.lab = 1.5)

pdf(file = "./FL Uncertainty/hardworking_hist.pdf")
hist(FL.sal0$hardworking, main = "", xlab = "Item salience for 'hard-working'",
     cex.axis = 1.2, cex.lab = 1.5)
dev.off()



## Comparing different methods for propagating uncertainty in Smith's S estimates

# Bootstrapping
set.seed(76567)
boot_samp <- list() # To store samples from each bootstrapping iteration
S_boot <- rep(NA, 1000) # To store each Smith's S estimate
for (i in 1:1000) {
  boot_samp[[i]] <- sample(FL.sal0$hardworking, size = nrow(FL.sal0), replace = TRUE)
  S_boot[i] <- mean(boot_samp[[i]])
}

head(boot_samp)
head(S_boot)

# Smith's S 95% percentile interval
quantile(S_boot, c(0.025, 0.5, 0.975))


# Linear model (default non-informative priors - 'student_t(3, 0, 2.5)' for both intercept and sigma terms, but with a lower bound of 0 for the sigma term)
linear_mod <- brm(formula = bf(
  hardworking ~ 1),
  data = FL.sal0,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5645,
  family = gaussian())

get_prior(linear_mod)

summary(linear_mod)

# Posterior predictions from this model
post_linear <- predict(linear_mod, summary = FALSE)
head(post_linear)

# Smith's S 95% credible interval
S_post_linear <- rep(NA, nrow(post_linear))
for (i in 1:nrow(post_linear)) {
  S_post_linear[i] <- mean(post_linear[i, ])
}

quantile(S_post_linear, c(0.025, 0.5, 0.975))


# Frequentist linear model
linear_mod_freq <- lm(hardworking ~ 1, data = FL.sal0)
summary(linear_mod_freq)

# Sample from this model 1,000 times, storing Smith's S for each iteration
samples_linear_freq <- simulate(linear_mod_freq, nsim = 1000, seed = 345)
str(samples_linear_freq)

# Smith's S 95% confidence interval
S_linear_freq <- rep(NA, ncol(samples_linear_freq))
for (i in 1:ncol(samples_linear_freq)) {
  S_linear_freq[i] <- mean(samples_linear_freq[, i])
}

quantile(S_linear_freq, c(0.025, 0.5, 0.975))



# ZOIB model (default non-informative priors; 'student_t(3, 0, 2.5)' for mean and phi, 'logistic(0, 1)' for zoi and coi)
zoib_mod <- brm(formula = bf(
  hardworking ~ 1, # Mean of beta distribution if between 0 and 1 (on logit scale)
  phi ~ 1, # Width of beta distribution if between 0 and 1 (on log scale)
  zoi ~ 1, # Probability of item salience being either 0 or 1 or not 0 or 1 (on logit scale)
  coi ~ 1), # Conditional inflation of 1s (i.e., probability of item salience being a 1, if either 0 or 1; on logit scale)
  data = FL.sal0,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5323,
  family = zero_one_inflated_beta())

get_prior(zoib_mod)

summary(zoib_mod)

# Posterior predictions from this model
post_zoib <- predict(zoib_mod, summary = FALSE)
head(post_zoib)

# Smith's S 95% credible interval
S_post_zoib <- rep(NA, nrow(post_zoib))
for (i in 1:nrow(post_zoib)) {
  S_post_zoib[i] <- mean(post_zoib[i, ])
}

quantile(S_post_zoib, c(0.025, 0.5, 0.975))


# Ordered Beta (default non-informative priors, where possible; 'student_t(3, 0, 2.5)' for mean and phi terms [note that if include a phi intercept model, have to manually specify a prior], flat for cut-points)
ordBeta_mod <- ordbetareg(
  formula = bf(hardworking ~ 1, # Mean of beta distribution if between 0 and 1 (on logit scale)
               phi ~ 1), # Width of beta distribution if between 0 and 1 (on log scale)
  data = FL.sal0,
  true_bounds = c(0, 1),
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5543)

get_prior(ordBeta_mod)

summary(ordBeta_mod)

# Posterior predictions from this model
post_ordBeta <- predict(ordBeta_mod, summary = FALSE)
head(post_ordBeta)

# Smith's S 95% credible interval
S_post_ordBeta <- rep(NA, nrow(post_ordBeta))
for (i in 1:nrow(post_ordBeta)) {
  S_post_ordBeta[i] <- mean(post_ordBeta[i, ])
}

quantile(S_post_ordBeta, c(0.025, 0.5, 0.975))



## Plot 50 predicted densities on top of observed density (plus bootstrapping)

pdf(file = "./FL Uncertainty/bootLinearZOIBOrdered.pdf", height = 12, width = 5)

# Four panel plot
par(mfrow = c(4, 1), mai = c(0.3, 0.5, 0.3, 0.5))

# Bootstrap results
plot(density(FL.sal0$hardworking), ylim = c(0, 4), xlim = c(0, 1), col = "blue", lwd = 4,
     main = "Bootstrap sampling", xlab = "")
for (i in 1:50) {
  lines(density(boot_samp[[i]]), col = rgb(0, 0, 0, 1/6))
}

# Linear results
plot(density(FL.sal0$hardworking), ylim = c(0, 4), xlim = c(0, 1), col = "blue", lwd = 4,
     main = "Linear model", xlab = "")
for (i in 1:50) {
  lines(density(post_linear[i, ]), col = rgb(0, 0, 0, 1/6))
}

# ZOIB results
plot(density(FL.sal0$hardworking), ylim = c(0, 4), xlim = c(0, 1), col = "blue", lwd = 4,
     main = "ZOIB", xlab = "")
for (i in 1:50) {
  lines(density(post_zoib[i, ]), col = rgb(0, 0, 0, 1/6))
}

# Ordered Beta results
plot(density(FL.sal0$hardworking), ylim = c(0, 4), xlim = c(0, 1), col = "blue", lwd = 4,
     main = "Ordered Beta", xlab = "")
for (i in 1:50) {
  lines(density(post_ordBeta[i, ]), col = rgb(0, 0, 0, 1/6))
}

dev.off()


## Plot comparing densities of bootstrap, ZOIB and ordered Beta results
par(mar = c(3, 5, 1, 2))

plot(density(S_boot), ylim = c(0, 10), xlim = c(0, 0.6), col = "blue", lwd = 4,
     main = "", xlab = "", cex.axis = 1.2, cex.lab = 1.5)
abline(v = median(S_boot), col = "blue", lwd = 2, lty = 2)
lines(density(S_post_zoib), col = "red", lwd = 4)
abline(v = median(S_post_zoib), col = "red", lwd = 2, lty = 2)
lines(density(S_post_ordBeta), col = "grey", lwd = 4)
abline(v = median(S_post_ordBeta), col = "grey", lwd = 2, lty = 2)
legend("topright", legend = c("Bootstrap", "ZOIB", "Ordered Beta"), 
       col = c("blue", "red", "grey"),
       pch = 15, pt.cex = 2)

# Save plot
pdf(file = "./FL Uncertainty/hardworking_comparison.pdf", height = 4, width = 6)

par(mar = c(3, 5, 1, 2))

plot(density(S_boot), ylim = c(0, 10), xlim = c(0, 0.6), col = "blue", lwd = 4,
     main = "", xlab = "", cex.axis = 1.2, cex.lab = 1.5)
abline(v = median(S_boot), col = "blue", lwd = 2, lty = 2)
lines(density(S_post_zoib), col = "red", lwd = 4)
abline(v = median(S_post_zoib), col = "red", lwd = 2, lty = 2)
lines(density(S_post_ordBeta), col = "grey", lwd = 4)
abline(v = median(S_post_ordBeta), col = "grey", lwd = 2, lty = 2)
legend("topright", legend = c("Bootstrap", "ZOIB", "Ordered Beta"), 
       col = c("blue", "red", "grey"),
       pch = 15, pt.cex = 2)

dev.off()



## Comparison of results if using more informative priors

# Linear model (We know that the mean item salience/Smith's S will be between 0 and 1, so use a normally-distributed prior for this centred on 0.5 with an SD of 0.5. Use an exponential distribution for the sigma term)
get_prior(linear_mod)

linear_mod_prior <- brm(formula = bf(
  hardworking ~ 1),
  data = FL.sal0,
  prior = c(prior(normal(0.5, 0.5), class = Intercept),
            prior(exponential(1), class = sigma)),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5465,
  family = gaussian())

summary(linear_mod_prior)

# Posterior predictions from this model
post_linear_prior <- predict(linear_mod_prior, summary = FALSE)
head(post_linear_prior)

# Smith's S 95% credible interval
S_post_linear_prior <- rep(NA, nrow(post_linear_prior))
for (i in 1:nrow(post_linear_prior)) {
  S_post_linear_prior[i] <- mean(post_linear_prior[i, ])
}

quantile(S_post_linear_prior, c(0.025, 0.5, 0.975))
quantile(S_post_linear, c(0.025, 0.5, 0.975))


# ZOIB model (specify tighter priors, centered on 0 with an SD of 1.5 for logit terms [beta mean, zoi and coi], and 0 and SD of 1 for log terms [phi])
get_prior(zoib_mod)

zoib_mod_prior <- brm(formula = bf(
  hardworking ~ 1, # Mean of beta distribution if between 0 and 1 (on logit scale)
  phi ~ 1, # Width of beta distribution if between 0 and 1 (on log scale)
  zoi ~ 1, # Probability of item salience being either 0 or 1 or not 0 or 1 (on logit scale)
  coi ~ 1), # Conditional inflation of 1s (i.e., probability of item salience being a 1, if either 0 or 1; on logit scale)
  data = FL.sal0,
  prior = c(prior(normal(0, 1.5), class = Intercept),
            prior(normal(0, 1), class = Intercept, dpar = phi),
            prior(normal(0, 1.5), class = Intercept, dpar = zoi),
            prior(normal(0, 1.5), class = Intercept, dpar = coi)),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 3235,
  family = zero_one_inflated_beta())

summary(zoib_mod_prior)

# Posterior predictions from this model
post_zoib_prior <- predict(zoib_mod_prior, summary = FALSE)
head(post_zoib_prior)

# Smith's S 95% credible interval
S_post_zoib_prior <- rep(NA, nrow(post_zoib_prior))
for (i in 1:nrow(post_zoib_prior)) {
  S_post_zoib_prior[i] <- mean(post_zoib_prior[i, ])
}

quantile(S_post_zoib_prior, c(0.025, 0.5, 0.975))
quantile(S_post_zoib, c(0.025, 0.5, 0.975))


# Ordered Beta (specify tighter priors, centered on 0 with an SD of 1.5 for logit terms [beta mean, and the two cut-points], and 0 and SD of 1 for log terms [phi])
get_prior(ordBeta_mod)

ordBeta_mod_prior <- ordbetareg(
  formula = bf(hardworking ~ 1, # Mean of beta distribution if between 0 and 1 (on logit scale)
               phi ~ 1), # Width of beta distribution if between 0 and 1 (on log scale)
  data = FL.sal0,
  true_bounds = c(0, 1),
  manual_prior = c(prior(normal(0, 1.5), class = Intercept),
    prior(normal(0, 1), class = Intercept, dpar = phi),
    prior(normal(0, 1), class = cutzero),
    prior(normal(0, 1), class = cutone)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 3455)

summary(ordBeta_mod_prior)

# Posterior predictions from this model
post_ordBeta_prior <- predict(ordBeta_mod_prior, summary = FALSE)
head(post_ordBeta_prior)

# Smith's S 95% credible interval
S_post_ordBeta_prior <- rep(NA, nrow(post_ordBeta_prior))
for (i in 1:nrow(post_ordBeta_prior)) {
  S_post_ordBeta_prior[i] <- mean(post_ordBeta_prior[i, ])
}

quantile(S_post_ordBeta_prior, c(0.025, 0.5, 0.975))
quantile(S_post_ordBeta, c(0.025, 0.5, 0.975))



### Comparing bootstrap, linear, ZOIB and ordered Beta in a variable with nearly 0 cultural salience (to demonstrate how linear models can give nonsensical results in edge cases like this)
summary(FL.sal0)

# Use 'talented', as only 2 non-zero values (smith's S = 0.01)
hist(FL.sal0$talented)
table(FL.sal0$talented)
summary(FL.sal0$talented)

# Bootstrap
set.seed(45678)
boot_samp <- list()
S_boot <- rep(NA, 1000)
for (i in 1:1000) {
  boot_samp[[i]] <- sample(FL.sal0$talented, size = nrow(FL.sal0), replace = TRUE)
  S_boot[i] <- mean(boot_samp[[i]])
}

summary(S_boot)
quantile(S_boot, c(0.025, 0.5, 0.975))


# Linear model
linear_mod <- brm(formula = bf(
  talented ~ 1),
  data = FL.sal0,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5645,
  family = gaussian())

summary(linear_mod)

# Posterior predictions from this model
post_linear <- predict(linear_mod, summary = FALSE)
head(post_linear)

S_post_linear <- rep(NA, nrow(post_linear))
for (i in 1:nrow(post_linear)) {
  S_post_linear[i] <- mean(post_linear[i, ])
}

summary(S_post_linear)
quantile(S_post_linear, c(0.025, 0.5, 0.975))


# Frequentist linear model
linear_mod_freq <- lm(talented ~ 1, data = FL.sal0)
summary(linear_mod_freq)

# Sample from this model 1,000 times, storing Smith's S for each iteration
samples_linear_freq <- simulate(linear_mod_freq, nsim = 1000, seed = 345)
str(samples_linear_freq)

# Smith's S 95% confidence interval
S_linear_freq <- rep(NA, ncol(samples_linear_freq))
for (i in 1:ncol(samples_linear_freq)) {
  S_linear_freq[i] <- mean(samples_linear_freq[, i])
}

summary(S_linear_freq)
quantile(S_linear_freq, c(0.025, 0.5, 0.975))


# ZOIB model
zoib_mod <- brm(formula = bf(
  talented ~ 1, # Mean of beta distribution if between 0 and 1 (on logit scale)
  phi ~ 1, # Width of beta distribution if between 0 and 1 (on log scale)
  zoi ~ 1, # Probability of item salience being either 0 or 1 or not 0 or 1 (on logit scale)
  coi ~ 1), # Conditional inflation of 1s (i.e., probability of item salience being a 1, if either 0 or 1; on logit scale)
  data = FL.sal0,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5323,
  family = zero_one_inflated_beta())

summary(zoib_mod)

# Posterior predictions from this model
post_zoib <- predict(zoib_mod, summary = FALSE)
head(post_zoib)

S_post_zoib <- rep(NA, nrow(post_zoib))
for (i in 1:nrow(post_zoib)) {
  S_post_zoib[i] <- mean(post_zoib[i, ])
}

summary(S_post_zoib)
quantile(S_post_zoib, c(0.025, 0.5, 0.975))


# Ordered Beta (quite a few divergence warnings if don't fix the 'cutone' parameter, because so are no '1s' to model - For more on this, see the code section below)
ordBeta_mod <- ordbetareg(
  formula = bf(talented ~ 1, 
               phi ~ 1, 
               cutone = 10), 
  data = FL.sal0,
  true_bounds = c(0, 1),
  manual_prior = c(prior(normal(0, 1), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5442)

summary(ordBeta_mod)

# Posterior predictions from this model
post_ordBeta <- predict(ordBeta_mod, summary = FALSE)
head(post_ordBeta)

S_post_ordBeta <- rep(NA, nrow(post_ordBeta))
for (i in 1:nrow(post_ordBeta)) {
  S_post_ordBeta[i] <- mean(post_ordBeta[i, ])
}

summary(S_post_ordBeta)
quantile(S_post_ordBeta, c(0.025, 0.5, 0.975))



#### Comparing models with no 0s, no 1s, or neither 0s or 1s.

## Edit data first, recoding the 'hard-working' variable
FL.sal0.work <- FL.sal0[, c("Subject", "hardworking")]
FL.sal0.work$hardworking_no1s <- ifelse(FL.sal0.work$hardworking == 1, 
                                        runif(n = nrow(FL.sal0.work[FL.sal0.work$hardworking == 1, ]), 
                                                             min = 0, max = 0.9),
                                        FL.sal0.work$hardworking)
FL.sal0.work$hardworking_no0s <- ifelse(FL.sal0.work$hardworking == 0, 
                                        runif(n = nrow(FL.sal0.work[FL.sal0.work$hardworking == 0, ]), 
                                              min = 0.1, max = 1),
                                        FL.sal0.work$hardworking)
FL.sal0.work$hardworking_no0s1s <- ifelse(FL.sal0.work$hardworking == 0 | FL.sal0.work$hardworking == 1, 
                                        runif(n = nrow(FL.sal0.work[FL.sal0.work$hardworking == 0 | 
                                                                      FL.sal0.work$hardworking_no1s == 1, ]), 
                                              min = 0.1, max = 0.9),
                                        FL.sal0.work$hardworking)
head(FL.sal0.work)
summary(FL.sal0.work)

hist(FL.sal0.work$hardworking)
hist(FL.sal0.work$hardworking_no1s)
hist(FL.sal0.work$hardworking_no0s)
hist(FL.sal0.work$hardworking_no0s1s)


### ZOIB models

## No 1s

# First, trying to estimate all terms (this works, but is a lot of uncertainty in the coi/one-inflated estimate)
zoib_mod_no1s <- brm(formula = bf(
  hardworking_no1s ~ 1, 
  phi ~ 1, 
  zoi ~ 1, 
  coi ~ 1), 
  data = FL.sal0.work,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 3234,
  family = zero_one_inflated_beta())

summary(zoib_mod_no1s)


# Next, with the coi/one-inflated term fixed to 0 (meaning 'of all the 0s or 1s, there are no 1s')
zoib_mod_no1s_fixed <- brm(formula = bf(
  hardworking_no1s ~ 1, 
  phi ~ 1, 
  zoi ~ 1, 
  coi = 0), 
  data = FL.sal0.work,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 32345,
  family = zero_one_inflated_beta())

summary(zoib_mod_no1s_fixed)


## No 0s

# First, trying to estimate all terms (this works, but is a lot of uncertainty in the coi/one-inflated estimate)
zoib_mod_no0s <- brm(formula = bf(
  hardworking_no0s ~ 1, 
  phi ~ 1, 
  zoi ~ 1, 
  coi ~ 1), 
  data = FL.sal0.work,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 6567,
  family = zero_one_inflated_beta())

summary(zoib_mod_no0s)

# Next, with the coi/one-inflated term fixed to 1 (meaning 'of all the 0s or 1s, there are no 0s')
zoib_mod_no0s_fixed <- brm(formula = bf(
  hardworking_no0s ~ 1, 
  phi ~ 1, 
  zoi ~ 1, 
  coi = 1), 
  data = FL.sal0.work,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 65678,
  family = zero_one_inflated_beta())

summary(zoib_mod_no0s_fixed)


## No 0s or 1s

# First, trying to estimate all terms (this works, but is a lot of uncertainty in the zoi/zero-one-inflated and coi/one-inflated estimates)
zoib_mod_no0s1s <- brm(formula = bf(
  hardworking_no0s1s ~ 1, 
  phi ~ 1, 
  zoi ~ 1, 
  coi ~ 1), 
  data = FL.sal0.work,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5647,
  family = zero_one_inflated_beta())

summary(zoib_mod_no0s1s)

# Next, with both the zoi/zero-one-inflated and coi/one-inflated terms fixed to 0 (meaning 'there are no 0s or 1s')
zoib_mod_no0s1s_fixed <- brm(formula = bf(
  hardworking_no0s1s ~ 1, 
  phi ~ 1, 
  zoi = 0, 
  coi = 0), 
  data = FL.sal0.work,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 56478,
  family = zero_one_inflated_beta())

summary(zoib_mod_no0s1s_fixed)


### Ordered Beta models

## No 1s

# First, if ignore that there are no 1s in the data (are lots divergent warnings and small effective sample sizes, as cannot model the one cut-point)
ordBeta_mod_no1s <- ordbetareg(
  formula = bf(hardworking_no1s ~ 1, 
               phi ~ 1), 
  data = FL.sal0.work,
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 1122)

summary(ordBeta_mod_no1s)

# Next, with the one cut-point fixed to 10 on the logit scale (meaning the probability of observing a 1 is miniscule, and effectively nil)
ordBeta_mod_no1s_fixed <- ordbetareg(
  formula = bf(hardworking_no1s ~ 1, 
               phi ~ 1,
               cutone = 10), 
  data = FL.sal0.work,
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 112233)

summary(ordBeta_mod_no1s_fixed)


## No 0s

# First, if ignore that there are no 0s in the data (are lots divergent warnings, large r-hat values and tiny effective sample sizes, as cannot model the zero cut-point)
ordBeta_mod_no0s <- ordbetareg(
  formula = bf(hardworking_no0s ~ 1, 
               phi ~ 1), 
  data = FL.sal0.work,
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 2233)

summary(ordBeta_mod_no0s)

# Next, with the zero cut-point fixed to -10 on the logit scale (meaning the probability of observing a 0 is miniscule, and effectively nil)
ordBeta_mod_no0s_fixed <- ordbetareg(
  formula = bf(hardworking_no0s ~ 1, 
               phi ~ 1,
               cutzero = -10), 
  data = FL.sal0.work,
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 223344)

summary(ordBeta_mod_no0s_fixed)


## No 0s or 1s

# First, if ignore that there are no 0s and no 1s in the data (are lots divergent warnings, large r-hat values and tiny effective sample sizes, as cannot model the zero or one cut-points)
ordBeta_mod_no0s1s <- ordbetareg(
  formula = bf(hardworking_no0s1s ~ 1, 
               phi ~ 1), 
  data = FL.sal0.work,
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 3344)

summary(ordBeta_mod_no0s1s)

# Next, with the zero cut-point fixed to -10 on the logit scale (meaning the probability of observing a 0 is miniscule, and effectively nil) and the one cut-point fixed to 10 on the logit scale (meaning the probability of observing a 1 is miniscule, and effectively nil)
ordBeta_mod_no0s1s_fixed <- ordbetareg(
  formula = bf(hardworking_no0s1s ~ 1, 
               phi ~ 1,
               cutzero = -10,
               cutone = 10), 
  data = FL.sal0.work,
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 334455)

summary(ordBeta_mod_no0s1s_fixed)




###############################################################################
### Describing and comparing multiple items

## Ordered Beta function for multiple variables (based on top X Smith's S salience values).
FLordBeta_SmithsS_multVar_top <- function(data, top = 8, cut_no0s = -10, cut_no1s = 10, print_model = TRUE, seed, chains = 4, iterations = 2000, warmup = 1000, cores = 4, IDs_first = TRUE) {
  
  # Make sure that 'brms' and 'ordbetareg' packages are loaded
  if (!"brms" %in% .packages()) {
    stop("Package brms is not loaded into library. Please load (or install then load, if necessary) and try again.")
  }
  if (!"ordbetareg" %in% .packages()) {
    stop("Package ordbetareg is not loaded into library. Please load (or install then load, if necessary) and try again.")
  }
  
  # Sanity check that the variables in the dataset have no missing data, or values < 0 and/or > 1. Exit if any do. Vary depending on whether first column is IDs or not
  if (IDs_first == TRUE) {
    for (i in 2:ncol(data)) { 
      var <- colnames(data)[i]
      if (any(is.na(data[[var]])) == TRUE) {
        stop("One or more specified variables have missing data - Function aborted.")
      }
      if (min(data[[var]]) < 0 | max(data[[var]]) > 1) {
        stop("One or more specified variables have values less than zero and/or greater than one - Function aborted.")
      }
    }
  } else if (IDs_first == FALSE) {
    for (i in 1:ncol(data)) { 
      var <- colnames(data)[i]
      if (any(is.na(data[[var]])) == TRUE) {
        stop("One or more specified variables have missing data - Function aborted.")
      }
      if (min(data[[var]]) < 0 | max(data[[var]]) > 1) {
        stop("One or more specified variables have values less than zero and/or greater than one - Function aborted.")
      }
    }
  }
  
  # Set seed
  set.seed(seed)
  
  # List to store results in
  res_list <- list()
  
  ## Select the top X variables in terms of Smith's S. Vary depending on whether IDs are in first column or not
  if (IDs_first == TRUE) {
    x <- colMeans(data[, -1])
  } else if (IDs_first == FALSE) {
    x <- colMeans(data)
  }
  
  x <- x[order(x, decreasing = TRUE)] # Order by Smith's S
  x <- x[1:top] # Take top X variables
  vars <- attributes(x)$names # Extract variable names
  
  # Loop over each of the top variables, first checking for 0s and 1s, then running the appropriate ordered Beta model
  for (i in 1:length(vars)) { 
    
    var <- vars[i] # Take each variable
    print(paste0("On variable ", i, ": ", var))
    
    # Make variable into a dataframe
    df_temp <- as.data.frame(cbind(v1 = data[[var]]))
    
    # See if 0s and 1s, to run appropriate ordered Beta model
    any_0s <- ifelse(min(df_temp$v1) == 0, TRUE, FALSE)
    any_1s <- ifelse(max(df_temp$v1) == 1, TRUE, FALSE)
    
    model_type <- ifelse(any_0s == TRUE & any_1s == TRUE, "Standard", 
                         ifelse(any_0s == TRUE & any_1s == FALSE, "No 1s",
                                ifelse(any_0s == FALSE & any_1s == TRUE, "No 0s", "No 0s or 1s")))
    
    ## Run appropriate ordered Beta model
    
    # Specify weakly-informative priors for phi intercept term
    priors <- c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi))
    
    if (model_type == "Standard") {
      
      print(paste0("There are both 0s and 1s for variable ", var, ". Running standard Ordered Beta."))
      
      ordBeta_mod <- ordbetareg(formula = bf(v1 ~ 1, phi ~ 1), 
                                data = df_temp,
                                true_bounds = c(0, 1),
                                manual_prior = priors,
                                phi_reg = "only",
                                chains = chains, iter = iterations, warmup = warmup, cores = cores)
      
    } else if (model_type == "No 1s") {
      
      print(paste0("There are 0s but no 1s for variable ", var, ". Running Ordered Beta with one cut-point term fixed to ", cut_no1s, " (i.e., a one value is almost impossible)."))
      
      ordBeta_mod <- ordbetareg(formula = bf(v1 ~ 1, phi ~ 1, cutone = cut_no1s), 
                                data = df_temp,
                                true_bounds = c(0, 1),
                                manual_prior = priors,
                                phi_reg = "only",
                                chains = chains, iter = iterations, warmup = warmup, cores = cores)
      
    } else if (model_type == "No 0s") {
      
      print(paste0("There are 1s but no 0s for variable ", var, ". Running Ordered Beta with zero cut-point term fixed to ", cut_no0s, " (i.e., a zero value is almost impossible)."))
      
      ordBeta_mod <- ordbetareg(formula = bf(v1 ~ 1, phi ~ 1, cutzero = cut_no0s), 
                                data = df_temp,
                                true_bounds = c(0, 1),
                                manual_prior = priors,
                                phi_reg = "only",
                                chains = chains, iter = iterations, warmup = warmup, cores = cores)
      
    } else if (model_type == "No 0s or 1s") {
      
      print(paste0("There are no 0s and no 1s for variable ", var, ". Running Ordered Beta with zero and one cut-point terms fixed to ", cut_no0s, " and ", cut_no1s, ", respectively (i.e., a standard beta model)."))
      
      ordBeta_mod <- ordbetareg(formula = bf(v1 ~ 1, phi ~ 1, cutzero = cut_no0s, cutone = cut_no1s), 
                                data = df_temp,
                                true_bounds = c(0, 1),
                                manual_prior = priors,
                                phi_reg = "only",
                                chains = chains, iter = iterations, warmup = warmup, cores = cores)
      
    }
    
    if (print_model == TRUE) {
      print(ordBeta_mod)
    }
    
    # Posterior predictions from this model
    post_ordBeta <- predict(ordBeta_mod, summary = FALSE)
    
    # Smith's S in each posterior sample
    S_post_ordBeta <- rep(NA, nrow(post_ordBeta))
    for (j in 1:nrow(post_ordBeta)) {
      S_post_ordBeta[j] <- mean(post_ordBeta[j, ])
    }
    
    # Store results in list
    res_list[[i]] <- S_post_ordBeta
    names(res_list)[i] <- var
  }
  return(res_list)
}


# Summarise top 8 items
top8 <- FLordBeta_SmithsS_multVar_top(data = FL.sal0, top = 8, seed = 54123)

# Summarise (first, create function to loop over stored results and extract credible intervals)
FL_SmithsS_summariseEstimates <- function(data, quantiles = c(0, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 1)) {
  res <- round(t(as.data.frame(lapply(data, quantile, probs = quantiles))), 2)
  return(res)
}

res <- FL_SmithsS_summariseEstimates(top8, quantiles = c(0.025, 0.5, 0.975))
res


## Flower plot summary of results

# First, standard flower plot, without uncertainty intervals
par(mar = c(0, 0, 0, 0))
FlowerPlot(S, "Good")      

# Edit row names of 'res' first, so they match the Smith's S data
S

rownames(res)[1] <- "hard-working"
rownames(res)[8] <- "good conscience"
res

# Now edit the Smith's S output to include intervals
S_uncert <- S
for (i in 1:length(rownames(res))) {
  name <- rownames(res)[i]
  S_uncert$CODE[S_uncert$CODE == name] <- paste0(name, "\n[", res[i, "2.5%"], "-", res[i, "97.5%"], "]")
}
S_uncert

# Flower plot
par(mar = c(0, 0, 0, 0))
FlowerPlot(S_uncert, "Good")

# Save as PDF
pdf(file = "./FL Uncertainty/flower_uncert.pdf", height = 6, width = 6)
par(mar = c(0, 0, 0, 0))
FlowerPlot(S_uncert, "Good")
dev.off()


## Displaying full distribution, along with CIs

# Convert list of posterior samples to data frame
res_full <- as.data.frame(top8)
names(res_full)[1] <- "hard working"
names(res_full)[8] <- "good conscience"
res_full

# Make long format
res_full <- res_full %>%
  pivot_longer(cols = everything(), names_to = "item", values_to = "S") %>%
  mutate(item = factor(item, levels = c("hard working", "kind", "helpful", "modest",
                                        "respectful", "honest", "intelligent", "good conscience"))) %>%
  arrange(item)
res_full

# Plot of results
(p_dist <- ggplot(res_full, aes(x = S, y = forcats::fct_rev(item), fill = forcats::fct_rev(item))) +
    stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_qi") +
    scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.9) +
    guides(fill = "none") +
    labs(x = expression(italic(hat(S)) ~ "cultural salience"), y = "Item") +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.line = element_line(colour = "black"),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()))

# Save as PDF
pdf(file = "./FL Uncertainty/top8_dist.pdf", height = 8, width = 8)
p_dist
dev.off()


### Comparison between items

## First on the absolute difference scale 

# Hard working vs kind
S_post_ordBeta_hardVsKind <- rep(NA, length(top8$hardworking))
for (i in 1:length(top8$hardworking)) {
  S_post_ordBeta_hardVsKind[i] <- top8$kind[i] - top8$hardworking[i]
}

summary(S_post_ordBeta_hardVsKind)
quantile(S_post_ordBeta_hardVsKind, c(0.025, 0.5, 0.975))
hist(S_post_ordBeta_hardVsKind)

# Hard working vs good conscience
S_post_ordBeta_hardVsconscience <- rep(NA, length(top8$hardworking))
for (i in 1:length(top8$hardworking)) {
  S_post_ordBeta_hardVsconscience[i] <- top8$`good conscience`[i] - top8$hardworking[i]
}

summary(S_post_ordBeta_hardVsconscience)
quantile(S_post_ordBeta_hardVsconscience, c(0.025, 0.5, 0.975))
hist(S_post_ordBeta_hardVsconscience)


## And on the ratio difference scale

# Hard working vs kind
S_post_ordBeta_hardVsKind_per <- rep(NA, length(top8$hardworking))
for (i in 1:length(top8$hardworking)) {
  S_post_ordBeta_hardVsKind_per[i] <- top8$hardworking[i] / top8$kind[i]
}

summary(S_post_ordBeta_hardVsKind_per)
quantile(S_post_ordBeta_hardVsKind_per, c(0.025, 0.5, 0.975))
hist(S_post_ordBeta_hardVsKind_per)

# Hard working vs good conscience
S_post_ordBeta_hardVsconscience_per <- rep(NA, length(top8$hardworking))
for (i in 1:length(top8$hardworking)) {
  S_post_ordBeta_hardVsconscience_per[i] <- top8$hardworking[i] / top8$`good conscience`[i]
}

summary(S_post_ordBeta_hardVsconscience_per)
quantile(S_post_ordBeta_hardVsconscience_per, c(0.025, 0.5, 0.975))
hist(S_post_ordBeta_hardVsconscience_per)



## Write a function to automate this and present results in a distribution plot - Using list from function above
FL_SmithsS_generateContrasts <- function(data, contrast = "absolute_diff") {
  
  # List to store results in
  res_list <- list()
  
  # Keep a counter of number of contrasts
  counter <- 0
  
  # Loop over each of the variables
  for (i in 1:length(names(data))) { 
    
    # Take each variable combination
    var_i <- names(data)[i] 
    for (j in 1:length(names(data))) {
      var_j <- names(data)[j] 
      
      # Skip if same pairing
      if (var_i == var_j) next
      
      print(paste0("On variables ", i, ": ", var_i, " and ", j, ": ", var_j))
      
      # Vector to store results in
      diff <- rep(NA, length(data[[var_i]]))
      
      # If want results in terms of difference in values vs difference in ratio of values
      for (k in 1:length(data[[var_i]])) {
        if (contrast == "absolute_diff") {
          diff[k] <- data[[var_j]][k] - data[[var_i]][k]
        } else if (contrast == "ratio_diff") {
          diff[k] <- data[[var_j]][k] / data[[var_i]][k]
        }
      }
      
      counter <- counter + 1
      
      # Store results in list
      res_list[[counter]] <- diff
      names(res_list)[counter] <- paste0(var_i, "X", var_j)
    }
  }
  return(res_list)
}


## On absolute difference scale
res_diff <- FL_SmithsS_generateContrasts(data = top8, contrast = "absolute_diff")
names(res_diff)


# Convert list of posterior samples to data frame
res_diff_full <- as.data.frame(res_diff)
res_diff_full

# Make long format and separate comparison levels
res_diff_full <- res_diff_full %>%
  pivot_longer(cols = everything(), names_to = "contrast", values_to = "S") %>%
  separate(contrast, c("base_level", "contrast_level"), sep = "X", remove = FALSE) %>%
  mutate(base_level = recode(base_level, "hardworking" = "hard working", 
                             "good.conscience" = "good conscience")) %>%
  mutate(contrast_level = recode(contrast_level, "hardworking" = "hard working", 
                                 "good.conscience" = "good conscience")) %>%
  mutate(base_level = factor(base_level, levels = c("hard working", "kind", "helpful", "modest",
                                                    "respectful", "honest", "intelligent", "good conscience"))) %>%
  mutate(contrast_level = factor(contrast_level, levels = c("hard working", "kind", "helpful", "modest",
                                                            "respectful", "honest", "intelligent", "good conscience"))) %>%
  arrange(base_level, contrast_level)
res_diff_full

# Keep just 'hard-working' as the base level
res_diff_work <- res_diff_full[res_diff_full$base_level == "hard working", ]
head(res_diff_work)

# Summary of comparisons
res_diff_work %>%
  group_by(contrast_level) %>%
  summarise(median = median(S), lci_2.5 = quantile(S, 0.025), uci_97.5 = quantile(S, 0.975),
            lci_12.5 = quantile(S, 0.125), uci_87.5 = quantile(S, 0.875),
            lci_25 = quantile(S, 0.25), uci_75 = quantile(S, 0.75))

# Plot of results
(p_dist_diff <- ggplot(res_diff_work, aes(x = S, y = forcats::fct_rev(contrast_level), 
                                          fill = forcats::fct_rev(contrast_level))) +
    stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_qi") +
    scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.9) +
    guides(fill = "none") +
    labs(x = expression("Difference in" ~ italic(hat(S)) ~ "cultural salience"), y = "Item contrast") +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.line = element_line(colour = "black"),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()))

# Save as PDF
pdf(file = "./FL Uncertainty/top8_contrasts_hardworking.pdf", height = 8, width = 8)
p_dist_diff
dev.off()



## On ratio difference scale
res_diff_ratio <- FL_SmithsS_generateContrasts(data = top8, contrast = "ratio_diff")
names(res_diff_ratio)


# Convert list of posterior samples to data frame
res_diff_full <- as.data.frame(res_diff_ratio)
res_diff_full

# Make long format and separate comparison levels
res_diff_full <- res_diff_full %>%
  pivot_longer(cols = everything(), names_to = "contrast", values_to = "S") %>%
  separate(contrast, c("base_level", "contrast_level"), sep = "X", remove = FALSE) %>%
  mutate(base_level = recode(base_level, "hardworking" = "hard working", 
                             "good.conscience" = "good conscience")) %>%
  mutate(contrast_level = recode(contrast_level, "hardworking" = "hard working", 
                                 "good.conscience" = "good conscience")) %>%
  mutate(base_level = factor(base_level, levels = c("hard working", "kind", "helpful", "modest",
                                                    "respectful", "honest", "intelligent", "good conscience"))) %>%
  mutate(contrast_level = factor(contrast_level, levels = c("hard working", "kind", "helpful", "modest",
                                                            "respectful", "honest", "intelligent", "good conscience"))) %>%
  arrange(base_level, contrast_level)
res_diff_full

# Keep just 'hard-working' as the base level
res_diff_work <- res_diff_full[res_diff_full$base_level == "hard working", ]
head(res_diff_work)

# Summary of comparisons
res_diff_work %>%
  group_by(contrast_level) %>%
  summarise(median = median(S), lci_2.5 = quantile(S, 0.025), uci_97.5 = quantile(S, 0.975))

# Plot of results
(p_dist_diff_ratio <- ggplot(res_diff_work, aes(x = S, y = forcats::fct_rev(contrast_level), 
                                          fill = forcats::fct_rev(contrast_level))) +
    stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_qi") +
    scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.9) +
    guides(fill = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), limits = c(0, 2.5)) +
    labs(x = expression("Ratio of difference in" ~ italic(hat(S)) ~ "cultural salience"), y = "Item contrast") +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          axis.line = element_line(colour = "black"),
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()))

# Save as PDF
pdf(file = "./FL Uncertainty/top8_contrasts_hardworking_ratio.pdf", height = 8, width = 8)
p_dist_diff_ratio
dev.off()



################################################################################
### Comparison across groups - Here, whether nominating the trait of 'kindness' differs by sex in Tyvans

## Prepare the data

# Keep just 'kindness' item salience
FL.sal0.kind <- FL.sal0[, c("Subject", "kind")]
head(FL.sal0.kind)

summary(FL.sal0.kind)
hist(FL.sal0.kind$kind)

# Read in demographic data and merge by ID
demo <- read.delim("demo.txt")
head(demo)
str(demo)

demo <- demo[, c("Subj", "Sex")]
colnames(demo)[1] <- "Subject"
head(demo)

kind.sex <- left_join(FL.sal0.kind, demo, by = "Subject")
head(kind.sex)
summary(kind.sex)

# Remove if sex is missing
kind.sex <- kind.sex[!is.na(kind.sex$Sex), ]
summary(kind.sex)
table(kind.sex$Sex)

# Descriptive stats and histograms of kindness by sex
by(kind.sex$kind, kind.sex$Sex, summary)

par(mar = c(5, 5, 1, 2))

hist(kind.sex$kind[kind.sex$Sex == 0], freq = FALSE, 
     ylim = c(0, 20), col = rgb(0,0,1,1/4), breaks = 20,
     main = "", xlab = "Item salience",
     cex.axis = 1.2, cex.lab = 1.5)
hist(kind.sex$kind[kind.sex$Sex == 1], freq = FALSE, col = rgb(1,0,0,1/4), breaks = 20, add = TRUE)
legend("topright", legend = c("Female", "Male"), 
       col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)),
       pch = 15, pt.cex = 2)

# Save this plot
pdf(file = "./FL Uncertainty/kindness_hist_bySex.pdf", height = 4, width = 6)

par(mar = c(5, 5, 1, 2))

hist(kind.sex$kind[kind.sex$Sex == 0], freq = FALSE, 
     ylim = c(0, 20), col = rgb(0,0,1,1/4), breaks = 20,
     main = "", xlab = "Item salience",
     cex.axis = 1.2, cex.lab = 1.5)
hist(kind.sex$kind[kind.sex$Sex == 1], freq = FALSE, col = rgb(1,0,0,1/4), breaks = 20, add = TRUE)
legend("topright", legend = c("Female", "Male"), 
       col = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)),
       pch = 15, pt.cex = 2)

dev.off()



#### First as a predictor in a regression context (no bootstrapping)

### ZOIB model

## ZOIB model with sex as predictor
zoib_sex <- brm(formula = bf(
  kind ~ Sex,
  phi ~ Sex, 
  zoi ~ Sex, 
  coi ~ Sex), 
  data = kind.sex,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 5748,
  family = zero_one_inflated_beta())

get_prior(zoib_sex)

summary(zoib_sex)

## Posterior predictions from this model by sex

# For males and females together
post_zoib_sex <- predict(zoib_sex, newdata = as.data.frame(cbind(Sex = kind.sex$Sex)), summary = FALSE)
head(post_zoib_sex)

# Calculate Smith's S for each posterior sample, separately by sex
S_post_zoib_male <- rep(NA, nrow(post_zoib_sex))
S_post_zoib_female <- rep(NA, nrow(post_zoib_sex))
S_post_zoib_sexDiff <- rep(NA, nrow(post_zoib_sex))
S_post_zoib_sexDiff_ratio <- rep(NA, nrow(post_zoib_sex))
for (i in 1:nrow(post_zoib_sex)) {
  S_post_zoib_male[i] <- mean(post_zoib_sex[i, ][kind.sex$Sex == 0])
  S_post_zoib_female[i] <- mean(post_zoib_sex[i, ][kind.sex$Sex == 1])
  S_post_zoib_sexDiff[i] <- S_post_zoib_female[i] - S_post_zoib_male[i]
  S_post_zoib_sexDiff_ratio[i] <- S_post_zoib_female[i] / S_post_zoib_male[i]
}

hist(S_post_zoib_male)
hist(S_post_zoib_female)
hist(S_post_zoib_sexDiff)
hist(S_post_zoib_sexDiff_ratio)

summary(S_post_zoib_male)
summary(S_post_zoib_female)
summary(S_post_zoib_sexDiff)
summary(S_post_zoib_sexDiff_ratio)

quantile(S_post_zoib_male, c(0.025, 0.5, 0.975))
quantile(S_post_zoib_female, c(0.025, 0.5, 0.975))
quantile(S_post_zoib_sexDiff, c(0.025, 0.5, 0.975))
quantile(S_post_zoib_sexDiff_ratio, c(0.025, 0.5, 0.975))


### Ordered Beta model - Note that as males and females differ in probabilities of 0s and 1s, have to include sex as a predictor for these cut-points (and manually specify priors), else results are biased.
sex_mod <- ordbetareg(
  formula = bf(kind ~ Sex, 
               phi ~ Sex,
               cutzero ~ Sex,
               cutone ~ Sex), 
  data = kind.sex,
  true_bounds = c(0, 1),
  manual_prior = c(prior(normal(0, 1.5), class = "b", dpar = "cutzero"),
                   prior(normal(0, 1.5), class = "b", dpar = "cutone")),
  phi_reg = "both",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 2576)

get_prior(sex_mod)

summary(sex_mod)

## Posterior predictions from this model by sex

# For males and females together
post <- predict(sex_mod, newdata = as.data.frame(cbind(Sex = kind.sex$Sex)), summary = FALSE)
head(post)

# Calculate Smith's S for each posterior sample, separately by sex
S_post_m <- rep(NA, nrow(post))
S_post_f <- rep(NA, nrow(post))
S_post_diff <- rep(NA, nrow(post))
S_post_ratio <- rep(NA, nrow(post))
for (i in 1:nrow(post)) {
  S_post_m[i] <- mean(post[i, ][kind.sex$Sex == 0])
  S_post_f[i] <- mean(post[i, ][kind.sex$Sex == 1])
  S_post_diff[i] <- S_post_f[i] - S_post_m[i]
  S_post_ratio[i] <- S_post_f[i] / S_post_m[i]
}

hist(S_post_m)
hist(S_post_f)
hist(S_post_diff)
hist(S_post_ratio)

summary(S_post_m)
summary(S_post_f)
summary(S_post_diff)
summary(S_post_ratio)

quantile(S_post_m, c(0.025, 0.5, 0.975))
quantile(S_post_f, c(0.025, 0.5, 0.975))
quantile(S_post_diff, c(0.025, 0.5, 0.975))
quantile(S_post_ratio, c(0.025, 0.5, 0.975))


#### Also in both groups separately (with bootstrapping)

### Bootstrapping

## Model for males
kind.male <- kind.sex[kind.sex$Sex == 0, ]

# Perform bootstrapping
set.seed(6543)
S_boot.m <- rep(NA, 1000)
for (i in 1:1000) {
  boot <- sample(kind.male$kind, size = nrow(kind.male), replace = TRUE)
  S_boot.m[i] <- mean(boot)
}

# 95% percentile intervals
quantile(S_boot.m, c(0.025, 0.5, 0.975))


## Model for females
kind.female <- kind.sex[kind.sex$Sex == 1, ]

# Perform bootstrapping
set.seed(3456)
S_boot.f <- rep(NA, 1000)
for (i in 1:1000) {
  boot <- sample(kind.female$kind, size = nrow(kind.female), replace = TRUE)
  S_boot.f[i] <- mean(boot)
}

# 95% percentile intervals
quantile(S_boot.f, c(0.025, 0.5, 0.975))


## Difference between male and female Smith's S values
S_boot.diff <- rep(NA, length(S_boot.m))
S_boot.diff_ratio <- rep(NA, length(S_boot.m))
for (i in 1:length(S_boot.m)) {
  S_boot.diff[i] <- S_boot.f[i] - S_boot.m[i]
  S_boot.diff_ratio[i] <- S_boot.f[i] / S_boot.m[i]
}

hist(S_boot.diff)
summary(S_boot.diff)
quantile(S_boot.diff, c(0.025, 0.5, 0.975))

hist(S_boot.diff_ratio)
summary(S_boot.diff_ratio)
quantile(S_boot.diff_ratio, c(0.025, 0.5, 0.975))



### ZOIB model

## Model for males
kind.male <- kind.sex[kind.sex$Sex == 0, ]

# ZOIB model
zoib_male <- brm(formula = bf(
  kind ~ 1,
  phi ~ 1, 
  zoi ~ 1, 
  coi ~ 1), 
  data = kind.male,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 2357,
  family = zero_one_inflated_beta())

summary(zoib_male)

# Posterior predictions from this model
post_zoib.m <- predict(zoib_male, summary = FALSE)
head(post_zoib.m)

# Calculate Smith's S for each posterior sample
S_post_zoib.m <- rep(NA, nrow(post_zoib.m))
for (i in 1:nrow(post_zoib.m)) {
  S_post_zoib.m[i] <- mean(post_zoib.m[i, ])
}

hist(S_post_zoib.m)
summary(S_post_zoib.m)
quantile(S_post_zoib.m, c(0.025, 0.5, 0.975))


## Model for females
kind.female <- kind.sex[kind.sex$Sex == 1, ]

# ZOIB model
zoib_female <- brm(formula = bf(
  kind ~ 1,
  phi ~ 1, 
  zoi ~ 1, 
  coi ~ 1), 
  data = kind.female,
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 863,
  family = zero_one_inflated_beta())

summary(zoib_female)

# Posterior predictions from this model
post_zoib.f <- predict(zoib_female, summary = FALSE)
head(post_zoib.f)

# Calculate Smith's S for each posterior sample
S_post_zoib.f <- rep(NA, nrow(post_zoib.f))
for (i in 1:nrow(post_zoib.f)) {
  S_post_zoib.f[i] <- mean(post_zoib.f[i, ])
}

hist(S_post_zoib.f)
summary(S_post_zoib.f)
quantile(S_post_zoib.f, c(0.025, 0.5, 0.975))


## Difference between male and female Smith's S values
S_post_zoib.diff <- rep(NA, length(S_post_zoib.m))
S_post_zoib.diff_ratio <- rep(NA, length(S_post_zoib.m))
for (i in 1:length(S_post_zoib.m)) {
  S_post_zoib.diff[i] <- S_post_zoib.f[i] - S_post_zoib.m[i]
  S_post_zoib.diff_ratio[i] <- S_post_zoib.f[i] / S_post_zoib.m[i]
}

hist(S_post_zoib.diff)
summary(S_post_zoib.diff)
quantile(S_post_zoib.diff, c(0.025, 0.5, 0.975))

hist(S_post_zoib.diff_ratio)
summary(S_post_zoib.diff_ratio)
quantile(S_post_zoib.diff_ratio, c(0.025, 0.5, 0.975))


### Ordered Beta model

## Model for males
kind.male <- kind.sex[kind.sex$Sex == 0, ]

# Ordered Beta model
ordBeta_male <- ordbetareg(
  formula = bf(kind ~ 1,
               phi ~ 1),
  data = kind.male,
  true_bounds = c(0, 1),
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 8568)

summary(ordBeta_male)

# Posterior predictions from this model
post_ordBeta.m <- predict(ordBeta_male, summary = FALSE)
head(post_ordBeta.m)

# Calculate Smith's S for each posterior sample
S_post_ordBeta.m <- rep(NA, nrow(post_ordBeta.m))
for (i in 1:nrow(post_ordBeta.m)) {
  S_post_ordBeta.m[i] <- mean(post_ordBeta.m[i, ])
}

hist(S_post_ordBeta.m)
summary(S_post_ordBeta.m)
quantile(S_post_ordBeta.m, c(0.025, 0.5, 0.975))


## Model for females
kind.female <- kind.sex[kind.sex$Sex == 1, ]

# Ordered Beta model
ordBeta_female <- ordbetareg(
  formula = bf(kind ~ 1,
               phi ~ 1),
  data = kind.female,
  true_bounds = c(0, 1),
  manual_prior = c(prior(student_t(3, 0, 2.5), class = Intercept, dpar = phi)),
  phi_reg = "only",
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 6445)

summary(ordBeta_female)

# Posterior predictions from this model
post_ordBeta.f <- predict(ordBeta_female, summary = FALSE)
head(post_ordBeta.f)

# Calculate Smith's S for each posterior sample
S_post_ordBeta.f <- rep(NA, nrow(post_ordBeta.f))
for (i in 1:nrow(post_ordBeta.f)) {
  S_post_ordBeta.f[i] <- mean(post_ordBeta.f[i, ])
}

hist(S_post_ordBeta.f)
summary(S_post_ordBeta.f)
quantile(S_post_ordBeta.f, c(0.025, 0.5, 0.975))


## Difference between male and female Smith's S values
S_post_ordBeta.diff <- rep(NA, length(S_post_ordBeta.m))
S_post_ordBeta.diff_ratio <- rep(NA, length(S_post_ordBeta.m))
for (i in 1:length(S_post_ordBeta.m)) {
  S_post_ordBeta.diff[i] <- S_post_ordBeta.f[i] - S_post_ordBeta.m[i]
  S_post_ordBeta.diff_ratio[i] <- S_post_ordBeta.f[i] / S_post_ordBeta.m[i]
}

hist(S_post_ordBeta.diff)
summary(S_post_ordBeta.diff)
quantile(S_post_ordBeta.diff, c(0.025, 0.5, 0.975))

hist(S_post_ordBeta.diff_ratio)
summary(S_post_ordBeta.diff_ratio)
quantile(S_post_ordBeta.diff_ratio, c(0.025, 0.5, 0.975))



###########################################################################
#### Example of applying methods above to free-list metrics other than Smith's S - First focusing on Jaccard's similarity and conceptual overlap (note that the data, code and example have been adapted here from Purzycki's 'Ethnographic Free-list Data' book, chapter 4 - For the original code, see https://github.com/bgpurzycki/free-list_QASS)

## Install and load the 'eulerr' package
#install.packages("eulerr")
library(eulerr)

## Read in the data
tyva <- read.delim("tyva_domains.txt", sep = "\t")

# This dataset contains 5 variables: Individual ID ("CERCID"), the order in which items were listed ("Order"), up to 5 things the big/moralising god dislikes ("BGD"; Buddha, in this context), up to 5 things the local god dislikes ("LGD"; spirit-masters, in this context), and up ot 5 things the police dislike ("POD"). These items have been coded into a general rubric (e.g., 'Morality', 'Virtue', 'Ecology', etc.) 
head(tyva)
str(tyva)


### For the example here, we are interested in the conceptual overlap of concerned regarding 'Morality' between big/moralising gods, local gods, and the police.

## First, have to do quite a bit of data processing to know, for each person, whether they selected 'Morality' or not

# Presence matrix, saying whether each person listed each domain or not
BGDbin <- FreeListTable(tyva, CODE = "BGD", Order = "Order", Subj = "CERCID", tableType = "PRESENCE")
LGDbin <- FreeListTable(tyva, CODE = "LGD", Order = "Order", Subj = "CERCID", tableType = "PRESENCE")
PODbin <- FreeListTable(tyva, CODE = "POD", Order = "Order", Subj = "CERCID", tableType = "PRESENCE")

# Rename 'Morality' variables to keep track of which target they refer to
BGDbin$BGDmor <- BGDbin$Morality # presence to merge later
LGDbin$LGDmor <- LGDbin$Morality # presence to merge later
PODbin$PODmor <- PODbin$Morality # presence to merge later

# Merge into one dataset
binmerge1 <- merge(BGDbin, LGDbin, by = "Subject", all = TRUE)
binmerge2 <- merge(binmerge1, PODbin, by = "Subject", all = TRUE)
lab0 <- c("Subject", "BGDmor", "LGDmor", "PODmor")
tyva_morality <- binmerge2[lab0]

# This is now a presence matrix, detailing whether each individual listed 'Morality' as a concern for each target.
head(tyva_morality)
summary(tyva_morality)


## Calculate Jaccard's similarity across the three domains

# Function to calculate Jaccard's similarity
jaccard <- function(var1, var2) {
  sums <- rowSums(cbind(var1, var2), na.rm = T)
  similarity <- length(sums[sums==2])
  total <- length(sums[sums==1]) + similarity
  similarity/total
}

# Calculate Jaccard's similarity manually
jaccard(tyva_morality$BGDmor, tyva_morality$LGDmor)
jaccard(tyva_morality$BGDmor, tyva_morality$PODmor)
jaccard(tyva_morality$LGDmor, tyva_morality$PODmor)

# Plot the conceptual overlap
fit <- euler(c(Buddha = 1, Spirits = 1, Police = 1, 
               "Buddha&Police" = .69, "Spirits&Police" = .13, "Buddha&Spirits" = .15))
fit$original.values # need to remove #'s we don't want in plot
valuestoplot <- c(NA, NA, NA, .15, .69, .13, NA)
plot(fit, fills = c("white", "darkgray", "lightgray"), quantities = valuestoplot)



####### Incorporating uncertainty into Jaccard's similarity

### Start with simplest method, bootstrapping

# Randomly sample 1,000 times from the observed item presence matrix (with replacement), calculating Jaccard's similarity for each - Also calculate differences in Jaccard's similarity as well
set.seed(45654)

# Placeholders for Jaccard's similarity
Jac_boot_BGvsLG <- rep(NA, 1000)
Jac_boot_BGvsPOL <- rep(NA, 1000)
Jac_boot_LGvsPOL <- rep(NA, 1000)

# Bootstrapping
for (i in 1:1000) {
  
  # Sample dataframe
  boot <- slice_sample(tyva_morality, n = nrow(tyva_morality), replace = TRUE)
  
  # Calcualte Jaccard's similarity
  Jac_boot_BGvsLG[i] <- jaccard(boot$BGDmor, boot$LGDmor)
  Jac_boot_BGvsPOL[i] <- jaccard(boot$BGDmor, boot$PODmor)
  Jac_boot_LGvsPOL[i] <- jaccard(boot$LGDmor, boot$PODmor)
}

# Results for overall Jaccard's similarity (95% percentile intervals)
quantile(Jac_boot_BGvsLG, c(0.025, 0.5, 0.975))
quantile(Jac_boot_BGvsPOL, c(0.025, 0.5, 0.975))
quantile(Jac_boot_LGvsPOL, c(0.025, 0.5, 0.975))


# Plot (with uncertainty values)
fit <- euler(c(Buddha = 1, Spirits = 1, Police = 1, 
               "Buddha&Police" = .69, "Spirits&Police" = .13, "Buddha&Spirits" = .14))
fit$original.values # need to remove #'s we don't want in plot
valuestoplot <- c(NA, NA, NA, "0.14\n[0.06, 0.25]", "0.69\n[0.58, 0.80]", "0.13\n[0.05, 0.21]", NA)
plot(fit, fills = c("white", "darkgray", "lightgray"), quantities = valuestoplot)



### More principled method: Multi-level Bayesian logistic regression

# Convert data to long format
tyva_morality_long <- tyva_morality %>%
  pivot_longer(cols = !Subject, names_to = "domain", values_to = "presence") %>%
  mutate(domain = factor(domain, levels = c("BGDmor", "LGDmor", "PODmor"))) %>%
  arrange(Subject)
tyva_morality_long


# Logistic model, with domain as random-effect within participants to allow this relationship to vary by participant.

# check the default priors (ideally these should be more informative, especially for a model with moderate complexity such as this one)
get_prior(bf(presence ~ domain + (domain | Subject)),
          data = tyva_morality_long,
          family = bernoulli())

# Running the model with more informative priors
logit_mod <- brm(
  formula = bf(presence ~ domain + (domain | Subject)),
  data = tyva_morality_long,
  prior = c(prior(normal(0, 1.5), class = b), # Weakly-regularising normal prior for the fixed-effects
            prior(normal(0, 1.5), class = Intercept), # Weakly-regularising normal prior for intercept
            prior(exponential(1), class = sd), # Exponential prior for variation in random effects
            prior(lkj(2), class = cor)), # LKJ prior for correlation between random effects
  chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 12345,
  family = bernoulli()
)

summary(logit_mod)


## One divergent transition noted, but no major warnings/issues and R-hat values seem okay (max = 1.01), but effective sample sizes for some of the random effects are a bit low (e.g., ~500), so, to be safe, will increase the number of sampling chains (from 1000 to 2000) and increase sampling depth (from default of 0.8 to 0.9)
logit_mod2 <- brm(
  formula = bf(presence ~ domain + (domain | Subject)),
  data = tyva_morality_long,
  prior = c(prior(normal(0, 1.5), class = b), # Weakly-regularising normal prior for the fixed-effects
            prior(normal(0, 1.5), class = Intercept), # Weakly-regularising normal prior for intercept
            prior(exponential(1), class = sd), # Exponential prior for variation in random effects
            prior(lkj(2), class = cor)), # LKJ prior for correlation between random effects
  chains = 4, iter = 3000, warmup = 1000, cores = 4, seed = 123456,
  control = list(adapt_delta = 0.9),
  family = bernoulli()
)

# This is looking better now, as no divergent transitions, all r-hat values are 1.00, and effective sample sizes are larger.
summary(logit_mod2)


## Posterior predictions from this model
df_pred <- tyva_morality_long[, c("Subject", "domain")]
post_logit <- predict(logit_mod2, newdata = df_pred, summary = FALSE)
head(post_logit)


## Calculate Jaccard's similarity for each posterior sample

# Placeholders for Jaccard's similarity
Jac_logit_BGvsLG <- rep(NA, nrow(post_logit))
Jac_logit_BGvsPOL <- rep(NA, nrow(post_logit))
Jac_logit_LGvsPOL <- rep(NA, nrow(post_logit))

# Looping over each posterior sample
for (i in 1:nrow(post_logit)) {
  
  # Append posterior predictions to dataframe, then convert to wide format
  temp <- df_pred
  temp$presence <- post_logit[i, ]
  
  temp_wide <- pivot_wider(temp, id_cols = "Subject", names_from = "domain", values_from = "presence")
  
  # Calculate Jaccard's similarity
  Jac_logit_BGvsLG[i] <- jaccard(temp_wide$BGDmor, temp_wide$LGDmor)
  Jac_logit_BGvsPOL[i] <- jaccard(temp_wide$BGDmor, temp_wide$PODmor)
  Jac_logit_LGvsPOL[i] <- jaccard(temp_wide$LGDmor, temp_wide$PODmor)
}

# Results for overall Jaccard's similarity (95% credible intervals)
quantile(Jac_logit_BGvsLG, c(0.025, 0.5, 0.975))
quantile(Jac_logit_BGvsPOL, c(0.025, 0.5, 0.975))
quantile(Jac_logit_LGvsPOL, c(0.025, 0.5, 0.975))


# Plot (with uncertainty values)
fit <- euler(c(Buddha = 1, Spirits = 1, Police = 1, 
               "Buddha&Police" = .60, "Spirits&Police" = .16, "Buddha&Spirits" = .18))
fit$original.values # need to remove #'s we don't want in plot
valuestoplot <- c(NA, NA, NA, "0.18\n[0.07, 0.31]", "0.60\n[0.46, 0.73]", "0.16\n[0.07, 0.26]", NA)
plot(fit, fills = c("white", "darkgray", "lightgray"), quantities = valuestoplot)


### Comparison between results

## Big/moralising gods vs local gods
jaccard(tyva_morality$BGDmor, tyva_morality$LGDmor) # Jaccard (without uncertainty)
quantile(Jac_boot_BGvsLG, c(0.025, 0.5, 0.975)) # Bootstrapping
quantile(Jac_logit_BGvsLG, c(0.025, 0.5, 0.975)) # Logistic MLM

## Big/moralising gods vs police
jaccard(tyva_morality$BGDmor, tyva_morality$PODmor) # Jaccard (without uncertainty)
quantile(Jac_boot_BGvsPOL, c(0.025, 0.5, 0.975)) # Bootstrapping
quantile(Jac_logit_BGvsPOL, c(0.025, 0.5, 0.975)) # Logistic MLM

## Local gods vs police
jaccard(tyva_morality$LGDmor, tyva_morality$PODmor) # Jaccard (without uncertainty)
quantile(Jac_boot_LGvsPOL, c(0.025, 0.5, 0.975)) # Bootstrapping
quantile(Jac_logit_LGvsPOL, c(0.025, 0.5, 0.975)) # Logistic MLM

# In this example, the logistic MLM performs okay, but bootstrapping does appear more accurate/less biased



#################
#### Further example of applying methods above to free-list metrics other than Smith's S - Now focusing on Cultural FST and the partitioning of variance between vs within societies (as above, note that the data, code and example have been adapted here from Purzycki's 'Ethnographic Free-list Data' book, chapter 4 - For the original code, see https://github.com/bgpurzycki/free-list_QASS)

## Read in the data
dat_fst <- read.csv("Cross-cultural_ERM1.csv", sep = ";")

# This dataset is of data from the Evolution of Religion and Morality (ERM) project. Here, we're only interested in free-list on what moralising gods dislike, so will make a presence/absence matrix of all categories from this question.
dat_fst_pres <- FreeListTable(dat_fst, CODE = "BGD", Order = "Order", 
                              Subj = "CERCID", tableType = "PRESENCE", 
                              GROUPING = "Culture")
head(dat_fst_pres)

# Exclude data if participant did not respond to this question
dat_fst_pres$freq <- rowSums(dat_fst_pres[,3:12])
dat_fst_pres <- dat_fst_pres[dat_fst_pres$freq != 0, ]  

# Keep just morality item
dat_fst_pres <- dat_fst_pres[, c("Subject", "Group", "Morality")]
head(dat_fst_pres)


## Function to calculate FST between two societies
FST <- function(ni, nj, xi, xj){ 
  pi <- xi/ni
  pj <- xj/nj
  pbar <- (xi + xj)/(ni + nj)
  num <- (ni/(ni + nj))*(pi - pbar)^2 + (nj/(ni + nj))*(pj - pbar)^2
  denom <- pbar*(1 - pbar)
  FST <- num/denom
  return(FST)
}

## Function to loop over multiple societies and create matrix of cultural FST estimates (based on summary table of data)
fstmatrix <- function(tslab, matrixtype = NULL) {
  m <- matrix(NA, nrow = nrow(tslab), ncol = nrow(tslab))
  rownames(m) <- colnames(m) <- rownames(tslab)
  for(i in 1:nrow(m)) {
    for(j in 1:ncol(m)) {
      m[i, j] <- FST(tslab[i, 1], tslab[j, 1], tslab[i, 2], tslab[j, 2])   
    }
  }
  if(!is.null(matrixtype)) {
    if(matrixtype == "upper") {
      m[lower.tri(m)] <- 0
    }
    if(matrixtype == "lower") {
      m[upper.tri(m)] <- 0
    }
  }
  return(m)
}

# Table summarising morality data for these societies
mortab <- table(dat_fst_pres$Morality, dat_fst_pres$Group)
mortab1 <- mortab[2,]
size <- colSums(mortab)
slab <- (as.data.frame(rbind(size, mortab1)))
tslab <- as.data.frame(t(slab))
tslab

# Matrix of cultural FST estimates, based on this summary table
fstmatrix(tslab)


#### Methods to propagate uncertainty 

### First, bootstrapping

## Function to generate bootstrap samples across a range of groups/societies
FST_boot_sampling <- function(data, group_var, target_var, iterations = 1000, seed) {
  
  # Set seed and list to store results in
  set.seed(seed)
  res_list <- list()
  
  for (i in 1:length(unique(data[[group_var]]))) {
    
    # Take each group in turn
    var <- unique(data[[group_var]])[i]
    dat_temp <- data[data[[group_var]] == var, ]
    
    print(paste0("On group ", i, ": ", var))
    
    # Perform bootstrapping on each site and store total number of times target mentioned
    boot_temp <- rep(NA, iterations)
    
    for (j in 1:iterations) {
      boot <- sample(dat_temp[[target_var]], size = nrow(dat_temp), replace = TRUE)
      boot_temp[j] <- sum(boot)
    }
    
    res_list[[i]] <- boot_temp
    names(res_list)[i] <- var
  }
  return(res_list)
}

# Run this function to generate estimates of the numbers of individuals selecting 'Morality' in each society
boot_samples <- FST_boot_sampling(data = dat_fst_pres, group_var = "Group", target_var = "Morality", iteration = 1000, seed = 123)
str(boot_samples)


## Function to calculate cultural FST between all groups, across these bootstrapped samples
fst_uncert <- function(dat_list, orig_data, group_var) {
  
  # List to store results in
  res_list <- list()
  
  # Keep a counter of number of FST calculations/lists to create
  counter <- 0
  
  # Loop over all combinations of groups (skipping if same group)
  for (i in 1:length(names(dat_list))) {
    
    var_i <- names(dat_list)[i]
    
    for (j in 1:length(names(dat_list))) {
      
      var_j <- names(dat_list)[j]
      
      # Skip if same pairing
      if (var_i == var_j) next
      
      print(paste0("On variables ", i, ": ", var_i, " and ", j, ": ", var_j))
      
      # Vector to store results in
      fst_temp <- rep(NA, length(dat_list[[var_i]]))
      
      # Calculate FST in each sample
      for (k in 1:length(dat_list[[var_i]])) {
        fst_temp[k] <- FST(ni = nrow(orig_data[orig_data[[group_var]] == var_i, ]), 
                           nj = nrow(orig_data[orig_data[[group_var]] == var_j, ]), 
                           xi = dat_list[[var_i]][k], 
                           xj = dat_list[[var_j]][k])
      }
      
      counter <- counter + 1
      
      # Store results in list
      res_list[[counter]] <- fst_temp
      names(res_list)[counter] <- paste0(var_i, "X", var_j)
    }
  }
  return(res_list)
}

# Run this function to calculate FST estimates
fst_boot <- fst_uncert(dat_list = boot_samples, orig_data = dat_fst_pres, group_var = "Group")
str(fst_boot)

# Summarise cultural FST for all societies
(fst_res <- round(t(as.data.frame(lapply(fst_boot, quantile, probs = c(0.025, 0.5, 0.975)))), 2))

# Put median Cultural FST estimates in a matrix (adapting Ben's earlier 'fstmatrix' function)
fstmatrix_uncert <- function(dat_res, dat_list, est = "50%", lower_quant = "2.5%", upper_quant = "97.5%", matrixtype = NULL) {
  
  # Extract variable names and re-convert '.'s back to spaces
  names <- names(dat_list)
  names <- gsub(" ", ".", names)
  
  # Populate matrix with cultural FST estimates, and uncertainty intervals
  m <- matrix(NA, nrow = length(names(dat_list)), ncol = length(names(dat_list)))
  rownames(m) <- colnames(m) <- names(dat_list)
  for(i in 1:nrow(m)) {
    for(j in 1:ncol(m)) {
      if (i == j) {
        m[i, j] <- NA
      } else if (i != j) {
        target <- paste0(names[i], "X", names[j])
        m[i, j] <- paste0(as.data.frame(dat_res)[[est]][rownames(dat_res) == target], " [",
                          as.data.frame(dat_res)[[lower_quant]][rownames(dat_res) == target], "-",
                          as.data.frame(dat_res)[[upper_quant]][rownames(dat_res) == target], "]")
      }
    }
  }
  if(!is.null(matrixtype)) {
    if(matrixtype == "upper") {
      m[lower.tri(m)] <- 0
    }
    if(matrixtype == "lower") {
      m[upper.tri(m)] <- 0
    }
  }
  return(m)
}

# FST matrices for different uncertainty intervals (2.5%, median/50% and 97.5%)
(fstmat <- fstmatrix_uncert(dat_res = fst_res, dat_list = boot_samples, 
                            est = "50%", lower_quant = "2.5%", upper_quant = "97.5%"))


### Alternative method using Bayesian logistic regression model (using 'brms'), followed by sampling posterios predictions for each group

# Note the '0 +' notation to exclude the traditional intercept and include use 'index' notation to estimate intercept separately for all levels of 'Group' (not just relative to a reference)
logit_moral <- brm(formula = bf(Morality ~ 0 + Group),
                   data = dat_fst_pres,
                   chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 12321,
                   family = bernoulli())

summary(logit_moral)

# Posterior predictions for each group
post_moral <- predict(logit_moral, newdata = as.data.frame(cbind(Group = dat_fst_pres$Group)), summary = FALSE)
head(post_moral)

# Hive these off into separate vectors in a list (so same structure as the bootstrapping output)
post_byGroup <- function(post_data, orig_data, group_var) {
  
  # Create a list to store results in
  dat_list <- list()
  
  # Loop over each group/society and extract posterior predictions, summarising total number within each group
  for (i in 1:length(unique(orig_data[[group_var]]))) {
    temp <- rep(NA, nrow(post_data))
    var <- unique(orig_data[[group_var]])[i]
    for (j in 1:nrow(post_data)) {
      temp[j] <- sum(post_data[j, ][orig_data[[group_var]] == var])
    }
    dat_list[[i]] <- temp
    names(dat_list)[i] <- var
  }
  return(dat_list)
}

# Run this function to split by group
post_list <- post_byGroup(post_data = post_moral, orig_data = dat_fst_pres, group_var = "Group")
str(post_list)

# Use 'fst_uncert' function above to calculate cultural FST between each group
fst_logit <- fst_uncert(dat_list = post_list, orig_data = dat_fst_pres, group_var = "Group")
str(fst_logit)

# Summarise cultural FST
(fst_res_logit <- round(t(as.data.frame(lapply(fst_logit, quantile, probs = c(0.025, 0.5, 0.975)))), 2))

# Summarise cultural FST in a matrix with uncertainty intervals
(fstmat_logit <- fstmatrix_uncert(dat_res = fst_res_logit, dat_list = post_list, 
                                  est = "50%", lower_quant = "2.5%", upper_quant = "97.5%"))

# Essentially the same results as with bootstrapping, but with wider uncertainty intervals
fstmat


### Alternative method using group/society as a random effect
logit_moral_re <- brm(formula = bf(Morality ~ 1 + (1 | Group)),
                      data = dat_fst_pres,
                      chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 12321,
                      family = bernoulli())

summary(logit_moral_re)

# Posterior predictions for each group
post_moral_re <- predict(logit_moral_re, newdata = as.data.frame(cbind(Group = dat_fst_pres$Group)), summary = FALSE)
head(post_moral_re)

# Run function to split posterior predictions by group
post_list_re <- post_byGroup(post_data = post_moral_re, orig_data = dat_fst_pres, group_var = "Group")
str(post_list_re)

# Use 'fst_uncert' function above to calculate cultural FST between each group
fst_logit_re <- fst_uncert(dat_list = post_list_re, orig_data = dat_fst_pres, group_var = "Group")
str(fst_logit_re)

# Summarise cultural FST
(fst_res_logit_re <- round(t(as.data.frame(lapply(fst_logit_re, quantile, probs = c(0.025, 0.5, 0.975)))), 2))

# Summarise cultural FST in a matrix with uncertainty intervals
(fstmat_logit_re <- fstmatrix_uncert(dat_res = fst_res_logit_re, dat_list = post_list_re, 
                                     est = "50%", lower_quant = "2.5%", upper_quant = "97.5%"))

# Compare with bootstrapping and fixed effects regression model
fstmat_logit
fstmat

# Again, broadly similar results to both (albeit with wider uncertainty intervals than bootstrapping again) - However, compared to bootstrapping and groups as fixed effects, in the random effects model some of the cultural FST estimates are closer to the null (e.g., for Coastal vs Inland Tanna cultural FST = 0.33 for bootstrapping and fixed effect model, but 0.28 for random effects model); Is because multi-level models are pooling/regularising, thus weakening the impact of more extreme values. Can see this when looking at the predictions from the random effects model, as the estimates are less extreme for the samples with more and fewer people listing 'Morality' (e.g., in Coastal Tanna 41/42 participants listed 'Morality'; in the bootstrap and fixed effects samples the median value was also 41, while in the random effects posterior sample it was 40 [with lower mean and minimum estimates as well]).
tslab
lapply(boot_samples, summary)
lapply(post_list, summary)
lapply(post_list_re, summary)
