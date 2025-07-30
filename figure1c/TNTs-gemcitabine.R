library(rethinking)
library(tidyverse)

library(readxl)
TNTs_gem <- read_excel("./figure1c/TNTs-gem.xlsx", 
                       skip = 1)

#filter data
data <- TNTs_gem %>% select(`gemcitabine (uM)`, `total count`) %>% rename(gem = `gemcitabine (uM)`, tot = `total count`)
#standardise data
stddata <- list(
  gemstd = standardize(data$gem),
  totstd = standardize(data$tot)
)
stddata

#estimator
bTNT <- ulam(
  alist(
    totstd ~ dnorm(mu, sd),
    mu <- a + b*gemstd,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sd ~ dexp(1)
  ), data=stddata, chains = 4, cores = 4
  
)

precis(bTNT)

#posterior distributions
post <- extract.samples(bTNT)
post
#graph standardised coefficients
dens(post$a)
dens(post$b)
#unstandardise coefficient b
b_ustd <- post$b*tot_sd/gem_sd
#graph unstandardised coefficient b
dens(b_ustd, show.zero = TRUE, #main="Posterior distibution of slope"
     ylab="", xlab="", cex.axis=1.5)

#check priors
prior <- extract.prior(bTNT)
mu <- link(bTNT, post=prior, data=list(gemstd=c(-2,2)))
plot(NULL, xlim=c(-2,2),ylim=c(-2,2))
for (i in 1:50) lines (c(-2,2), mu[i,], col=col.alpha("black",0.4))

#plot(stddata$totstd~stddata$gemstd)

#posterior prediction
plot(data$tot~data$gem,xlab="", ylab="", cex.lab = 1.5, 
     pch = 19, col = rgb(0, 0, 0, 0.5), #main = "Dose-dependent increase in TNTs (89 % credible interval)", cex.main = 2, 
     cex.axis = 1.5
     )
#set x-axis values(gemcitabine)
ns <- 100
gem <- seq(from=-1.15, to = 1.5, length.out=ns)
#destandardise gem
gem_mean <- mean(TNTs_gem$`gemcitabine (uM)`)
gem_sd <- sd(TNTs_gem$`gemcitabine (uM)`)
gem_destand <- gem*gem_sd+gem_mean
gem_destand
#link model
mu <- link(bTNT, data=list(gemstd=gem))
tot_stand_mu <- apply(mu, 2, mean)
tot_stand_pi <- apply(mu, 2, PI)
#destandardise output values
tot_sd <- sd(TNTs_gem$`total count`)
tot_mean <- mean(TNTs_gem$`total count`)
tot_mu <- tot_stand_mu*tot_sd+tot_mean
tot_pi <- tot_stand_pi*tot_sd+tot_mean
#plot prediction lines
lines(gem_destand, tot_mu, lty=2, lwd=1.5)
shade(tot_pi, gem_destand, xpd=TRUE)
