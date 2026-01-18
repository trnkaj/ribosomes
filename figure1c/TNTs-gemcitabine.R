library(rethinking)
library(tidyverse)

library(readxl)

TNTs_gem <- read_excel("~/Downloads/ribosomes-main/figure1c/TNTs-gemcitabine.xlsx", 
                       skip = 1)

#filter data
data <- TNTs_gem %>% select(`gemcitabine (uM)`, `total count`) %>% rename(gem = `gemcitabine (uM)`, tot = `total count`)
#standardise data
stddata <- list(
  gemstd = standardize(data$gem),
  totstd = standardize(data$tot)
)
#stddata

#estimator
bTNT <- ulam(
  alist(
    totstd ~ dnorm(mu, sd),
    mu <- a + b*gemstd,
    a ~ dnorm(0,1),
    b ~ dnorm(0,1),
    sd ~ dexp(1)
  ), data=stddata, chains = 4, cores = 4, iter=4000
  
)

precis(bTNT)

#posterior distributions
post <- extract.samples(bTNT)
#post
#graph standardised coefficients
dens(post$a)
dens(post$b)
#unstandardise coefficient b
tot_sd <- sd(TNTs_gem$`total count`)
tot_mean <- mean(TNTs_gem$`total count`)
b_ustd <- post$b*tot_sd/gem_sd
b_ustd_mean <-  mean(b_ustd)
b_ustd_cri95 <- HPDI(b_ustd, prob=0.95)

print(paste("Estimated increase in TNTs per μM Gemcitabine:", round(b_ustd_mean, 2), 
            " (95% CrI:", round(b_ustd_cri95[1], 2), ",", round(b_ustd_cri95[2], 2), ")"))

#graph unstandardised coefficient b
dens(b_ustd, show.zero = TRUE, main="Posterior distribution of regression slope",
     ylab="", xlab="Change in TNTs per μM Gemcitabine", cex.axis=1, adj=1.5, xlim=c(-20,300))
shade(density(b_ustd, adj=1.5), HPDI(b_ustd, prob=0.89))

#check priors
prior <- extract.prior(bTNT)
mu <- link(bTNT, post=prior, data=list(gemstd=c(-2,2)))
plot(NULL, xlim=c(-2,2),ylim=c(-2,2))
for (i in 1:50) lines (c(-2,2), mu[i,], col=col.alpha("black",0.4))

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
#gem_destand
#link model
mu <- link(bTNT, data=list(gemstd=gem))
tot_stand_mu <- apply(mu, 2, mean)
tot_stand_pi <- apply(mu, 2, PI)
#destandardise output values
tot_mu <- tot_stand_mu*tot_sd+tot_mean
tot_pi <- tot_stand_pi*tot_sd+tot_mean
#plot prediction lines
lines(gem_destand, tot_mu, lty=2, lwd=1.5)
shade(tot_pi, gem_destand, xpd=TRUE)
mtext(paste("Estimated increase in TNTs per μM Gemcitabine:", round(b_ustd_mean, 2), 
      " (95% CrI:", round(b_ustd_cri95[1], 2), ",", round(b_ustd_cri95[2], 2), ")"), side=3, line=0.5, cex=1.2)