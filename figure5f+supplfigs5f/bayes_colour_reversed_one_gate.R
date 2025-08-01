library(rethinking)
library(tidyverse)


####load experimental data, csv produced by scripts from figure5de+supplfigs5de
a <- read.csv("./figure5de+supplfigs5de/colour-reversed-exp23_27_28_one_gate.csv")



#modify experimental dataframe a (from green_sirna sciript)
a$treat <- ifelse(a$dataset == "sirna" | a$dataset == "mix_sirna", 1, 0)
a$mix_control <- ifelse(a$dataset == "mix_control", 1, 0)
a$mix_sirna <- ifelse(a$dataset == "mix_sirna", 1, 0)
#create coluomn exp_no for indexing "experiment" as indeces 1 to n
a$exp_no <- as.numeric(as.factor(a$experiment)) #convert experiment to numeric factor for indexing
a$treatment <- ifelse(a$dataset == "control", 1, ifelse(a$dataset == "sirna", 2, ifelse(a$dataset == "mix_control", 3, 4))) #create treatment column with 1 for control, 2 for sirna, 3 for mix_control and 4 for mix_sirna
#standardize log10_alexa values
a$log10_alexa_std <- standardize(a$log10_alexa) 

###check the data
print(mean(subset(a, subset = treat == 0 & mix_control == 0 & exp_no == 1, select = log10_alexa)$log10_alexa)) #mean log10_alexa for control cells in experiment 1
print(mean(subset(a, subset = treat == 0 & mix_control == 0 & exp_no == 2, select = log10_alexa)$log10_alexa)) #mean log10_alexa for control cells in experiment 2
print(mean(subset(a, subset = treat == 0 & mix_control == 0 & exp_no == 3, select = log10_alexa)$log10_alexa)) #mean log10_alexa for control cells in experiment 3

#select relevant columns for the model
data_a <- a %>%
  select(log10_alexa_std, exp_no, treat, mix_control, mix_sirna, treatment) # %>%
#  mutate(across(c(treat, mix_control, mix_sirna), as.inte)) #convert binary indicators to numeric


#####varying intercepts
m1_reverse_colour_one_gate <- ulam(
  alist(
    log10_alexa_std ~ dnorm(mu, sigma),
    mu <- a[exp_no] + b_treat*treat + b_mix_control*mix_control + b_mix_sirna*mix_sirna, #b_treat is the effect of sirna, treat = 1 for sirna, 0 for control; b_mix_control is the effect of mixing on control, mix_coontrol = 1 for mixed control cells, 0 for mixed sirna cells or non_mixed cells; b_mix_sirna is the effect of mixing on sirna, mix_sirna = 1 for mixed sirna cells, 0 for mixed control cells or non-mixed cells
    a[exp_no] ~ dnorm(a_bar, sigma_a),
    b_treat ~ dnorm(0,1),
    b_mix_control ~ dnorm(0,1),
    b_mix_sirna ~ dnorm(0,1),
    a_bar ~ dnorm(0, 5),
    sigma_a ~ dexp(1),
    sigma ~  dexp(1)
  ), 
  data = data_a, 
  chains = 4, 
  cores = 4
  # control = list(adapt_delta = 0.99)
  #start = list(
  #  a = rep(3,3),
  #  a_bar = 3, 
  #  b_treat = 0,
  #  b_mix_control = 0,
  #  b_mix_sirna = 0,
  #  sigma_a = 0.1,
  #  sigma = 0.1
  #)
)

plot(m1_reverse_colour_one_gate, depth = 2)
precis(m1_reverse_colour_one_gate, depth = 2)
#trankplot(ribo_m6_green_sirna, depth = 2)

post_m1_rcog <- extract.samples(m1_reverse_colour_one_gate)

#destandardise samples
post_m1_rcog$a_bar_norm <- post_m1_rcog$a_bar * sd(a$log10_alexa) + mean(a$log10_alexa)
#post_m1_rcog$a <- post_m1_rcog$a * sd(a$log10_alexa) + mean(a$log10_alexa)
post_m1_rcog$b_treat_norm <- post_m1_rcog$b_treat * sd(a$log10_alexa)
post_m1_rcog$b_mix_control_norm <- post_m1_rcog$b_mix_control * sd(a$log10_alexa)
post_m1_rcog$b_mix_sirna_norm <- post_m1_rcog$b_mix_sirna * sd(a$log10_alexa)

#predict values for measurements with treat = 0, mix_control = 0, mix_sirna = 0, i.e. control
mu_control <- 10^(post_m1_rcog$a_bar_norm + post_m1_rcog$b_treat_norm*0 + post_m1_rcog$b_mix_control_norm*0 + post_m1_rcog$b_mix_sirna_norm*0)
mu_sirna <- 10^(post_m1_rcog$a_bar_norm + post_m1_rcog$b_treat_norm*1 + post_m1_rcog$b_mix_control_norm*0 + post_m1_rcog$b_mix_sirna_norm*0)
mu_control_mix <- 10^(post_m1_rcog$a_bar_norm + post_m1_rcog$b_treat_norm*0 + post_m1_rcog$b_mix_control_norm*1 + post_m1_rcog$b_mix_sirna_norm*0)
mu_sirna_mix <-10^(post_m1_rcog$a_bar_norm + post_m1_rcog$b_treat_norm*1 + post_m1_rcog$b_mix_control_norm*0 + post_m1_rcog$b_mix_sirna_norm*1)

#plot density of the posterior distributions for control, sirna, control_mix and sirna_mix
#save to pdf
pdf("./figure5f+supplfigs5f/bayes_colour_reversed_one_gate_posterior.pdf", width = 8, height = 6)
dens(mu_control, col = "#60AF1A", lwd = 2, #xlab = "Alexa647 Intensity", ylab = "Density", main = "Posterior Distributions of Alexa647 Intensity, m1_reverse_colour_one_gate", 
     xlim = c(1, 2100), ylim = c(0, 0.023))
dens(mu_sirna, col = "#0B7D9E", lwd = 2, add = TRUE)
dens(mu_control_mix, col = "#A8DC79", lwd = 2, add = TRUE)
dens(mu_sirna_mix, col = "#00C0BA", lwd = 2, add = TRUE)
dev.off()

#plot subtracted values for mu_control - mu_control_mix and mu_sirna - mu_sirna_mix
#save to pdf
pdf("./figure5f+supplfigs5f/bayes_colour_reversed_one_gate_posterior_difference.pdf", width = 8, height = 6)
dens(mu_control_mix - mu_control, col = "#A8DC79", lwd = 2, #xlab = "Difference in Alexa647 Intensity", ylab = "Density", main = "Posterior Distributions of Differences in Alexa647 Intensity, m1_reverse_colour_one_gate", 
     xlim = c(-800, 300), ylim = c(0, 0.03))
dens(mu_sirna_mix - mu_sirna, col = "#00C0BA", lwd = 2, add = TRUE)
abline(v = 0, col = "black", lwd = 2, lty = 2)
dev.off()