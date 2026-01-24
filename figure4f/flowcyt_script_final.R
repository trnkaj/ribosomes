#############################################
# RUN ONLY ONCE
install.packages("BiocManager")
BiocManager::install("flowCore")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggcyto")
#############################################

library(flowCore)
library(ggcyto)
library(tidyverse)

# 1) Set gates from controls

# 1.1) Set gates from control cells with only DAPI staining

# load the data
file_dapi_control <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp1_blue_time_24_BL24_.fcs", 
                                transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_dapi_control <- exprs(file_dapi_control)

# create data frame from the extracted numbers
df_dapi_control <- as.data.frame(exprs_dapi_control)


### SINGLETS GATING
# define maximal FSC-W
singlets_gate <- 140000

# set filter1 based on maximal FSC-W and sort cells into two groups with smaller and bigger FSC-W
df_dapi_control <- df_dapi_control %>% 
  mutate(filter1 = if_else(`FSC-W` > singlets_gate, "out", "in"))

# plot scatter plot with the singlets gate
df_dapi_control %>% 
  ggplot(aes(x = `FSC-A`,
             y = `FSC-W`,
             colour = filter1))+
  geom_point(size = 0.2)+
  theme_classic()+
  geom_hline(aes(yintercept = singlets_gate), colour = "red")+
  scale_colour_manual(values = c("black", "red"))+
  theme(legend.position = "none")

# plot PE of cells that were left after singlets and debris gating
df_dapi_control %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `FSC-H`,
             y = `PE-A`))+
  geom_point(size = 0.2)+
  theme_classic()

# create a new column with added minimal value of PE-A + 1 and then log-transform it
df_dapi_control <- df_dapi_control %>% 
  mutate(`log PE-A` = log(df_dapi_control$`PE-A` + abs(min(df_dapi_control$`PE-A`)) + 1))

# filter only cells after singlets and debris gating and within 5 sd from mean PE-A
gate_pe <- df_dapi_control %>% 
  dplyr::filter(filter1 == "in") %>% ####### TADY MĚNIT FILTER1/FILTER2
  dplyr::filter(`PE-A` < (exp(mean(df_dapi_control$`log PE-A`) + 10*sd(df_dapi_control$`log PE-A`)) - abs(min(df_dapi_control$`PE-A`)) - 1))

# set the gate as maximal value of PE-A
gate_pe <- max(gate_pe$`PE-A`)



# 1.2) Set gates from control cells with only PE staining

# load the data
file_pe_control <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp1_TMR_time_24_TMR24_.fcs", 
                              transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_pe_control <- exprs(file_pe_control)

# create data frame from the extracted numbers
df_pe_control <- as.data.frame(exprs_pe_control)


### SINGLETS GATING
# define maximal FSC-W
singlets_gate <- 140000

# set filter1 based on maximal FSC-W and sort cells into two groups with smaller and bigger FSC-W
df_pe_control <- df_pe_control %>% 
  mutate(filter1 = if_else(`FSC-W` > singlets_gate, "out", "in"))

# plot scatter plot with the singlets gate
df_pe_control %>% 
  ggplot(aes(x = `FSC-A`,
             y = `FSC-W`,
             colour = filter1))+
  geom_point(size = 0.2)+
  theme_classic()+
  geom_hline(aes(yintercept = singlets_gate), colour = "red")+
  scale_colour_manual(values = c("black", "red"))+
  theme(legend.position = "none")

# plot DAPI of cells that were left after singlets and debris gating
df_pe_control %>% 
  dplyr::filter(filter1 == "in") %>% ####### TADY MĚNIT FILTER1/FILTER2
  ggplot(aes(x = `FSC-H`,
             y = `DAPI-A`))+
  geom_point(size = 0.2)+
  theme_classic()

# create a new column with added minimal value of DAPI-A + 1 and then log-transform it
df_pe_control <- df_pe_control %>% 
  mutate(`log DAPI-A` = log(df_pe_control$`DAPI-A` + abs(min(df_pe_control$`DAPI-A`)) + 1))

# filter only cells after singlets and debris gating and within 5 sd from mean FITC-A
gate_dapi <- df_pe_control %>% 
  dplyr::filter(filter1 == "in") %>% ####### TADY MĚNIT FILTER1/FILTER2
  dplyr::filter(`DAPI-A` < (exp(mean(df_pe_control$`log DAPI-A`) + 10*sd(df_pe_control$`log DAPI-A`)) - abs(min(df_pe_control$`DAPI-A`)) - 1))

# set the gate as maximal value of FITC-A
gate_dapi <- max(gate_dapi$`DAPI-A`)





# 2) Set gates for experimental samples

# 2.1) DSHaloRPS9Blue-20241023-exp1_HaloRPS9_blue_NT_time_24_E24_.fcs

# load the data
file_E1 <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp1_HaloRPS9_blue_NT_time_24_E24_.fcs", 
                    transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_E1 <- exprs(file_E1)

# create data frame from the extracted numbers
df_E1 <- as.data.frame(exprs_E1)

# set sequence for singlets gating
singlets_maxima_grid <- seq(90000, 205000, 5000)
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create a new data frame with repeated rows for each y-intercept of singlets_maxima
df_E1_expanded <- df_E1 %>%
  crossing(singlets_maxima = singlets_maxima_grid)

# create the plot with facets
df_E1_expanded %>% 
  ggplot(aes(x = `FSC-A`, y = `FSC-W`)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_hline(aes(yintercept = singlets_maxima), colour = "red") +
  facet_wrap(~ singlets_maxima, ncol = 4)  # adjust ncol for the desired layout

# create an empty data frame to store the results
analysis_E1 <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_E1 <- df_E1 %>% 
    mutate(filter1 = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_E1$filter1 == "in") 
  num_double <- sum(df_E1$filter1 == "in" & df_E1$`DAPI-A` > gate_dapi & df_E1$`PE-A` > gate_pe) 
  # Add the results to the data frame
  analysis_E1 <- rbind(analysis_E1, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

# plot the results of sensitivity analysis
chosen_gate <- 140000
analysis_E1 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(x = "singlets gate",
       y = "double positive cells / all gated cells")+
  scale_y_continuous(limits = c(0, 0.01))+
  scale_x_continuous(breaks = seq(90000, 205000, 10000))+
  geom_vline(xintercept = chosen_gate)


analysis_E1 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp1_HaloRPS9_blue_NT_time_24_E24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 80))

df_E1 <- df_E1 %>% 
  mutate(filter1 = if_else(`FSC-W` > chosen_gate, "out", "in"))

# plot singlets gating with the chosen gate
df_E1 %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp1_HaloRPS9_blue_NT_time_24_E24_")

# plot the resulting graph (with all gating)
df_E1 %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.3)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  theme(legend.position = "none")







# 2.2) DSHaloRPS9Blue-20241023-exp2_HaloRPS9_blue_NT_time_24_E24_.fcs

# load the data
file_E2 <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp2_HaloRPS9_blue_NT_time_24_E24_.fcs", 
                    transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_E2 <- exprs(file_E2)

# create data frame from the extracted numbers
df_E2 <- as.data.frame(exprs_E2)

# set sequence for singlets gating
singlets_maxima_grid <- seq(90000, 205000, 5000)
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create a new data frame with repeated rows for each y-intercept of singlets_maxima
df_E2_expanded <- df_E2 %>%
  crossing(singlets_maxima = singlets_maxima_grid)

# create the plot with facets
df_E2_expanded %>% 
  ggplot(aes(x = `FSC-A`, y = `FSC-W`)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_hline(aes(yintercept = singlets_maxima), colour = "red") +
  facet_wrap(~ singlets_maxima, ncol = 4)  # adjust ncol for the desired layout

# create an empty data frame to store the results
analysis_E2 <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_E2 <- df_E2 %>% 
    mutate(filter1 = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_E2$filter1 == "in") 
  num_double <- sum(df_E2$filter1 == "in" & df_E2$`DAPI-A` > gate_dapi & df_E2$`PE-A` > gate_pe) 
  # Add the results to the data frame
  analysis_E2 <- rbind(analysis_E2, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

# plot the results of sensitivity analysis
analysis_E2 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(x = "singlets gate",
       y = "double positive cells / all gated cells")+
  scale_y_continuous(limits = c(0, 0.003))+
  scale_x_continuous(breaks = seq(90000, 205000, 10000))+
  geom_vline(xintercept = 140000)

analysis_E2 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp2_HaloRPS9_blue_NT_time_24_E24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 60))

chosen_gate <- 140000

df_E2 <- df_E2 %>% 
  mutate(filter1 = if_else(`FSC-W` > chosen_gate, "out", "in"))

# plot singlets gating with the chosen gate
df_E2 %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp2_HaloRPS9_blue_NT_time_24_E24_")

# plot the resulting graph (with all gating)
df_E2 %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.3)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  theme(legend.position = "none")






# 2.3) DSHaloRPS9Blue-20241023-exp3_HaloRPS9_blue_NT_time_24_E24_.fcs

# load the data
file_E3 <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp3_HaloRPS9_blue_NT_time_24_E24_.fcs", 
                    transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_E3 <- exprs(file_E3)

# create data frame from the extracted numbers
df_E3 <- as.data.frame(exprs_E3)

# set sequence for singlets gating
singlets_maxima_grid <- seq(90000, 205000, 5000)
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create a new data frame with repeated rows for each y-intercept of singlets_maxima
df_E3_expanded <- df_E3 %>%
  crossing(singlets_maxima = singlets_maxima_grid)

# create the plot with facets
df_E3_expanded %>% 
  ggplot(aes(x = `FSC-A`, y = `FSC-W`)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_hline(aes(yintercept = singlets_maxima), colour = "red") +
  facet_wrap(~ singlets_maxima, ncol = 4)  # adjust ncol for the desired layout

# create an empty data frame to store the results
analysis_E3 <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_E3 <- df_E3 %>% 
    mutate(filter1 = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_E3$filter1 == "in") 
  num_double <- sum(df_E3$filter1 == "in" & df_E3$`DAPI-A` > gate_dapi & df_E3$`PE-A` > gate_pe)
  # Add the results to the data frame
  analysis_E3 <- rbind(analysis_E3, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

# plot the results of sensitivity analysis
analysis_E3 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(x = "singlets gate",
       y = "double positive cells / all gated cells")+
  scale_y_continuous(limits = c(0, 0.01))+
  scale_x_continuous(breaks = seq(90000, 205000, 10000))+
  geom_vline(xintercept = chosen_gate)

analysis_E3 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp3_HaloRPS9_blue_NT_time_24_E24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 65))

chosen_gate <- 140000

df_E3 <- df_E3 %>% 
  mutate(filter1 = if_else(`FSC-W` > chosen_gate, "out", "in"))

# plot singlets gating with the chosen gate
df_E3 %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp3_HaloRPS9_blue_NT_time_24_E24_")

# plot the resulting graph (with all gating)
df_E3 %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.3)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  theme(legend.position = "none")









# 2.4) DSHaloRPS9Blue-20241023-exp1_HaloRPS9_blue_GEM_time_24_F24_.fcs

# load the data
file_F1 <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp1_HaloRPS9_blue_GEM_time_24_F24_.fcs", 
                    transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_F1 <- exprs(file_F1)

# create data frame from the extracted numbers
df_F1 <- as.data.frame(exprs_F1)

# set sequence for singlets gating
singlets_maxima_grid <- seq(90000, 205000, 5000)
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create a new data frame with repeated rows for each y-intercept of singlets_maxima
df_F1_expanded <- df_F1 %>%
  crossing(singlets_maxima = singlets_maxima_grid)

# create the plot with facets
df_F1_expanded %>% 
  ggplot(aes(x = `FSC-A`, y = `FSC-W`)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_hline(aes(yintercept = singlets_maxima), colour = "red") +
  facet_wrap(~ singlets_maxima, ncol = 4)  # adjust ncol for the desired layout

# create an empty data frame to store the results
analysis_F1 <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_F1 <- df_F1 %>% 
    mutate(filter1 = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_F1$filter1 == "in") 
  num_double <- sum(df_F1$filter1 == "in" & df_F1$`DAPI-A` > gate_dapi & df_F1$`PE-A` > gate_pe) 
  # Add the results to the data frame
  analysis_F1 <- rbind(analysis_F1, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

# plot the results of sensitivity analysis
analysis_F1 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(x = "singlets gate",
       y = "double positive cells / all gated cells")+
  scale_y_continuous(limits = c(0, 0.006))+
  scale_x_continuous(breaks = seq(90000, 205000, 10000))+
  geom_vline(xintercept = chosen_gate)

analysis_F1 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp1_HaloRPS9_blue_GEM_time_24_F24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 90))

chosen_gate <- 140000

df_F1 <- df_F1 %>% 
  mutate(filter1 = if_else(`FSC-W` > chosen_gate, "out", "in"))

# plot singlets gating with the chosen gate
df_F1 %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp1_HaloRPS9_blue_GEM_time_24_F24_")

# plot the resulting graph (with all gating)
df_F1 %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.3)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  theme(legend.position = "none")






# 2.5) DSHaloRPS9Blue-20241023-exp2_HaloRPS9_blue_GEM_time_24_F24_.fcs

# load the data
file_F2 <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp2_HaloRPS9_blue_GEM_time_24_F24_.fcs", 
                    transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_F2 <- exprs(file_F2)

# create data frame from the extracted numbers
df_F2 <- as.data.frame(exprs_F2)

# set sequence for singlets gating
singlets_maxima_grid <- seq(90000, 205000, 5000)
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create a new data frame with repeated rows for each y-intercept of singlets_maxima
df_F2_expanded <- df_F2 %>%
  crossing(singlets_maxima = singlets_maxima_grid)

# create the plot with facets
df_F2_expanded %>% 
  ggplot(aes(x = `FSC-A`, y = `FSC-W`)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_hline(aes(yintercept = singlets_maxima), colour = "red") +
  facet_wrap(~ singlets_maxima, ncol = 4)  # adjust ncol for the desired layout

# create an empty data frame to store the results
analysis_F2 <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_F2 <- df_F2 %>% 
    mutate(filter1 = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_F2$filter1 == "in") 
  num_double <- sum(df_F2$filter1 == "in" & df_F2$`DAPI-A` > gate_dapi & df_F2$`PE-A` > gate_pe)
  # Add the results to the data frame
  analysis_F2 <- rbind(analysis_F2, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

# plot the results of sensitivity analysis
analysis_F2 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(x = "singlets gate",
       y = "double positive cells / all gated cells")+
  scale_y_continuous(limits = c(0, 0.01))+
  scale_x_continuous(breaks = seq(90000, 205000, 10000))+
  geom_vline(xintercept = chosen_gate)

analysis_F2 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp2_HaloRPS9_blue_GEM_time_24_F24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 70))

chosen_gate <- 140000

df_F2 <- df_F2 %>% 
  mutate(filter1 = if_else(`FSC-W` > chosen_gate, "out", "in"))

# plot singlets gating with the chosen gate
df_F2 %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp2_HaloRPS9_blue_GEM_time_24_F24_")

# plot the resulting graph (with all gating)
df_F2 %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.3)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  theme(legend.position = "none")







# 2.6) DSHaloRPS9Blue-20241023-exp3_HaloRPS9_blue_GEM_time_24_F24_.fcs

# load the data
file_F3 <- read.FCS("~/Downloads/ribosomes-main/figure3f/DSHaloRPS9Blue-20241023-exp3_HaloRPS9_blue_GEM_time_24_F24_.fcs", 
                    transformation = FALSE, truncate_max_range = FALSE)
# extract numbers from fcs file
exprs_F3 <- exprs(file_F3)

# create data frame from the extracted numbers
df_F3 <- as.data.frame(exprs_F3)

# set sequence for singlets gating
singlets_maxima_grid <- seq(90000, 205000, 5000)
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create a new data frame with repeated rows for each y-intercept of singlets_maxima
df_F3_expanded <- df_F3 %>%
  crossing(singlets_maxima = singlets_maxima_grid)

# create the plot with facets
df_F3_expanded %>% 
  ggplot(aes(x = `FSC-A`, y = `FSC-W`)) +
  geom_point(size = 0.2) +
  theme_classic() +
  geom_hline(aes(yintercept = singlets_maxima), colour = "red") +
  facet_wrap(~ singlets_maxima, ncol = 4)  # adjust ncol for the desired layout

# create an empty data frame to store the results
analysis_F3 <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_F3 <- df_F3 %>% 
    mutate(filter1 = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_F3$filter1 == "in") 
  num_double <- sum(df_F3$filter1 == "in" & df_F3$`DAPI-A` > gate_dapi & df_F3$`PE-A` > gate_pe) 
  # Add the results to the data frame
  analysis_F3 <- rbind(analysis_F3, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

chosen_gate <- 140000
# plot the results of sensitivity analysis
analysis_F3 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(x = "singlets gate",
       y = "double positive cells / all gated cells")+
  scale_y_continuous(limits = c(0, 0.008))+
  scale_x_continuous(breaks = seq(90000, 205000, 10000))+
  geom_vline(xintercept = chosen_gate)

analysis_F3 %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp3_HaloRPS9_blue_GEM_time_24_F24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 120))


df_F3 <- df_F3 %>% 
  mutate(filter1 = if_else(`FSC-W` > chosen_gate, "out", "in"))

# plot singlets gating with the chosen gate
df_F3 %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp3_HaloRPS9_blue_GEM_time_24_F24_")

# plot the resulting graph (with all gating)
df_F3 %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.3)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  theme(legend.position = "none")












# 3) Back check gates in controls

# 3.1) DAPI control

# set sequence for singlets gating
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create an empty data frame to store the results
analysis_dapi <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_dapi_control <- df_dapi_control %>% 
    mutate(filter_check = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_dapi_control$filter_check == "in") 
  num_double <- sum(df_dapi_control$filter_check == "in" & df_dapi_control$`DAPI-A` > gate_dapi & df_dapi_control$`PE-A` > gate_pe) 
  # Add the results to the data frame
  analysis_dapi <- rbind(analysis_dapi, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

chosen_gate <- 140000

# plot the results of sensitivity analysis
analysis_dapi %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "PE exp1_blue_time_24_BL24_",
       x = "singlets gate",
       y = "double positive cells/all gated cells")+
  scale_y_continuous(limits = c(0, 0.011))+
  geom_vline(aes(xintercept = chosen_gate), colour = "red")


analysis_dapi %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "PE exp1_blue_time_24_BL24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 10))+
  geom_vline(aes(xintercept = chosen_gate), colour = "red")


# plot singlets gating with the chosen gate
df_dapi_control %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp1_blue_time_24_BL24_")


# plot the resulting graph (with all gating)
df_dapi_control %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.2)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  labs(title = "exp1_blue_time_24_BL24_")








# 3.2) PE control

# set sequence for singlets gating
singlets_maxima_analysis <- seq(90000, 205000, 500)

# create an empty data frame to store the results
analysis_pe <- data.frame(singlets_gate = numeric(), cells_gated = numeric(), cells_double = numeric())

# use a for loop to calculate "in" and "out" for each singlets maximum value
for (gate in singlets_maxima_analysis) {
  # create the filter1 column based on the current singlets maximum value
  df_pe_control <- df_pe_control %>% 
    mutate(filter_check = if_else(`FSC-W` > gate, "out", "in"))
  num_gated <- sum(df_pe_control$filter_check == "in") 
  num_double <- sum(df_pe_control$filter_check == "in" & df_pe_control$`DAPI-A` > gate_dapi & df_pe_control$`PE-A` > gate_pe) 
  # Add the results to the data frame
  analysis_pe <- rbind(analysis_pe, data.frame(singlets_gate = gate, cells_gated = num_gated, cells_double = num_double))
}

chosen_gate <- 140000

# plot the results of sensitivity analysis
analysis_pe %>% 
  ggplot(aes(x = singlets_gate, y = cells_double/cells_gated))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp1_TMR_time_24_TMR24_",
       x = "singlets gate",
       y = "double positive cells/all gated cells")+
  scale_y_continuous(limits = c(0, 0.01))+
  geom_vline(aes(xintercept = chosen_gate), colour = "red")

analysis_pe %>% 
  ggplot(aes(x = singlets_gate, y = cells_double))+
  geom_point(size = 0.2)+
  theme_bw()+
  scale_x_continuous(breaks = singlets_maxima_grid)+
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5))+
  labs(title = "exp1_TMR_time_24_TMR24_",
       x = "singlets gate",
       y = "double positive cells")+
  scale_y_continuous(limits = c(0, 10))+
  geom_vline(aes(xintercept = chosen_gate), colour = "red")


# plot singlets gating with the chosen gate
df_pe_control %>% 
  ggplot(aes(x = `FSC-A`, 
             y = `FSC-W`, 
             colour = filter1)) +
  geom_point(size = 0.2) +
  theme_classic() +
  scale_colour_manual(values = c("black", "red"))+
  geom_hline(aes(yintercept = chosen_gate), colour = "red")+
  geom_hline(aes(yintercept = 90000), colour = "red", linetype = "dashed")+
  geom_hline(aes(yintercept = 205000), colour = "red", linetype = "dashed")+
  theme(legend.position = "none")+
  labs(title = "exp1_TMR_time_24_TMR24_")

# plot the resulting graph (with all gating)
df_pe_control %>% 
  dplyr::filter(filter1 == "in") %>% 
  ggplot(aes(x = `DAPI-A`,
             y = `PE-A`))+
  geom_point(size = 0.2)+
  theme_classic()+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_vline(aes(xintercept = gate_dapi))+
  geom_hline(aes(yintercept = gate_pe))+
  labs(title = "exp1_TMR_time_24_TMR24_")




# 4) Number of double positive cells

df_list <- list(df_E1, df_E2, df_E3, df_F1, df_F2, df_F3)

# Applying the condition to each data frame
result <- sapply(df_list, function(df) {
  as.numeric(sum(df$filter1 == "in" & df$`DAPI-A` > gate_dapi & df$`PE-A` > gate_pe))
})

doublestained <- data.frame(num_of_cells = result, 
                               exp = c("E1", "E2", "E3", "F1", "F2", "F3"), 
                               group = rep(c("E", "F"), each = 3),
                               gemcitabine = rep(c("NT", "gemcitabine"), each = 3))


doublestained %>% 
  ggplot(aes(x = group, 
             y = num_of_cells,
             shape = fct_relevel(gemcitabine, c("NT", "gemcitabine"))))+
  geom_point(size = 3)+
  theme_classic()+
  scale_y_continuous(limits = c(0, 80))+
  scale_shape_manual(values = c(21, 16))+
  scale_x_discrete(labels = c("E" = "control", "F" = "0.5 μM gemcitabine"))+
  labs(x = "",
       y = "number of double positive cells")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10))+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 11,
    shape = 95,
    fill = "black")

#STATISTICAL ANALYSIS OF MEANS
#create data frame of numbers of doublestained cells for bayesian model (modified from "result" on line 911)
doublestained_b <- data.frame(num_of_cells = result, 
                            exp = c("E1", "E2", "E3", "F1", "F2", "F3"), 
                            group = rep(c(1, 2), each = 3),
                            gemcitabine = rep(c(0, 0.5), each = 3))

#check the variance of data to select a suitable prior for sigma
print(sd(doublestained_b$num_of_cells[doublestained_b$group == 1]))
print(sd(doublestained_b$num_of_cells[doublestained_b$group == 2]))

#bayesian model for means comparison
comp_means <- ulam(
  alist(
    num_of_cells ~ dnorm(mu, sigma),
    mu <- a[group],
    a[group] ~ dnorm(50, 100),
    sigma ~ dhalfnorm(0,10)
  ),
  data = doublestained_b,
  chains = 4,
  cores = 4,
  iter = 4000
)

precis(comp_means, depth=2)

#extract posterior samples
post_means <- extract.samples(comp_means)
#plot the means with credible intervals
dens(post_means$a[,1], 
     xlab = "number of double positive cells", 
     main = "Posterior distributions of the means", adj=1.5, xlim=c(0,80))
dens(post_means$a[,2], add = TRUE, adj=1.5)


#calculate posterior difference between means
diff_means <- post_means$a[,2] - post_means$a[,1]
#calculate 95% credible interval
ci_diff <- HPDI(diff_means, prob = 0.95)
print(paste("The estimated difference in means is", round(mean(diff_means),2),"(95% CrI:", round(ci_diff[1],2),",",round(ci_diff[2],2),")"))
#plot posterior distribution of the difference with shaded credible interval
dens(diff_means, 
     xlab = "difference in means (number of double positive cells)", 
     main = "Posterior distribution of the difference in means", adj=1.5)
shade(density(diff_means, adj=1.5),ci_diff)
abline(v = 0, lty = 2)

# CHECK MODEL PREDICTIVE PERFORMANCE WITH DIFFERENT PRIORS FOR SIGMA
# --- 1. Define the models with different sigma priors ---

# Model A: Tight Prior
# Assumes noise is small (~1)
m_tight <- ulam(
  alist(
    num_of_cells ~ dnorm(mu, sigma),
    mu <- a[group],
    a[group] ~ dnorm(50, 100),
    sigma ~ dhalfnorm(0,1) 
  ),
  data = doublestained_b, chains = 4, cores = 4, iter = 4000,
  log_lik = TRUE
)

# Model B: Medium Prior
# Assumes noise is moderate (~7-10).
m_medium <- ulam(
  alist(
    num_of_cells ~ dnorm(mu, sigma),
    mu <- a[group],
    a[group] ~ dnorm(50, 100),
    sigma ~ dhalfnorm(0, 10)
  ),
  data = doublestained_b, chains = 4, cores = 4, iter = 4000,
  log_lik = TRUE
)

# Model C: Loose Prior
# Assumes noise could be huge (~20+)
m_loose <- ulam(
  alist(
    num_of_cells ~ dnorm(mu, sigma),
    mu <- a[group],
    a[group] ~ dnorm(50, 100),
    sigma ~ dhalfnorm(0, 20)
  ),
  data = doublestained_b, chains = 4, cores = 4, iter = 4000,
  log_lik = TRUE
)

# --- 2. Compare them ---
# You can use func=WAIC (default) or func=LOO (better for outliers)
comparison <- compare(m_tight, m_medium, m_loose, func = WAIC)

print(comparison)
