library(flowCore)  # opening FSC files
library(ggplot2)
library(scales)    # in plots labels = comma
library(dplyr)

# Define function to load and extract raw data from FCS files
# Convert to dataframe

LoadFc <- function(path){
  file <- read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)  # Načtení souboru
  exprs <- exprs(file)  # Extrakce surových dat
  df <- as.data.frame(exprs)  # Převedení data frame
  return(df)
}

# Apply LoadFc
# Copy and paste fill path of each individual .fcs file

neg_exp19 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250523-exp19_NEG_.fcs")
background_647_exp19 <- LoadFc("./figure5de+supplfigs5de/siRNA-HPGclick-20250523-exp19_NEGSEK_1_.fcs")
sicon_sicon_exp19 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-202550523-exp19_siCONsiCON_1_.fcs")
sirna_sicon_exp19 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-202550523-exp19_siRNAsiCON_1_.fcs")
sirna_sirna_exp19 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-202550523-exp19_siRNAsiRNA_1_.fcs")

neg_exp20 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250606-exp20_NEG_.fcs")
background_647_exp20 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250606-exp20_BACK647_.fcs")
sicon_sicon_exp20 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250606-exp20_siCONsiCON_.fcs")
sirna_sicon_exp20 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250606-exp20_siRNAsiCON_.fcs")
sirna_sirna_exp20 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250606-exp20_siRNAsiRNA_.fcs")

neg_exp22 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250613-exp22_NEG_.fcs")
background_647_exp22 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250613-exp22_BACK647_.fcs")
sicon_sicon_exp22 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250613-exp22_siCONsiCON_.fcs")
sirna_sicon_1_exp22 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250613-exp22_siRNAsiCON_1_.fcs")
sirna_sicon_2_exp22 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250613-exp22_siRNAsiCON_2_.fcs")
sirna_sirna_exp22 <- LoadFc ("./figure5de+supplfigs5de/siRNA-HPGclick-20250613-exp22_siRNAsiRNA_.fcs")
sirna_sicon_exp22 <- rbind(sirna_sicon_1_exp22, sirna_sicon_2_exp22)

#------------------------------- exp19 data processing---------------------------------

# Visualisation of single-cell gating based on FSC-A vs FSC-H
# Red lines show gate used to exclude debris and doublets

ggplot(sirna_sirna_exp19, aes(x = `FSC-A`, y = `FSC-H`)) +
  geom_point(alpha = 0.1, size = 1) + 
  geom_vline(xintercept = 100000, color = "red")+
  geom_vline(xintercept = 260000, color = "red")+
  geom_hline(yintercept = 50000, color = "red")+
  geom_hline(yintercept = 150000, color = "red")+
  scale_x_continuous(labels = comma)+
  scale_y_continuous(labels = comma)+
  theme_classic()

# Define FITC gates
# Events are filtrated as single-cell gate & FITC > 0 for standard deviation calculation
single_cells_neg_exp19 <- sirna_sirna_exp19 %>% filter (`FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000 & `FITC-H` > 0)

# FITC gate threshold based on negative population (control cells)
# Gate is defined as mean (FICT-H) plus five times standard deviation (less strict gate)
fitc_gate_exp19 <- round(10^(mean(log10(single_cells_neg_exp19$`FITC-H`)) + 5 * sd(log10(single_cells_neg_exp19$`FITC-H`))))

# Use FITC gates as variables for mixed sample
fitc_mix_con_exp19 <- fitc_gate_exp19
fitc_mix_sirna_exp19 <- fitc_gate_exp19

# FITC gates visualisation
ggplot(single_cells_neg_exp19, aes(x = `FSC-H`, y = `FITC-H`))+
  geom_point(size = 1, alpha = 0.3)+
  geom_hline(yintercept = fitc_gate_exp19)+
  scale_y_log10()+
  theme_classic()

# Create datasets 
# Filter events based on single-cell gating
# Add new column to label each dataset
# Add log-transformed Alexa 647 values for further analysis
control <- sicon_sicon_exp19 %>% 
  filter(`FITC-H`> 0 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "control")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

sirna <- sirna_sirna_exp19  %>% 
  filter(`FITC-H`> 0 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000) %>%
  mutate(dataset = "sirna")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

mix_control <- sirna_sicon_exp19 %>% 
  filter(`FITC-H`>= fitc_mix_con_exp19 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "mix_control")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

mix_sirna <- sirna_sicon_exp19 %>% 
  filter(`FITC-H`< fitc_mix_sirna_exp19 & `FITC-H`>0 &`Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "mix_sirna")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

# Combine all datasets into one
# Add column to specify experimnet number
fitc_alexa_exp19 <- rbind(control, sirna, mix_control, mix_sirna)%>%
  mutate(experiment = "exp19")


#------------------------------- exp20 data processing---------------------------------

# Visualisation of single-cell gating based on FSC-A vs FSC-H
# Red lines show gate used to exclude debris and doublets

ggplot(sirna_sirna_exp20, aes(x = `FSC-A`, y = `FSC-H`)) +
  geom_point(alpha = 0.1, size = 1) + 
  geom_vline(xintercept = 100000, color = "red")+
  geom_vline(xintercept = 260000, color = "red")+
  geom_hline(yintercept = 50000, color = "red")+
  geom_hline(yintercept = 150000, color = "red")+
  scale_x_continuous(labels = comma)+
  scale_y_continuous(labels = comma)+
  theme_classic()

# Define FITC gates
# Events are filtrated as single-cell gate & FITC > 0 for standard deviation calculation
single_cells_neg_exp20 <- sirna_sirna_exp20 %>% filter (`FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000 & `FITC-H` > 0)

# FITC gate threshold based on negative population (control cells)
# Gate is defined as mean (FICT-H) plus five times standard deviation (less strict gate)
fitc_gate_exp20 <- round(10^(mean(log10(single_cells_neg_exp20$`FITC-H`)) + 5 * sd(log10(single_cells_neg_exp20$`FITC-H`))))

# Use FITC gates as variables for mixed sample
fitc_mix_con_exp20 <- fitc_gate_exp20
fitc_mix_sirna_exp20 <- fitc_gate_exp20

# FITC gates visualisation
ggplot(single_cells_neg_exp20, aes(x = `FSC-H`, y = `FITC-H`))+
  geom_point(size = 1, alpha = 0.3)+
  geom_hline(yintercept = fitc_gate_exp20)+
  scale_y_log10()+
  theme_classic()

# Create datasets 
# Filter events based on single-cell gating
# Add new column to label each dataset
# Add log-transformed Alexa 647 values for further analysis
control <- sicon_sicon_exp20 %>% 
  filter(`FITC-H`> 0 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "control")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

sirna <- sirna_sirna_exp20  %>% 
  filter(`FITC-H`> 0 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000) %>%
  mutate(dataset = "sirna")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

mix_control <- sirna_sicon_exp20 %>% 
  filter(`FITC-H`>= fitc_mix_con_exp20 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "mix_control")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

mix_sirna <- sirna_sicon_exp20 %>% 
  filter(`FITC-H`< fitc_mix_sirna_exp20 & `FITC-H`>0 &`Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "mix_sirna")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

# Combine all datasets into one
# Add column to specify experimnet number
fitc_alexa_exp20 <- rbind(control, sirna, mix_control, mix_sirna)%>%
  mutate(experiment = "exp20")

#------------------------------- exp22 data processing---------------------------------

# Visualisation of single-cell gating based on FSC-A vs FSC-H
# Red lines show gate used to exclude debris and doublets

ggplot(sirna_sirna_exp22, aes(x = `FSC-A`, y = `FSC-H`)) +
  geom_point(alpha = 0.1, size = 1) + 
  geom_vline(xintercept = 100000, color = "red")+
  geom_vline(xintercept = 260000, color = "red")+
  geom_hline(yintercept = 50000, color = "red")+
  geom_hline(yintercept = 150000, color = "red")+
  scale_x_continuous(labels = comma)+
  scale_y_continuous(labels = comma)+
  theme_classic()

# Define FITC gates
# Events are filtrated as single-cell gate & FITC > 0 for standard deviation calculation
single_cells_neg_exp22 <- sirna_sirna_exp22 %>% filter (`FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000 & `FITC-H` > 0)

# FITC gate threshold based on negative population (control cells)
# Gate is defined as mean (FICT-H) plus five times standard deviation (less strict gate)
fitc_gate_exp22 <- round(10^(mean(log10(single_cells_neg_exp22$`FITC-H`)) + 5 * sd(log10(single_cells_neg_exp22$`FITC-H`))))

# Use FITC gates as variables for mixed sample
fitc_mix_con_exp22 <- fitc_gate_exp22
fitc_mix_sirna_exp22 <- fitc_gate_exp22

# FITC gates visualisation
ggplot(single_cells_neg_exp22, aes(x = `FSC-H`, y = `FITC-H`))+
  geom_point(size = 1, alpha = 0.3)+
  geom_hline(yintercept = fitc_gate_exp22)+
  scale_y_log10()+
  theme_classic()

# Create datasets 
# Filter events based on single-cell gating
# Add new column to label each dataset
# Add log-transformed Alexa 647 values for further analysis
control <- sicon_sicon_exp22 %>% 
  filter(`FITC-H`> 0 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "control")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

sirna <- sirna_sirna_exp22  %>% 
  filter(`FITC-H`> 0 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000) %>%
  mutate(dataset = "sirna")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

mix_control <- sirna_sicon_exp22 %>% 
  filter(`FITC-H`> fitc_mix_con_exp22 & `Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "mix_control")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

mix_sirna <- sirna_sicon_exp22 %>% 
  filter(`FITC-H`< fitc_mix_sirna_exp22 & `FITC-H`>0 &`Alexa 647-H`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "mix_sirna")%>%
  mutate(log10_alexa = log10(`Alexa 647-H`))

# Combine all datasets into one
# Add column to specify experiment number
fitc_alexa_exp22 <- rbind(control, sirna, mix_control, mix_sirna)%>%
  mutate(experiment = "exp22")

# Combine all experiments into one dataset
a <- rbind(fitc_alexa_exp19, fitc_alexa_exp20, fitc_alexa_exp22)
#export all data into a csv file
write.csv(a, "./figure5de+supplfigs5de/colour-reversed-exp23_27_28_one_gate.csv", row.names = FALSE)

#------------------------------- Experiments analysis ---------------------------------

# Choose number of experiment for analysis

exp <- a$experiment == "exp22"

# Variable- number of bins used for analysis
bins <- 80

# Geometric mean values for each sample
geom_mean_con <- exp(mean(log(a$`Alexa 647-H`[exp & a$dataset == "control"])))
geom_mean_mix_control <- exp(mean(log(a$`Alexa 647-H`[exp & a$dataset == "mix_control" & a$`Alexa 647-H` > 0])))
geom_mean_mix_sirna <- exp(mean(log(a$`Alexa 647-H`[exp & a$dataset == "mix_sirna" & a$`Alexa 647-H` > 0])))
geom_mean_sirna <- exp(mean(log(a$`Alexa 647-H`[exp & a$dataset == "sirna" & a$`Alexa 647-H` > 0])))

# Calculate histogram data from log-trasfrormed Alexa 647 values
control_hist <- hist(a$log10_alexa[exp & (a$dataset == "control")], breaks = bins, plot = FALSE)
mix_sirna_hist <- hist(a$log10_alexa[exp & a$dataset == "mix_sirna"], breaks = bins, plot = FALSE)
mix_control_hist <- hist(a$log10_alexa[exp & a$dataset == "mix_control"], breaks = bins, plot = FALSE)
sirna_hist <- hist(a$log10_alexa[exp & a$dataset == "sirna"], breaks = bins, plot = FALSE)

# Convert calculated histogram data into dataframe
# Define function
hist_to_df <- function(h, dataset_name) {
  data.frame(
    x = h$mids,
    count = h$counts,
    density = h$density) %>%
    mutate(dataset = dataset_name)}

# Apply function
control_df      <- hist_to_df(control_hist, "control")
sirna_df        <- hist_to_df(sirna_hist, "sirna")
mix_control_df  <- hist_to_df(mix_control_hist, "mix_control")
mix_sirna_df    <- hist_to_df(mix_sirna_hist, "mix_sirna")

# Combine dataframes info one for further analysis
all_hist_df <- bind_rows(control_df, sirna_df, mix_control_df, mix_sirna_df)

#------------------------------- Final plots ---------------------------------

# Overall translation histogram plot

ggplot(all_hist_df, aes(x = x, y = count)) +
  geom_col(aes(fill = dataset), width = 0.045, alpha = 0.2, position = "identity") +
  geom_line(aes(colour = dataset),  size = 0.4) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("10", "100", "1000", "10000"))+
  theme_classic(base_size = 7)+
  scale_color_manual(values = c(
    "control" = "#60AF1A",
    "sirna" = "#0B7D9E",
    "mix_control" = "#A8DC79",
    "mix_sirna" = "#00C0BA"))+
  scale_fill_manual(values = c(
    "control" = "#60AF1A",
    "sirna" = "#0B7D9E",
    "mix_control" = "#A8DC79",
    "mix_sirna" = "#00C0BA"))+
  theme(legend.position = "none",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank())

# Define experiment number for plotting purposes
exp_num <- "22"
fitc_gate <- get(paste0("fitc_gate_exp", exp_num))
exp <- a[a$experiment == paste0("exp", exp_num), ]

# Scatter plot of FITC-Alexa 647 distribution

ggplot(data = exp, aes(y = `Alexa 647-H`, x = `FITC-H`, colour = dataset))+
  geom_point(size = 1, alpha = 0.3)+
  scale_x_log10(labels = comma)+
  scale_y_log10(labels = comma)+
  scale_color_manual(values = c(
    "control" = "#60AF1A",
    "sirna" = "#0B7D9E",
    "mix_control" = "#A8DC79",
    "mix_sirna" = "#00C0BA"))+
  theme_classic(base_size = 7)+
  theme(legend.position = "none",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5),
        #axis.text = element_blank(),
        axis.title = element_blank())

# Histogram plot of FITC distribution

ggplot(data = exp, aes(x = `FITC-H`, fill = dataset)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 120) +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma)+
  theme_classic(base_size = 7) +
  geom_vline(xintercept = fitc_gate, linewidth = 0.3)+
  scale_fill_manual(values = c(
    "control" = "#60AF1A",
    "sirna" = "#0B7D9E",
    "mix_control" = "#A8DC79",
    "mix_sirna" = "#00C0BA"))+
  theme(legend.position = "none",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank())


















