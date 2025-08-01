library(flowCore)   # opening FSC files
library(ggplot2)
library(scales)     # in plots labels = comma
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

cycloheximide_exp3_1 <- LoadFc("./figure5b/siRNA-HPGclick-20250214-exp3_CHX_time_5_1_.fcs")
cycloheximide_exp3_2 <- LoadFc("./figure5b/siRNA-HPGclick-20250214-exp3_CHX_time_5_2_.fcs")
sirna_control_exp3_1 <- LoadFc("./figure5b/siRNA-HPGclick-20250214-exp3_CON_time_5_1_.fcs")
sirna_control_exp3_2 <- LoadFc("./figure5b/siRNA-HPGclick-20250214-exp3_CON_time_5_2_.fcs")

sirna_control_exp3 <- rbind(sirna_control_exp3_1, sirna_control_exp3_2)
cycloheximide_exp3 <- rbind(cycloheximide_exp3_1, cycloheximide_exp3_2)

#------------------------------- exp6 data processing ---------------------------------

# Visualisation of single-cell gating based on FSC-A vs FSC-H
# Red lines show gate used to exclude debris and doublets

ggplot(sirna_control_exp3, aes(x = `FSC-A`, y = `FSC-H`)) +
  geom_point(alpha = 0.1, size = 1) + 
  geom_vline(xintercept = 100000, color = "red")+
  geom_vline(xintercept = 260000, color = "red")+
  geom_hline(yintercept = 50000, color = "red")+
  geom_hline(yintercept = 170000, color = "red")+
  scale_x_continuous(labels = comma)+
  scale_y_continuous(labels = comma)+
  theme_classic()

# Create datasets 
# Filter events based on single-cell gating
# Add new column to label each dataset
# Add log-transformed FITC-H values for further analysis
control <- sirna_control_exp3 %>% 
  filter(`FITC-A`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 150000)%>%
  mutate(dataset = "control")%>%
  mutate(log10_fitc = log10(`FITC-A`))

cycloheximide <- cycloheximide_exp3 %>% 
  filter(`FITC-A`> 0 & `FSC-A` > 100000 & `FSC-A` <260000 & `FSC-H` > 50000 & `FSC-H` < 170000)%>%
  mutate(dataset = "cycloheximide")%>%
  mutate(log10_fitc = log10(`FITC-A`))

# Combine all datasets into one
# Add column to specify experiment number
t <- bind_rows(control,cycloheximide) %>%
  mutate(experiment = "exp6")

#------------------------------- Experiments analysis ---------------------------------

# Choose number of experiment for analysis
exp <- t$experiment == "exp6"

# Variable- number of bins used for analysis
bins <- 80

# Geometric mean values for each sample
geom_mean_con <- exp(mean(log(t$`FITC-A`[exp & t$dataset == "control"])))
geom_mean_chx <- exp(mean(log(t$`FITC-A`[exp & t$dataset == "cycloheximide"])))

# Calculate histogram data from log-trasfrormed FITC-H values
control_hist <- hist(t$log10_fitc[exp & (t$dataset == "control")], breaks = bins, plot = FALSE)
chx_hist <- hist(t$log10_fitc[exp & t$dataset == "cycloheximide"], breaks = bins, plot = FALSE)

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
chx_df          <- hist_to_df(chx_hist, "chx")

# Combine dataframes info one for further analysis
all_hist_df <- bind_rows(control_df, chx_df)

#------------------------------- Final plot ---------------------------------

# Overall translation histogram plot

ggplot(all_hist_df, aes(x = x, y = count)) +
  geom_col(aes(fill = dataset), width = 0.045, alpha = 0.2, position = "identity") +
  geom_line(aes(colour = dataset),  size = 0.4) +
  theme_classic(base_size = 7)+
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("10", "100", "1000", "10000"))+
  scale_color_manual(values = c(
    "control" = "#60AF1A",
    "sirna" = "#0B7D9E",
    "chx" = "grey"))+
  scale_fill_manual(values = c(
    "control" = "#60AF1A",
    "sirna" = "#0B7D9E",
    "chx" = "grey"))+
  theme(legend.position = "none",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5),
        axis.title = element_blank())




