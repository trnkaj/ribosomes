library(tidyverse)

#230323
cellglo230323 <- X230323_CellGlo_gemcit %>% 
  pivot_longer(2:13,
               names_to = "column", 
               values_to = "luminescence") %>% 
  rename(row = "<>")

cellglo230323$treatment <- c(rep("control_gemcit", 3), 0, rep("blank", 3), rep(0, 5), rep("control_gemcit", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("triton", 3), rep(0, 9)) 
cellglo230323 <- cellglo230323 %>% filter(treatment !=0)
cellglo230323$concentration <- c(rep(0, 9), rep(0.01, 3), rep(0.1, 3), rep(1, 3), rep(10, 3), rep(100, 3), rep(0, 3))
cellglo230323$date <- c(rep("230323", 27))
cellglo230323 <- cellglo230323 %>% relocate(date, .before = row)
cellglo230323$number_of_cells <- c(rep(300, 27))

blank230323 <- cellglo230323 %>% 
  select(treatment, luminescence) %>% 
  filter(treatment == "blank")
blank230323_value <- mean(blank230323$luminescence)
cellglo230323$lum_minus_blank <- cellglo230323$luminescence - blank230323_value

control230323 <- cellglo230323 %>% 
  select(treatment, lum_minus_blank) %>% 
  filter(treatment == "control_gemcit")
control230323_value <- mean(control230323$lum_minus_blank)

#230331
cellglo230331 <- X230331_CellGlo_gemcit %>% 
  pivot_longer(2:13,
               names_to = "column", 
               values_to = "luminescence") %>% 
  rename(row = "<>")

cellglo230331$treatment <-  c(rep("blank", 4), rep(0, 8), rep("control_gemcit", 4), rep(0, 8), rep("gemcitabine", 4), rep(0, 8), rep("gemcitabine", 4), rep(0, 8), rep("gemcitabine", 4), rep(0, 8), rep("gemcitabine", 4), rep(0, 8), rep("gemcitabine", 4), rep(0, 8), rep("triton", 4), rep(0, 8))
cellglo230331 <- cellglo230331 %>% filter(treatment !=0)
cellglo230331$concentration <- c(rep(0,8), rep(0.01, 4), rep(0.1, 4), rep(1, 4), rep(10, 4), rep(100, 4), rep(0, 4))
cellglo230331$date <- c(rep("230331", 32))
cellglo230331 <- cellglo230331 %>% relocate(date, .before = row)
cellglo230331$number_of_cells <- c(rep(300, 32))

blank230331 <- cellglo230331 %>% 
  select(treatment, luminescence) %>% 
  filter(treatment == "blank")
blank230331_value <- mean(blank230331$luminescence)

cellglo230331$lum_minus_blank <- cellglo230331$luminescence - blank230331_value

control230331 <- cellglo230331 %>% 
  select(treatment, lum_minus_blank) %>% 
  filter(treatment == "control_gemcit")
control230331_value <- mean(control230331$lum_minus_blank)

#230409
cellglo230409 <- X230409_CellGlo_gemcit %>% 
  pivot_longer(2:13,
               names_to = "column", 
               values_to = "luminescence") %>% 
  rename(row = "<>")

cellglo230409$treatment <-  c(rep("blank", 3), rep(0, 10), rep("control_gemcit", 2), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep(0, 9), rep("gemcitabine", 3), rep("triton", 3), rep(0, 6)) 
cellglo230409 <- cellglo230409 %>% filter(treatment !=0)
cellglo230409$concentration <- c(rep(0,5), rep(0.01, 3), rep(0.1, 3), rep(1, 3), rep(10, 3), rep(100, 3), rep(200, 3), rep(0, 3))
cellglo230409$date <- c(rep("230409", 26))
cellglo230409 <- cellglo230409 %>% relocate(date, .before = row)
cellglo230409$number_of_cells <- c(rep(300, 26))

blank230409 <- cellglo230409 %>% 
  select(treatment, luminescence) %>% 
  filter(treatment == "blank")
blank230409_value <- mean(blank230409$luminescence)

cellglo230409$lum_minus_blank <- cellglo230409$luminescence - blank230409_value

control230409 <- cellglo230409 %>% 
  select(treatment, lum_minus_blank) %>% 
  filter(treatment == "control_gemcit")
control230409_value <- mean(control230409$lum_minus_blank)

#merge
cellglo_data <- merge(cellglo230323, cellglo230331, all.x = TRUE, all.y = TRUE)
cellglo_data <- merge(cellglo_data, cellglo230409, all.x = TRUE, all.y = TRUE)

control_gemcit_total <- mean(c(control230323_value, control230331_value, control230409_value))


library(drc)

drm_gemcit <- cellglo_data %>% 
  dplyr::filter(treatment %in% c("control_gemcit", "gemcitabine")) %>% 
  dplyr::select(lum_minus_blank, concentration)
drm_gemcit$adj_lum <- drm_gemcit$lum_minus_blank/1000

model_gemcit <- drm(drm_gemcit$adj_lum ~ drm_gemcit$concentration,
                    fct = LL.4(fixed = c(NA, 0, control_gemcit_total/1000, NA), names = c("slope", "min", "max", "EC50")))

summary(model_gemcit)

plot(model_gemcit, type = "all", 
     main = "Dose-response curve", 
     xlab = "gemcitabine (μM)", 
     ylab = "luminescence (in thousands)", 
     ylim = c(0,500),
     cex.lab = 1.4,
     pch = 19,
     col = rgb(0, 0, 0, 0.5))
text(10, 400, "IC50 = 1.6 μM", cex = 1.4)

