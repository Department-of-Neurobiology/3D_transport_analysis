###############################################
##Nataliya Trushina, 2022##
##Analysis of Kif1A 3D trajectories##
#Requirements: csv tables with collected information
###############################################

# Libraries
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(plotly)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read all files with selected name ending and extension
filenames <- Sys.glob(file.path("imaris_collect*.csv"))
filenames
# make data frame to collect data into through the loop
data_merge_info <- data.frame()
#loop through all files
for (i in filenames){  
  x <- read.table(i, sep = ";",header = TRUE)
  name <- gsub("imaris_collect_Imaris_","",i)
  name <- gsub(".csv","",name)
  x$Condition <- name
  print(head(x))
  data_merge_info <- rbind(data_merge_info, x)
}
# Order the conditions manually
data_merge_info$Condition <- factor(data_merge_info$Condition, levels=c("B6", "TKO"))
data_merge_info$Name <- as.factor(data_merge_info$Name)

### Plots collected results ###
cbp <- c("#000000", "#967491")

set.seed(2022-02-16)

# # Aggregate and plot mean velocities per track
# p <- ggplot(data_merge_info,aes(x=Condition,y=Velocity))+
#   geom_boxplot(width=0.5, color ="#555555") +
#   geom_jitter(size=0.5, width=0.3, aes(color=Condition)) +
#   scale_color_manual(values = cbp) +
#   theme_classic() +
#   ylab("Velocity by TrackID") +
#   scale_y_continuous(breaks = pretty(data_merge_info$Velocity, n = 6)) +
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "B6")
# p
# ggsave(p1, filename = paste("Velocity.png", sep=""), type = "cairo", width = 4, height = 5)

# Aggregate and plot mean velocities per cell
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Velocity))
colnames(df_for_plot_mean) <- c("Name", "Velocity")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p1 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Velocity)) + #, text = Name
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Velocity by dendrite") +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$Velocity, n = 5)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p1, filename = paste("mean_velocity_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)
write.table(df_for_plot_mean, "mean_velocities_per_dendrite.csv", sep = ";",dec = '.', row.names = FALSE, col.names = TRUE)

# plot <- ggplot(data_merge_info,aes(x=Condition,y=Processivity))+
#   geom_boxplot(width=0.5, color ="#555555") +
#   geom_jitter(size=0.5, width=0.3, aes(color=Condition)) +
#   scale_color_manual(values = cbp) +
#   theme_classic() +
#   scale_y_continuous(breaks = pretty(data_merge_info$Processivity, n = 6)) +
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "B6")

# Aggregate and plot mean processivity per dendrite
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Processivity))
colnames(df_for_plot_mean) <- c("Name", "Processivity")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p2 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Processivity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Processivity by dendrite") +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$Processivity, n = 5)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p2, filename = paste("mean_processivity_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)
write.table(df_for_plot_mean, "mean_processivities_per_dendrite.csv", sep = ";",dec = '.', row.names = FALSE, col.names = TRUE)

# Aggregate and plot mean speed per dendrite
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Speed))
colnames(df_for_plot_mean) <- c("Name", "Speed")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p3 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Speed))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Speed by dendrite") +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$Speed, n = 5)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p3, filename = paste("mean_speed_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)
write.table(df_for_plot_mean, "mean_speed_per_dendrite.csv", sep = ";",dec = '.', row.names = FALSE, col.names = TRUE)


dev.off()

all_plot <- ggarrange(p3, p1, p2,
                      labels = c("A", "B", "C"),
                      ncol = 3, nrow = 1,
                      common.legend = TRUE, legend = "bottom")
all_plot
ggsave(all_plot, filename = paste("mean_params_per_dendrite.png", sep=""), type = "cairo", width = 12, height = 6)


# my_plot <- ggplotly(p1, tooltip = "Name")
# my_plot
# htmlwidgets::saveWidget(as_widget(my_plot), "mean_velocity_per_dendrite.png.html")

# ggplot(data_merge_info,aes(x=Condition,y=Track.Straightness))+
#   geom_boxplot(width=0.5, color ="#555555") +
#   geom_jitter(width=0.2, aes(color=Condition)) +
#   scale_color_manual(values = cbp) +
#   theme_classic() +
#   scale_y_continuous(breaks = pretty(data_merge_info$Track.Straightness, n = 6)) +
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "B6")


### Stop analysis ### - no difference
data_merge_info_real_stops <- data_merge_info[data_merge_info$Processivity >= 0.2,]
df_for_plot_mean <- data_merge_info_real_stops %>% group_by(Name) %>% dplyr::summarise(mean(Stop_count))
colnames(df_for_plot_mean) <- c("Name", "Stop_count")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p_stop <- ggplot(df_for_plot_mean,aes(x=Condition,y=Stop_count)) + #, text = Name
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), width=0.5, color ="#555555") +
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Real stop count by dendrite") +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$Stop_count, n = 10)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
p_stop
ggsave(p_stop, filename = paste("mean_stop_count_for_0.2processive_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)

### Stop time analysis ### - no difference for processive tracks, only if comparing all - more time spent stationary for TauKO
data_merge_info_stop_time <- data_merge_info[data_merge_info$Processivity >= 0.2,]
df_for_plot_mean <- data_merge_info_stop_time %>% group_by(Name) %>% dplyr::summarise(mean(Time_stopped))
colnames(df_for_plot_mean) <- c("Name", "Time_stopped")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p_stop_time <- ggplot(df_for_plot_mean,aes(x=Condition,y=Time_stopped)) + #, text = Name
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), width=0.5, color ="#555555") +
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Stationary time for processive tracks by dendrite [s]") +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$Time_stopped, n = 10)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
p_stop_time
ggsave(p_stop_time, filename = paste("mean_stop_time_for_0.2processive_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)
