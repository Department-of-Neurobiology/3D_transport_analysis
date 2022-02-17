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
# Aggregate and plot mean velocities per cell
p1 <- ggplot(data_merge_info,aes(x=Condition,y=Velocity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=0.5, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Velocity by TrackID") +
  scale_y_continuous(breaks = pretty(data_merge_info$Velocity, n = 6)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
p1
ggsave(p1, filename = paste("Velocity.png", sep=""), type = "cairo", width = 4, height = 5)


# Aggregate and plot mean velocities per cell
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Velocity))
colnames(df_for_plot_mean) <- c("Name", "Velocity")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p1 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Velocity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Velocity by dendrite") +
  scale_y_continuous(breaks = pretty(data_merge_info$Velocity, n = 10)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p1, filename = paste("mean_velocity_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)


# Aggregate and plot mean Velocity per dendrite
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Track.Displacement.Length))
df_for_plot_mean_2 <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Track.Duration))
colnames(df_for_plot_mean) <- c("Name", "Track.Displacement.Length")
colnames(df_for_plot_mean_2) <- c("Name", "Track.Duration")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
df_for_plot_mean <- merge(df_for_plot_mean, df_for_plot_mean_2, by ="Name", all.x = TRUE)
df_for_plot_mean$Velocity <- df_for_plot_mean$Track.Displacement.Length/df_for_plot_mean$Track.Duration

p2 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Velocity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("mean(Track.Displacement.Length)/mean(Track.Duration) per dendrite") +
  scale_y_continuous(limits = c(0,3),breaks = pretty(data_merge_info$Velocity, n = 4)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p2, filename = paste("meanTrack.Displacement.Length_by_meanTrack.Duration_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)


plot <- ggplot(data_merge_info,aes(x=Condition,y=Processivity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=0.5, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  scale_y_continuous(breaks = pretty(data_merge_info$Processivity, n = 6)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")



# Aggregate and plot mean processivity per dendrite
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Processivity))
colnames(df_for_plot_mean) <- c("Name", "Processivity")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p3 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Processivity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Processivity by dendrite") +
  scale_y_continuous(limits = c(0,1.2),breaks = pretty(data_merge_info$Processivity, n = 4)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p3, filename = paste("mean_processivity_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)


# Aggregate and plot mean processivity per dendrite
df_for_plot_mean <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Track.Displacement.Length))
df_for_plot_mean_2 <- data_merge_info %>% group_by(Name) %>% dplyr::summarise(mean(Track.Length))
colnames(df_for_plot_mean) <- c("Name", "Track.Displacement.Length")
colnames(df_for_plot_mean_2) <- c("Name", "Track.Length")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_info[c("Name", "Condition")]), by ="Name", all.x = TRUE)
df_for_plot_mean <- merge(df_for_plot_mean, df_for_plot_mean_2, by ="Name", all.x = TRUE)
df_for_plot_mean$Processivity <- df_for_plot_mean$Track.Displacement.Length/df_for_plot_mean$Track.Length

p4 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Processivity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("mean(Track.Displacement.Length)/mean(Track.Length) per dendrite") +
  scale_y_continuous(limits = c(0,1.2),breaks = pretty(data_merge_info$Processivity, n = 4)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
ggsave(p4, filename = paste("meanTrack.Displacement.Length_by_meanTrack.Length_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)

dev.off()

library(ggpubr)
all_plot <- ggarrange(p1, p2, p3, p4, 
                      labels = c("A", "B", "C", "D"),
                      ncol = 4, nrow = 1,
                      common.legend = TRUE, legend = "bottom")
all_plot
ggsave(all_plot, filename = paste("mean_calculation_comparison_per_dendrite.png", sep=""), type = "cairo", width = 12, height = 6)

# ggplot(data_merge_info,aes(x=Condition,y=Track.Straightness))+
#   geom_boxplot(width=0.5, color ="#555555") +
#   geom_jitter(width=0.2, aes(color=Condition)) +
#   scale_color_manual(values = cbp) +
#   theme_classic() +
#   scale_y_continuous(breaks = pretty(data_merge_info$Track.Straightness, n = 6)) +
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "B6")
