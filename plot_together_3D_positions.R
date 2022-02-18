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

set.seed(2022-02-16)

# read all files with selected name ending and extension
filenames <- Sys.glob(file.path("positions*.csv"))
filenames
# make data frame to collect data into through the loop
data_merge_pos <- data.frame()
#loop through all files
for (i in filenames){  
  x <- read.table(i, sep = ";",header = TRUE)
  name <- gsub("positions_Imaris_","",i)
  name <- gsub(".csv","",name)
  x$Condition <- name
  print(head(x))
  data_merge_pos <- rbind(data_merge_pos, x)
}
#To order the conditions manually
data_merge_pos$Condition <- factor(data_merge_pos$Condition, levels=c("B6", "TKO"))
#plot histograms for collected filtered results
cbp <- c("#000000", "#967491")
cols <- c("Moving"="#adb5bd", "Stationary"="#000000")

data_merge_pos_orig <- data_merge_pos
data_merge_pos <- na.omit(data_merge_pos)

### Count points with respective phase
p_stop_counts_1 <- ggplot(data_merge_pos, aes(x = factor(Condition)))+
  #facet_grid(~condition) +
  geom_bar(aes(fill = Movement, color = Condition), position="fill", width=0.5) +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab("Percent of respective points") +
  xlab("")
p_stop_counts_1

ggsave(p_stop_counts_1, filename = paste("stop_percent.png", sep=""), type = "cairo", width = 4, height = 5)

p_stop_counts_2 <- ggplot(data_merge_pos, aes(x = factor(Condition)))+
  #facet_grid(~condition) +
  geom_bar(aes(fill = Movement, color = Condition), width=0.5) +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab("Count of respective points") +
  xlab("")
p_stop_counts_2

# Time spent stationary vs moving
p_time_phase_1 <- ggplot(data_merge_pos, aes(x = Condition, y = delta_time))+
  #facet_grid(~Condition) +
  geom_bar(stat = "identity", aes(fill = Movement), position="fill", width=0.5) +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab("Percent of time spent in a phase") +
  xlab("")
p_time_phase_1

p_time_phase_2 <- ggplot(data_merge_pos, aes(x = Condition, y = delta_time))+
  #facet_grid(~Condition) +
  geom_bar(stat = "identity", aes(fill = Movement), width=0.5) +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab("Time spent in a phase") +
  xlab("")
p_time_phase_2

ggsave(p_time_phase_2, filename = paste("stop_time.png", sep=""), type = "cairo", width = 4, height = 5)

all_plot <- ggarrange(p_stop_counts_2, p_stop_counts_1, p_time_phase_2, p_time_phase_1,
                      labels = c("A", "B", "C", "D"),
                      ncol = 4, nrow = 1,
                      common.legend = TRUE, legend = "bottom")
all_plot
ggsave(all_plot, filename = paste("stop_params_together.png", sep=""), type = "cairo", width = 12, height = 6)


# Try no aggregation for velocity comparison
data_merge_pos$TrackID <- as.factor(data_merge_pos$TrackID)
data_merge_pos$Condition <- as.factor(data_merge_pos$Condition)
data_merge_pos$Movement <- as.factor(data_merge_pos$Movement)
data_merge_pos_mov <- data_merge_pos[data_merge_pos$Movement == "Moving",]
df_for_plot_mean <- data_merge_pos_mov %>% group_by(Name) %>% dplyr::summarise(mean(velocity))
colnames(df_for_plot_mean) <- c("Name", "Velocity")
df_for_plot_mean <- merge(df_for_plot_mean, unique(data_merge_pos_mov[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p1 <- ggplot(df_for_plot_mean,aes(x=Condition,y=Velocity))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=0.5, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$Velocity, n = 6)) + #limits = c(0,1),
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
p1
ggsave(p1, filename = paste("Velocity_moving_means_per_dendrites.png", sep=""), type = "cairo", width = 4, height = 5)




# Fun
p <- ggplot(data_merge_pos[data_merge_pos$Movement == "Moving",],aes(x=Position.X,y=Position.Y))+
  geom_point(aes(color=Condition)) +
  geom_point(shape=21, data = data_merge_pos[data_merge_pos$Movement == "Stationary",], size=3, aes(color=Condition), fill="red") +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cols) +
  theme_classic() 
p


df_list <- list()
df_list <- split(data_merge_pos_orig, data_merge_pos_orig$TrackID)
for(i in 1:length(df_list)){
  start_x <- df_list[[i]]$Position.X[[1]]
  start_y <- df_list[[i]]$Position.Y[[1]]
  start_z <- df_list[[i]]$Position.Z[[1]]
  
  x_turn <- df_list[[i]]$Position.X[[1]]
  
  
  for (j in 1:nrow(df_list[[i]])) {
    df_list[[i]]$Position.X[[j]] <- df_list[[i]]$Position.X[[j]]-start_x
    df_list[[i]]$Position.Y[[j]] <- df_list[[i]]$Position.Y[[j]]-start_y
    df_list[[i]]$Position.Z[[j]] <- df_list[[i]]$Position.Z[[j]]-start_z
  }
}
df_transform <- do.call(rbind, df_list)

p <- ggplot(df_transform[df_transform$Movement == "Moving",],aes(x=Position.X,y=Position.Y))+
  geom_point(aes(color=Condition)) +
  geom_point(shape=21, data = df_transform[df_transform$Movement == "Stationary",], size=3, aes(color=Condition), fill="red") +
  scale_color_manual(values = cbp) +
  scale_fill_manual(values = cols) +
  theme_classic() 
p


p_3D <- plot_ly(data_merge_pos_orig, x=~Position.X, y=~Position.Y, z=~Position.Z, 
             type="scatter3d", mode="markers",
             color = ~Condition, colors = cbp)

htmlwidgets::saveWidget(as_widget(p_3D), "all_tracks_together.html")
