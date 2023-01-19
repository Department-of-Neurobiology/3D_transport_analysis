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
# library(coords)
library(mosaic)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

set.seed(2022-02-16)

# read all files with selected name ending and extension
filenames <- Sys.glob(file.path("df_displ*.csv"))
filenames
# make data frame to collect data into through the loop
data_merge_pos <- data.frame()
#loop through all files
for (i in filenames){  
  x <- read.table(i, sep = ";", header = TRUE)
  x$TrackID <- as.factor(x$TrackID)
  name <- gsub("df_displ_Imaris_","",i)
  name <- gsub(".csv","",name)
  x$Condition <- name
  print(head(x))
  data_merge_pos <- rbind(data_merge_pos, x)
}
data_merge_pos <- data_merge_pos[order(data_merge_pos$Time_s),]
data_merge_pos <- data_merge_pos[order(data_merge_pos$TrackID),]
#To order the conditions manually
data_merge_pos$Condition <- factor(data_merge_pos$Condition, levels=c("B6", "TKO"))
#plot histograms for collected filtered results
cbp <- c("#000000", "#967491")
#cols <- c("Moving"="#adb5bd", "Stationary"="#000000")
cols <- c("backward" = "#ff9955", "stationary" = "#000000", "forward" = "#87decd")

data_merge_pos_orig <- data_merge_pos
processive_IDs <- read.table("processive_IDs.csv", sep = ";", header = TRUE)
data_merge_pos <- data_merge_pos[data_merge_pos$TrackID %in% processive_IDs$x,]
data_merge_pos <- na.omit(data_merge_pos)
data_merge_pos <- droplevels(data_merge_pos)

# Problem with some tracks being renamed because of large numbers - tried to fix it in collecting script
# data_merge_pos$TrackID <- sub("^", "ID", data_merge_pos$TrackID)
# data_merge_pos$TrackID <- as.factor(data_merge_pos$TrackID)
# options(scipen=999)
# data_merge_pos$TrackID <- gsub("0000","",data_merge_pos$TrackID)


df_list <- list()
df_list <- split(data_merge_pos, data_merge_pos$TrackID)
for(i in 1:length(df_list)){
  #i=1
  df <- df_list[[i]]
  #start_mat <- as.matrix(df[c("Displacement.Delta.X","Displacement.Delta.Y","Displacement.Delta.Z")][1,])

  mat <-  df[c("Displacement.Delta.X", "Displacement.Delta.Y", "Displacement.Delta.Z")]
  mat$x0 <- 0
  mat$y0 <- 0
  mat$z0 <- 0
  mat <- as.matrix(mat)
  df_list[[i]]$distance <- apply(mat, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  
  for (j in 1:nrow(df_list[[i]])) {
      #print(j)
      if (j==1) {
        df_list[[i]]$distance_per_time[[j]] <- NA
        next
      } else {
        df_list[[i]]$distance_per_time[[j]] <- df_list[[i]]$distance[[j]]/(df_list[[i]]$Time_s[[j]]-df_list[[i]]$Time_s[[j-1]])
      }
    }
  
  
  last_vec <- c(as.matrix(df[c("Displacement.Delta.X","Displacement.Delta.Y","Displacement.Delta.Z")][nrow(df),]))
  
  sum_vec <-c(0,0,0)
  for (k in 1:nrow(df_list[[i]])) {
    sum_vec <- sum_vec + c(as.matrix(df[c("Displacement.Delta.X","Displacement.Delta.Y","Displacement.Delta.Z")][k,]))
  }
  
  for (j in 1:nrow(df_list[[i]])) {
    row_vec <- c(as.matrix(df[c("Displacement.Delta.X","Displacement.Delta.Y","Displacement.Delta.Z")][j,]))
    
    #project(y1, x1, type='length')
    sign_of_vector <- sign(project(row_vec ~ sum_vec))
  
    # Distance is calculated from vector length and it is where stops are calculated from
    # But can happen that the projection is zero and length not (perpendicular) => make sure to take those into account
    if(df_list[[i]]$distance[[j]] != 0) {
      if(sign_of_vector == 0) {
        print("Oh no!")
      }
    }
    df_list[[i]]$distance_sign[[j]] <- sign_of_vector*df_list[[i]]$distance[[j]]
    df_list[[i]]$distance_sign_per_time[[j]] <- (sign_of_vector*df_list[[i]]$distance_per_time[[j]])
    df_list[[i]]$motion[[j]] <- ifelse( df_list[[i]]$distance_sign_per_time[[j]] < -0.02, "backward", ifelse( df_list[[i]]$distance_sign_per_time[[j]] > 0.02, "forward", "stationary"))
  }
}
df_transform_1 <- do.call(rbind, df_list)
df_transform_1 <- df_transform_1[order(df_transform_1$Time_s),]
df_transform_1 <- df_transform_1[order(df_transform_1$TrackID),]

df_transform_1 <-  as.data.frame(lapply(df_transform_1, unlist))
df_transform_1$distance_sign <- as.numeric(df_transform_1$distance_sign)
df_transform_1$motion <- as.factor(df_transform_1$motion)
df_transform_1$disp_length <- NA

#create an empty list to store the split dataframes
df_list <- list()
#split the dataframe df_transform_1 by the column 'TrackID' and store it in the df_list
df_list <- split(df_transform_1, df_transform_1$TrackID)
#loop through the list of dataframes
for(i in 1:length(df_list)){
  #loop through each row of the current dataframe
  for (j in 1:nrow(df_list[[i]])) {
    #print(j)
    #check if the current row is the first row
    if (j==1) {
      #if it is the first row, assign NA to the 'disp_length' column
      df_list[[i]]$disp_length[[j]] <- NA
      next
    }
    #check if the current row is the second row
    if (j==2) {
      #if it is the second row, assign the value of 'distance_sign' to the 'disp_length' column
      df_list[[i]]$disp_length[[j]] <- df_list[[i]]$distance_sign[[j]]
    } else {
      #if it is not the first or second row, add the value of 'distance_sign' to the previous value in 'disp_length' column
      df_list[[i]]$disp_length[[j]] <- df_list[[i]]$disp_length[[j-1]] + df_list[[i]]$distance_sign[[j]]
    }
  }
}
#recombine the list of dataframes into a single dataframe
df_transform <- do.call(rbind, df_list)

df_transform <- na.omit(df_transform)
df_transform <- droplevels(df_transform)


p <- ggplot(df_transform, aes(x = disp_length, y = Time_s)) + #, text = paste("TrackID:",TrackID)
  #geom_path(data = df_for_plot, aes(x = `Displacement X Reference Frame`, y = Time, color = TrackID), linetype = "dashed") + #has to have color = TrackID to be shown as separate
  facet_grid(~Condition) +
  geom_path(aes(group=TrackID)) +
  geom_point(aes(color=motion), size=3) + #, shape=21
  scale_color_manual(values = cols) +
  #dark_theme_gray() +
  theme_classic() +
  ggtitle("Displacement representation of processive tracks") +
  #xlab("Displacement") +
  #ylab("Time") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks = element_line(colour="black"),
        axis.text.y = element_text(color = "black", size = 10)) +
  scale_y_reverse(breaks = seq(0, 200, by = 60)) #+ #limits = c(60, 0), breaks = seq(0, 60, by = 10)
  #scale_y_continuous() +
  #scale_x_continuous(limits = c(-9, 9), breaks = seq(-8, 8, by = 2)) #limits = c(-2, 4), breaks = seq(-2, 4, by = 2)
p
ggsave(p, filename = paste("displacement_plot_0.2processive_per_dendrite.png", sep=""), type = "cairo", width = 15, height = 10)


# # Density plot
# ggplot(df_transform[df_transform$motion != "stationary",], aes(x = distance_sign)) +
#   facet_wrap(~Condition, ncol = 1, strip.position = "bottom") +
#   geom_density(aes(fill = Condition, color = motion), alpha = 0.4) +
#   #geom_histogram(aes(fill = motion), alpha = 0.4,bins=500)+
#   #theme_classic() +
#   theme_classic() +
#   scale_color_manual(values=cols) +
#   scale_fill_manual(values=cbp) +
#   ggtitle("Trajectories, mobile fraction only") +
#   xlab("Delta_x") +
#   ylab("Density")
# 
ggplot(df_transform,aes(x=Condition))+
  geom_bar(width=0.5, aes(fill=motion), alpha=0.8, position="fill") +
  scale_fill_manual(values=cols) +
  ylab("Percentage of directional movement")

ggsave(p, filename = paste("displacement_plot_0.2processive_per_dendrite.png", sep=""), type = "cairo", width = 15, height = 10)

# 
# ggplot(df_transform, aes(x = as.numeric(distance_sign), fill = Condition, color = Condition))+
#   #facet_grid(~condition) +
#   #geom_bar(aes(fill = Condition, color = Condition), position="fill", width=0.5) +
#   geom_density(alpha=0.4, fill="transparent") +
#   scale_color_manual(values = cbp) +
#   scale_fill_manual(values = cbp) +
#   theme_classic() +
#   ylab("Density of length displacement") +
#   xlab("")
# 
# 
# plot_transform <- ggplot(df_transform[df_transform$Movement == "Moving",],aes(x=Position.X,y=Position.Y))+
#   geom_point(aes(color=Condition)) +
#   geom_path(aes(group = 1)) +
#   geom_point(shape=21, data = df_transform[df_transform$Movement == "Stationary",], size=3, aes(color=Condition), fill="red") +
#   scale_color_manual(values = cbp) +
#   scale_fill_manual(values = cols) +
#   theme_classic() 
# 
# 
# plot_wo_transformation
# 
# 
# ggarrange(plot_wo_transformation, plot_transform,
#                       labels = c("A", "B"),
#                       ncol = 1, nrow = 2,
#                       common.legend = TRUE, legend = "bottom")
# 
# plot_ly(df_transform[df_transform$Movement == "Moving",], x=~Position.X, y=~Position.Y, z=~Position.Z, 
#                 type="scatter3d", mode="markers",
#                 color = ~Condition, colors = cbp)
# 
# 
# p_3D <- plot_ly(data_merge_pos, x=~Displacement.Delta.X, y=~Displacement.Delta.Y, z=~Displacement.Delta.Z,
#              type="scatter3d", mode="markers+lines",
#              color = ~Condition, colors = cbp) 
#htmlwidgets::saveWidget(as_widget(p_3D), "all_tracks_together.html")



#count state fractions
df_list_3 <- list()
df_list_3 <- split(df_transform, df_transform$TrackID)

fractions_datalist <- list()
fractions_df <- create_empty_table(3,1)
changes_datalist <- list()
state_changes_df <- create_empty_table(13,1)
#for each trajectory creates data lists for further assembly into tables
for(i in 1:length(df_list_3)){
  #i=1
  #set row names later as trajectory IDs
  name <- as.character(df_list_3[[i]]$TrackID[1])
  #state changes for each trajectory
  retr_stat=0
  ant_stat=0
  retr_ant=0
  ant_retr=0
  stat_retr=0
  stat_ant=0
  
  #additional columns to count all changes to and from stops
  all_stat=0
  stat_all=0
  
  if (nrow(df_list_3[[i]])-1 > 1) { #24...
    #state percentage for each trajectory - IS IT REASONABLE? then further on they will be taken as mean for each cell
    fractions_df[1,] <- 100*summary(df_list_3[[i]]$motion)[1]/sum(summary(df_list_3[[i]]$motion))
    fractions_df[2,] <- 100*summary(df_list_3[[i]]$motion)[2]/sum(summary(df_list_3[[i]]$motion))
    fractions_df[3,] <- 100*summary(df_list_3[[i]]$motion)[3]/sum(summary(df_list_3[[i]]$motion))
    fractions_datalist[name] <-  as.list(fractions_df)
    for (j in 2:(nrow(df_list_3[[i]])-1)) {
      state_change <- ifelse(df_list_3[[i]]$motion[[j]]==df_list_3[[i]]$motion[[j-1]],"same","different")
      if (state_change=="different") {
        if (df_list_3[[i]]$motion[[j]]=="stationary" & df_list_3[[i]]$motion[[j-1]]=="backward") {
          retr_stat=retr_stat+1
        }
        if (df_list_3[[i]]$motion[[j]]=="stationary" & df_list_3[[i]]$motion[[j-1]]=="forward") {
          ant_stat=ant_stat+1
        }
        if (df_list_3[[i]]$motion[[j]]=="forward" & df_list_3[[i]]$motion[[j-1]]=="backward") {
          retr_ant=retr_ant+1
        }
        if (df_list_3[[i]]$motion[[j]]=="backward" & df_list_3[[i]]$motion[[j-1]]=="forward") {
          ant_retr=ant_retr+1
        }
        if (df_list_3[[i]]$motion[[j]]=="backward" & df_list_3[[i]]$motion[[j-1]]=="stationary") {
          stat_retr=stat_retr+1
        }
        if (df_list_3[[i]]$motion[[j]]=="forward" & df_list_3[[i]]$motion[[j-1]]=="stationary") {
          stat_ant=stat_ant+1
        }
        #optional: to count all changes to and from stops
        #if ((df_list_3[[i]]$motion[[j]]=="forward" | df_list_3[[i]]$motion[[j]]=="backward") & df_list_3[[i]]$motion[[j-1]]=="stationary") {
        #stat_all=stat_all+1
        #}
      }  
    }
    #relative to trajectory length, before was not divided, beware of the gaps time!=tracked_spots*0.2, can be higher, better take from Tibco table
    n <- (nrow(df_list_3[[i]]))
    if (n != 0) {
      track_duration <- df_list_3[[i]]$Time_s[[n]]-df_list_3[[i]]$Time_s[[1]]
      #print(track_duration)
      state_changes_df[1,] <- retr_stat/track_duration
      state_changes_df[2,] <- ant_stat/track_duration
      state_changes_df[3,] <- retr_ant/track_duration
      state_changes_df[4,] <- ant_retr/track_duration
      state_changes_df[5,] <- stat_retr/track_duration
      state_changes_df[6,] <- stat_ant/track_duration
      state_changes_df[7,] <- (retr_stat+ant_stat)/track_duration
      state_changes_df[8,] <- (stat_retr+stat_ant)/track_duration
      state_changes_df[9,] <- (retr_stat+retr_ant)/track_duration
      state_changes_df[10,] <- (ant_stat+ant_retr)/track_duration
      state_changes_df[11,] <- (stat_retr+stat_ant+retr_stat+retr_ant+ant_stat+ant_retr)/track_duration
      state_changes_df[12,] <- as.character(df_list_3[[i]]$Condition[1])
      state_changes_df[13,] <- as.character(df_list_3[[i]]$Name[1])
      changes_datalist[name] <-  as.list(state_changes_df)
    }
  }
}

assembled_state_changes <- as.data.frame(do.call(rbind, changes_datalist))
assembled_state_changes <- setNames(assembled_state_changes,c("RS","AS","RA","AR","SR","SA","all_to_stop","stop_to_all","retr_to_stop_ant","ant_to_stop_retr","all_changes","Condition", "Name"))
assembled_state_changes$Condition <- factor(assembled_state_changes$Condition, levels=c("B6", "TKO"))

assembled_state_changes[c("RS","AS","RA","AR","SR","SA","all_to_stop","stop_to_all","retr_to_stop_ant","ant_to_stop_retr","all_changes")] <- lapply(assembled_state_changes[c("RS","AS","RA","AR","SR","SA","all_to_stop","stop_to_all","retr_to_stop_ant","ant_to_stop_retr","all_changes")],as.numeric)

ggplot(assembled_state_changes,aes(x=Condition,y=all_changes))+
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  #ylab("all_to_stop") +
  #scale_y_continuous(breaks = pretty(assembled_state_changes$all_to_stop, n = 5)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")

df_for_plot_mean <- assembled_state_changes %>% group_by(Name) %>% dplyr::summarise(mean(all_changes))
colnames(df_for_plot_mean) <- c("Name", "all_changes")
df_for_plot_mean <- merge(df_for_plot_mean, unique(assembled_state_changes[c("Name", "Condition")]), by ="Name", all.x = TRUE)
p_all_changes <- ggplot(df_for_plot_mean,aes(x=Condition,y=all_changes)) + #, text = Name
  #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), width=0.5, color ="#555555") +
  geom_boxplot(width=0.5, color ="#555555") +
  geom_jitter(size=2, width=0.3, aes(color=Condition)) +
  scale_color_manual(values = cbp) +
  theme_classic() +
  ylab("Number of all direction changes by dendrite") +
  scale_y_continuous(breaks = pretty(df_for_plot_mean$all_changes, n = 5)) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "B6")
p_all_changes

ggsave(p_all_changes, filename = paste("all_changes_0.2processive_per_dendrite.png", sep=""), type = "cairo", width = 4, height = 5)

