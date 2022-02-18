###############################################
##Nataliya Trushina, 2022##
##3D transport analysis##

# Requirements: xlsx table with Displacement delta Reference Frame
# When the Excel files come from a different OneDrive account or from direct Imaris output, sometimes the files can not be opened with read_excel()
# open the file in Excel and save it again or use linux subsystem to automatically write them into new files to be able to open them
# cd /mnt/c/Users/NITru/OneDrive/Documents/PhD_work/Projects/Marina/3D_transport_analysis
# unoconv -f xlsx *.xls
# for dir in *; do [ -d "$dir" ] && unoconv -f xlsx "$dir"/*.xls; done

# Set the correct path to the directory with condition folders in line 48.
# The script requires folders with xlsx tables, folder names should start from "Imaris_".
# Enable lines 118-130 to plot 3D visualisation of tracks.
###############################################

#Libraries
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(stringr)
library(plotly)
library(zoo)
library(ggpubr)
library(rlist)

# Set colors
cols <- c("Moving" = "#118ab2", "Stationary" = "#ef476f")
#choose multiple colors for gradient
col_gradient <- colorRampPalette(c("#390099", "#9E0059", "#FF0054", "#FF5400", "#ee9b00"))

# Create required functions
create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
###############################################

condition_folders <- Sys.glob(file.path("Imaris*")) #check that only required folders start with Tau*

for (condition in condition_folders) {
  print(condition)
  setwd(paste("C:/Users/NITru/OneDrive/Documents/PhD_work/Projects/Marina/3D_transport_analysis/",condition,sep=""))
  #setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  filenames <- Sys.glob(file.path("*.xlsx"))
  total_full <- data.frame()
  df_imaris_collect_full <- data.frame()
  
  for (x in filenames){  
    #x <-"250220_m1w1s2_2_imaris_seg2.xlsx"
    name <- gsub(".xlsx","",x)
    print(name)
    toMatch_1 <- c("021019", "071019", "120919","231019","201019")
    pattern_1 <- paste(toMatch_1,collapse="|")
    toMatch_2 <- c("200529")
    pattern_2 <- paste(toMatch_2,collapse="|")
    toMatch_3 <- c("200226", "200303", "200304")
    pattern_3 <- paste(toMatch_3,collapse="|")
    toMatch_4 <- c("250220")
    pattern_4 <- paste(toMatch_4,collapse="|")
    
    if (grepl(pattern_1, name)) {
      frame_rate <- 600
    } else if (grepl(pattern_2, name)){
      frame_rate <- 850
    } else if (grepl(pattern_3, name)){
      frame_rate <- 890
    } else if (grepl(pattern_4, name)){
      frame_rate <- 900
    } else {
      print("Another date of acquisition!")
      break
    }
    
    df <- read_excel(x, sheet = "Position", skip = 1)
    df$TrackID <- as.factor(df$TrackID)
    df_subset <- df 
    df_subset$Time_s <- df_subset$Time*frame_rate/1000
    df_subset <- df_subset[c("TrackID", "Position X", "Position Y", "Position Z", "Time", "Time_s")]

    df_delta <- read_excel(x, sheet = "Displacement Delta Length", skip = 1)
    df_delta$Time_s <- df_delta$Time*frame_rate/1000
    df_delta <- df_delta[c("TrackID", "Displacement Delta Length", "Time", "Time_s")]
    total <- merge(df_subset,df_delta,by=c("TrackID", "Time","Time_s"))
    total <- total[order(total$Time_s),]
    total <- total[order(total$TrackID),]
    
    total$delta_time <- NA
    total$velocity <- NA
    df_list <- list()
    df_list <- split(total, total$TrackID)
    for(i in 1:length(df_list)){
      #print(df_list[[i]]$TrackID[1])
      if (nrow(df_list[[i]]) > 1) {
        for (j in 2:nrow(df_list[[i]])) {
          #IMPORTANT time is not always 1 apart, delta_x = delta_disp/delta_t
          df_list[[i]]$delta_time[[j]] <- (df_list[[i]]$Time_s[[j]]-df_list[[i]]$Time_s[[j-1]])
          # print((df_list[[i]]$`Displacement Delta Length`[[j]])/
          #         (df_list[[i]]$Time_s[[j]]-df_list[[i]]$Time_s[[j-1]]))
          df_list[[i]]$velocity[[j]] <- (df_list[[i]]$`Displacement Delta Length`[[j]])/
            (df_list[[i]]$Time_s[[j]]-df_list[[i]]$Time_s[[j-1]])
        }
      }  
    }
    total <- do.call(rbind, df_list)

    total$Movement <- ifelse(total$velocity <= 0.02, "Stationary", "Moving")
    total$Movement <- as.factor(total$Movement)
    total$Movement <- factor(total$Movement, levels=c("Moving", "Stationary"))
    total$Name <- name
    
    total_full <- rbind(total_full,total)
    
    # Stop number count to add to info
    df_stop_count <- na.omit(total[total$Movement=="Stationary",]) %>% group_by(TrackID) %>% dplyr::summarise(count(Movement))
    df_stop_count <- df_stop_count[c("TrackID", "freq")]
    colnames(df_stop_count) <- c("ID", "Stop_count")
    
    # Time spent not moving to add to info
    df_stop_time <- na.omit(total[total$Movement=="Stationary",]) %>% group_by(TrackID) %>% dplyr::summarise(sum(delta_time))
    colnames(df_stop_time) <- c("ID", "Time_stopped")
    
    # cbp <- col_gradient(nlevels(total$TrackID)) #gradient colors for tracks
    # cbp <- append(cbp, c("#adb5bd", "#000000")) #add colors for "Moving" and "Stationary"
    # p <- plot_ly(total, x=~`Position X`, y=~`Position Y`, z=~`Position Z`,  # can be added together but problematic coloring: marker = list(color = ~Movement, colors = cols))
    #              hoverinfo = 'text',
    #              text=~paste('</br> Time: ', Time, '[s]',
    #                          '</br> Delta time: ', sprintf("%0.2f", round(delta_time, digits = 2)), '[s]',
    #                          '</br> Delta movement: ', sprintf("%0.2f", round(velocity, digits = 2)), '[\U003BCm/s]'),
    #              type="scatter3d", mode="lines",
    #              color = ~TrackID, colors = cbp)  %>%
    #   add_markers(total, x=~`Position X`, y=~`Position Y`, z=~`Position Z`, color = ~Movement, colors = cols)
    # htmlwidgets::saveWidget(as_widget(p), paste("stops_",name,".html", sep=""))
    
    ##################
    # Extra information collect
    ##################
    
    df_duration <- read_excel(x, sheet = "Track Duration", skip = 1)
    df_duration <- df_duration[c("ID", "Track Duration")]
    df_displ <- read_excel(x, sheet = "Track Displacement Length", skip = 1)
    df_displ <- df_displ[c("ID", "Track Displacement Length")]
    df_straight <- read_excel(x, sheet = "Track Straightness", skip = 1)
    df_straight <- df_straight[c("ID", "Track Straightness")]
    df_length <- read_excel(x, sheet = "Track Length", skip = 1)
    df_length <- df_length[c("ID", "Track Length")]
    df_imaris_collect <- merge(df_duration,df_displ,by="ID")
    df_imaris_collect <- merge(df_imaris_collect,df_straight,by="ID")
    df_imaris_collect <- merge(df_imaris_collect,df_length,by="ID")

    df_imaris_collect$ID <- as.factor(df_imaris_collect$ID)
    df_imaris_collect$Velocity <- df_imaris_collect$`Track Displacement Length`/df_imaris_collect$`Track Duration`
    df_imaris_collect$Processivity <- df_imaris_collect$`Track Displacement Length`/df_imaris_collect$`Track Length`
    df_imaris_collect$Speed <- df_imaris_collect$`Track Length`/df_imaris_collect$`Track Duration`
    df_imaris_collect$Name <- name
    df_imaris_collect <- merge(df_imaris_collect,df_stop_count,by="ID") 
    df_imaris_collect <- merge(df_imaris_collect,df_stop_time,by="ID") 
    
    df_imaris_collect_full <- rbind(df_imaris_collect_full,df_imaris_collect)
  }
  write.table(total_full, paste("../positions_", condition, ".csv", sep=""), sep = ";",dec = '.', row.names = FALSE, col.names = TRUE)
  write.table(df_imaris_collect_full, paste("../imaris_collect_", condition, ".csv", sep=""), sep = ";",dec = '.', row.names = FALSE, col.names = TRUE)
}
  
  
