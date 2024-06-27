## Core - bubble plot - jki_seq2_seq4
#### Adriana Giongo
#### (21.05.2023) 

#Load packages
library("ggplot2")
library("reshape2")


####### Example --> CORE
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
core_BS = read.csv("~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Bubble_plot_jki_seq2_seq4/jki_seq2_seq4_core_BS.csv", header = TRUE)
core_BS

#convert data frame from a "wide" format to a "long" format
core_BS_melt <- melt(core_BS, id = c("Sample", "Field"))
core_BS_melt

#colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
#             "#F6AE2D","#86BBD8")
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57")

core_BS_melt$Sample <- factor(core_BS_melt$Sample,levels=unique(core_BS_melt$Sample))

plot_core_BS_melt = ggplot(core_BS_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Field), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 50), range = c(1,20), breaks = c(1,5,10,15,20)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Field")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
  #scale_y_discrete(limits = rev(levels(core_BS_melt$variable))) 

plot_core_BS_melt

ggsave("plot_core_BS_melt.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Core_bubble_plot_jki_seq2_seq4/", width = 30, height = 80, units = "cm", dpi = 300, device = "tiff")


### CORE_jki_seq2_seq4_Go_2020

#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
core_Go_2020_BS_RH_RP = read.csv("~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Bubble_plot_jki_seq2_seq4/core_Go_2020_BS_RH_RP.csv", header = TRUE)
core_Go_2020_BS_RH_RP

#convert data frame from a "wide" format to a "long" format
core_Go_2020_BS_RH_RP_melt <- melt(core_Go_2020_BS_RH_RP, id = c("Sample", "Field"))
core_Go_2020_BS_RH_RP_melt

colours = c( "#e7e5cc",  "#9cc184", "#1e3d14")

core_Go_2020_BS_RH_RP_melt$Sample <- factor(core_Go_2020_BS_RH_RP_melt$Sample,levels=unique(core_Go_2020_BS_RH_RP$Sample))

bubble_plot_core_Go_2020_BS_RH_RP_melt = ggplot(core_Go_2020_BS_RH_RP_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Field), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 50), range = c(1,50), breaks = c(0.5,1,2.5)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Field")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 13), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.title = element_text(size = 13, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
bubble_plot_core_Go_2020_BS_RH_RP_melt

ggsave("bubble_plot_core_Go_2020_BS_RH_RP_melt.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Bubble_plot_jki_seq2_seq4/", width = 28, height = 22, units = "cm", dpi = 200, device = "tiff")


### CORE_jki_seq2_seq4_Ki_2020

#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
core_Ki_2020_BS_RH_RP = read.csv("~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Bubble_plot_jki_seq2_seq4/core_Ki_2020_BS_RH_RP.csv", header = TRUE)
core_Ki_2020_BS_RH_RP

#convert data frame from a "wide" format to a "long" format
core_Ki_2020_BS_RH_RP_melt <- melt(core_Ki_2020_BS_RH_RP, id = c("Sample", "Field"))
core_Ki_2020_BS_RH_RP_melt

colours = c( "#e7e5cc",  "#447243")

core_Ki_2020_BS_RH_RP_melt$Sample <- factor(core_Ki_2020_BS_RH_RP_melt$Sample,levels=unique(core_Ki_2020_BS_RH_RP$Sample))

bubble_plot_core_Ki_2020_BS_RH_RP_melt = ggplot(core_Ki_2020_BS_RH_RP_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Field), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 50), range = c(1,50), breaks = c(0.5,1,2.5)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Field")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 13), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.title = element_text(size = 13, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
bubble_plot_core_Ki_2020_BS_RH_RP_melt

ggsave("bubble_plot_core_Ki_2020_BS_RH_RP_melt.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Bubble_plot_jki_seq2_seq4/", width = 26, height = 10, units = "cm", dpi = 200, device = "tiff")




