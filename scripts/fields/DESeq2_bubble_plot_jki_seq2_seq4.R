## Core - bubble plot - jki_seq5

#Load packages
library("ggplot2")
library("reshape2")

### T1
####### 
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
bubble_sigtab_T1 = read.csv("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/input_boxplot_deseq2_T1.csv", header = TRUE)
bubble_sigtab_T1

#convert data frame from a "wide" format to a "long" format
bubble_sigtab_T1_melt <- melt(bubble_sigtab_T1, id = c("Sample", "Rotation"))
bubble_sigtab_T1_melt

colours = c("#d7d3a9", "#74a553")

bubble_sigtab_T1_melt$Sample <- factor(bubble_sigtab_T1_melt$Sample,levels=unique(bubble_sigtab_T1_melt$Sample))

plot_bubble_sigtab_T1_melt = ggplot(bubble_sigtab_T1_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Rotation, stroke = 0.1), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 5), range = c(0,10), breaks = c(0,0.1,0.5,1,5)) + 
  labs(x = "", y = "", size = "(%)", fill = "Rotation")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12,  angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
plot_bubble_sigtab_T1_melt

ggsave("plot_bubble_sigtab_T1_melt.tiff", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Bubble_plot_jki_seq5/", width = 18, height = 10, units = "cm", dpi = 300, device = "tiff")


### T2
####### 
#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
bubble_sigtab_T2 = read.csv("~/Documents/R_analysis/jki_seq5/output_jki_seq5/Tables_jki_seq5/input_boxplot_deseq2_T2.csv", header = TRUE)
bubble_sigtab_T2

#convert data frame from a "wide" format to a "long" format
bubble_sigtab_T2_melt <- melt(bubble_sigtab_T2, id = c("Sample", "Rotation"))
bubble_sigtab_T2_melt

colours = c("#d7d3a9", "#74a553")

bubble_sigtab_T2_melt$Sample <- factor(bubble_sigtab_T2_melt$Sample,levels=unique(bubble_sigtab_T2_melt$Sample))

plot_bubble_sigtab_T2_melt = ggplot(bubble_sigtab_T2_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  geom_point(aes(size = value, fill = Rotation, stroke = 0.1), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 3), range = c(0,7), breaks = c(0,0.1,0.5,1,3)) + 
  labs(x = "", y = "", size = "(%)", fill = "Rotation")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12,  angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black",  size = 11), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours, guide = guide_legend(override.aes = list(size=5)))
plot_bubble_sigtab_T2_melt

ggsave("plot_bubble_sigtab_T2_melt.tiff", path = "~/Documents/R_analysis/jki_seq5/output_jki_seq5/Bubble_plot_jki_seq5/", width = 12, height = 10, units = "cm", dpi = 300, device = "tiff")

