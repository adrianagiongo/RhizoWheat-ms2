## Bubble plot
#### Adriana Giongo

#Load packages
library("ggplot2")
library("reshape2")

#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
isolates_jki_seq2_seq4 = read.csv("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/jki_seq2_seq4_isolates_bubble.csv", header = TRUE)
isolates_jki_seq2_seq4

#convert data frame from a "wide" format to a "long" format
isolates_jki_seq2_seq4_melt <- melt(isolates_jki_seq2_seq4, id = c("Sample", "Microhabitat", "Soil"))
isolates_jki_seq2_seq4_melt

colours = c("#b05644", "#d9b967", "#57896a")

isolates_jki_seq2_seq4_melt$Sample <- factor(isolates_jki_seq2_seq4_melt$Sample,levels=unique(isolates_jki_seq2_seq4$Sample))

bubbleplot_isolates_jki_seq2_seq4_melt = ggplot(isolates_jki_seq2_seq4_melt, aes(x = Sample, y = variable)) + 
  theme_bw() +
  facet_grid(~Soil, scales = "free_x", space = "free_x") +
  geom_point(aes(size = value, fill = Microhabitat), alpha = 0.75, shape = 21) + 
  guides(fill = guide_legend(override.aes = list(size=5))) +
  scale_size_continuous(limits = c(0.1, 15), range = c(0.1,15), breaks = c(1, 3, 6, 9, 12)) + 
  labs(x= "", y = "", size = "Number of isolates", fill = "Microhabitat")  + 
  theme(legend.key = unit(1, units = "cm"))+
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 12, angle = 0, vjust = 0.3, hjust = 1),
        axis.text.y = element_text(colour = "black", size = 12, face = "italic"),
        legend.text = element_text(size = 12, colour ="black"),
        legend.title = element_text(size = 12),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +
  scale_fill_manual(values = colours)
bubbleplot_isolates_jki_seq2_seq4_melt

ggsave("bubbleplot_isolates_jki_seq2_seq4_melt.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 22, height = 12, units = "cm", dpi = 300, device = "png")



### The end! Have fun : )