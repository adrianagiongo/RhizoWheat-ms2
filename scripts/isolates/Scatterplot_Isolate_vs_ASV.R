#load libraries
library("ggplot2")
library("tidyr")
library("readxl")
library("dplyr")
library("RColorBrewer")

## ----> RH and RP only
#load data 
Isolate_vs_ASV<-read_excel("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/Isolate_vs_ASV.xlsx")
Isolate_vs_ASV


# Define color palette
color.site <- c("Sandy" = "#A34828",
               "Silty" = "#1d5b65",
               "Silty_Sandy" = "grey")

# Reorder the Taxa variable according to custom order and reverse the order
custom_order <- c('Arthrobacter', 'Bacillus', 'Chryseobacterium', 'Ensifer',
                  'Flavobacterium', 'Microbacterium', 'Paenibacillus',
                  'Pedobacter', 'Plantibacter', 'Pseudomonas',
                  'Rhodococcus', 'Sphingobacterium', 'Sphingomonas',
                  'Stenotrophomonas', 'Variovorax')

Isolate_vs_ASV$Taxa <- factor(Isolate_vs_ASV$Taxa, levels = rev(custom_order))

# Create the barplot
Isolate_vs_ASV_scatterplot <- ggplot(data = Isolate_vs_ASV, aes(x = Similarity, y = Taxa, fill = Sample, color = Sample)) +
  geom_point(position = position_jitter(width = 0.3, height = 0.3), size = 3, alpha = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16, color = "black", face = "italic"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  #scale_x_continuous(limits = c(94, 100), breaks = seq(93, 100, by = 3), label = c("94", "97", "100")) +
  scale_color_manual(values = color.site, breaks = c('Sandy', 'Silty', 'Silty_Sandy')) +
  xlab("Similarity (%)")+
  ylab("")

Isolate_vs_ASV_scatterplot

ggsave("Isolate_vs_ASV_scatterplot.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 20, height = 12, units = "cm",dpi = 300)


## ----> BS, RH and RP 
#load data 
Isolate_vs_ASV_BS_RH_RP<-read_excel("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/Isolate_vs_ASV_BS_RH_RP.xlsx")
Isolate_vs_ASV_BS_RH_RP


# Define color palette
color.site <- c("Sandy" = "#A34828",
                "Silty" = "#1d5b65",
                "Silty_Sandy" = "grey")

# Reorder the Taxa variable according to custom order and reverse the order
custom_order <- c('Arthrobacter', 'Bacillus', 'Chryseobacterium', 'Ensifer',
                  'Flavobacterium', 'Microbacterium', 'Paenibacillus',
                  'Pedobacter', 'Plantibacter', 'Pseudomonas',
                  'Rhodococcus', 'Sphingobacterium', 'Sphingomonas',
                  'Stenotrophomonas', 'Variovorax')

Isolate_vs_ASV_BS_RH_RP$Taxa <- factor(Isolate_vs_ASV_BS_RH_RP$Taxa, levels = rev(custom_order))

# Create the barplot
Isolate_vs_ASV_BS_RH_RP_scatterplot <- ggplot(data = Isolate_vs_ASV_BS_RH_RP, aes(x = Similarity, y = Taxa, fill = Sample, color = Sample)) +
  geom_point(position = position_jitter(width = 0.3, height = 0.3), size = 3, alpha = 0.8) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16, color = "black", face = "italic"),
        axis.title.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black")) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  #scale_x_continuous(limits = c(94, 100), breaks = seq(93, 100, by = 3), label = c("94", "97", "100")) +
  scale_color_manual(values = color.site, breaks = c('Sandy', 'Silty', 'Silty_Sandy')) +
  xlab("Similarity (%)")+
  ylab("")

Isolate_vs_ASV_BS_RH_RP_scatterplot

ggsave("Isolate_vs_ASV_BS_RH_RP_scatterplot.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 20, height = 12, units = "cm",dpi = 300)


