#Load packages
library("ggplot2")

# Shoot #80c698
# Roots #c6ba9c
# Disease reduction #549ccb

#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
Isolates_GH07_B = read.csv("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/Isolates_GH07_B.csv", header = TRUE)
Isolates_GH07_B


# Define custom order for treatments
custom_order <- c('B_mycoides_57_24', 'B_pumilus_45_39', 
                  'P_extremorientalis_48_38', 'P_proteolytica_68_48',
                  'S_panacis_40_16')

# Convert Treatment to factor with custom order
Isolates_GH07_B$Treatment <- factor(Isolates_GH07_B$Treatment, levels = custom_order)


# Create the barplot
Isolates_GH07_B_barplot <- ggplot(data = Isolates_GH07_B, aes(x = Treatment, y = Result, fill = Test)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_bar(stat = "identity", color="black", position = position_dodge(0.9)) +
  theme_bw() +
  #scale_y_continuous(limits = c(-10, 60), breaks = seq(-10, 60, by = 10,  label = c("-10", "10", "30", "20", "40", "60"))) +
  ylab("(%)") +
  scale_x_discrete("Treatment", labels = c('B_mycoides_57_24' = "B. mycoides 57.24",
                                           'B_pumilus_45_39' = "B. pumilus 45.39",
                                           'P_brassicacearum_37_15' = "P. brassicacearum 37.15",
                                           'P_extremorientalis_48_38' = "P. extremorientalis 48.38",
                                           'P_proteolytica_68_48' = "P. proteolytica 68.48",
                                           'S_panacis_40_16' = "S. panacis 40.16")) + 
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5, hjust = 1)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=18, color="black")) +
  theme(strip.text.x = element_text(size = 18, face="bold")) +
  scale_fill_manual(values = c("#c6ba9c", "#80c698"))
Isolates_GH07_B_barplot

ggsave("Isolates_GH07_B_barplot.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 20, height = 20, units = "cm",dpi = 300)







#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
Isolates_GH10_C = read.csv("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/Isolates_GH10_C.csv", header = TRUE)
Isolates_GH10_C


# Define custom order for treatments
custom_order <- c('B_mycoides_57_24', 'B_pumilus_45_39', 'E_morelensis_45_46', 
                  'P_brassicacearum_37_15', 'P_extremorientalis_48_38', 'P_lini_49_51', 
                  'P_proteolytica_68_48', 'P_kribbensis_57_20', 'S_panacis_40_16', 'Consortium_1')

# Convert Treatment to factor with custom order
Isolates_GH10_C$Treatment <- factor(Isolates_GH10_C$Treatment, levels = custom_order)


# Create the barplot
Isolates_GH10_C_barplot <- ggplot(data = Isolates_GH10_C, aes(x = Treatment, y = Result, fill = Test)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_bar(stat = "identity", color="black", position = position_dodge(0.9)) +
  theme_bw() +
  #scale_y_continuous(limits = c(-10, 60), breaks = seq(-10, 60, by = 10,  label = c("-10", "10", "30", "20", "40", "60"))) +
  ylab("(%)") +
  scale_x_discrete("Treatment", labels = c('B_mycoides_57_24' = "B. mycoides 57.24",
                                           'B_pumilus_45_39' = "B. pumilus 45.39",
                                           'E_morelensis_45_46' = "E. morelensis 45.46",
                                           'P_brassicacearum_37_15' = "P. brassicacearum 37.15",
                                           'P_extremorientalis_48_38' = "P. extremorientalis 48.38",
                                           'P_lini_49_51' = "P. lini 49.51",
                                           'P_proteolytica_68_48' = "P. proteolytica 68.48",
                                           'P_kribbensis_57_20' = "P. kribbensis 57.20",
                                           'S_panacis_40_16' = "S. panacis 40.16",
                                           'Consortium_1' = "Consortium 1")) + 
  theme(axis.text.x = element_text(angle = 90, size=18, color="black", vjust = 0.5, hjust = 1)) +
  theme(axis.text.y = element_text(size=18, color="black")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=18, color="black")) +
  theme(strip.text.x = element_text(size = 18, face="bold")) +
  scale_fill_manual(values = c("#549ccb", "#c6ba9c"))
Isolates_GH10_C_barplot

ggsave("Isolates_GH10_C_barplot.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 20, height = 20, units = "cm",dpi = 300)




