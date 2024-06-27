#Load packages
library("ggplot2")

#upload your data to R - exchange "Your_csv_file.csv" with the name of your csv file
Climate = read.csv("~/Documents/R_analysis/jki_seq2_seq4/data_jki_seq2_seq4/Climate1.csv", header = TRUE)
Climate


# Assuming you have a data frame named 'Climate' with columns 'Soil', 'Precipitation', and 'Season'

# Convert Soil to a factor with specific order
Climate$Soil <- factor(Climate$Soil, levels = c("Silty", "Sandy"))

# Create the barplot
Climate_barplot <- ggplot(data = Climate, aes(x = Soil, y = Precipitation, fill = Season)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  geom_bar(stat = "identity", color="black", position = position_dodge(0.9)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 55), breaks = seq(0, 55, by = 10), label = c("0", "10", "20", "30", "40", "50")) +
  ylab("Precipitation (mm)") +
  theme(axis.text.x = element_text(angle = 0, size=22, color="black", vjust = 0.5, hjust = 0.5)) +
  theme(axis.text.y = element_text(size=22, color="black")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=22, color="black")) +
  theme(strip.text.x = element_text(size = 22, face="bold")) +
  scale_fill_manual(values = c("#bc9c90", "#95dcf8"))
Climate_barplot

ggsave("Climate_barplot.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/", width = 10, height = 9, units = "cm",dpi = 300, device = "svg")

