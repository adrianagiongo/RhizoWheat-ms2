###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")
library("vegan")
library("ggpubr")
library("ggplot2")

#W1= "#e7e5cc", W2= "#9cc184", W3= "#447243", WM= "#1e3d14"

###Select groups and subset data from the cleaned dataset

##Creating rarefaction curve of non-rarefied samples (files separated by year were created in the Ordination script)
#Using rarecurve()
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2020_filt)), step=20, cex=0.5, col = "blue")
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2021_filt)), step=20, cex=0.5, col = "blue")
# 
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2020_filt)), step=20, cex=0.5, col = "green")
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2021_filt)), step=20, cex=0.5, col = "green")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq2_seq4_Go_2020_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Go_2020_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Go_2020_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Go_2020_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Go_2020_filt_rarefied)

psO_jki_seq2_seq4_Go_2021_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Go_2021_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Go_2021_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Go_2021_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Go_2021_filt_rarefied)

psO_jki_seq2_seq4_Ki_2020_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Ki_2020_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Ki_2020_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Ki_2020_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Ki_2020_filt_rarefied)

psO_jki_seq2_seq4_Ki_2021_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Ki_2021_filt, rngseed=2022, sample.size = 22000, trimOTUs=TRUE) #sample A291KIT2W1R3RH was removed. 
psO_jki_seq2_seq4_Ki_2021_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Ki_2021_filt_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2020_filt_rarefied)), step=20, cex=0.5, col = "black")
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2021_filt_rarefied)), step=20, cex=0.5, col = "black")
# 
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2020_filt_rarefied)), step=20, cex=0.5, col = "black")
# rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2021_filt_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Go_2020_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Go_2021_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Ki_2020_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Ki_2021_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Go_2020_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Go_2021_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Ki_2020_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Ki_2021_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Go_2020_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Go_2021_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Ki_2020_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Ki_2021_filt_rarefied, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity
#Prepare file using meta function
psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Go_2020_filt_rarefied)
psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Go_2021_filt_rarefied)

psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Ki_2020_filt_rarefied)
psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Ki_2021_filt_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for Go
psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_1$observed
psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Go_2020_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Go_2020_filt_rarefied_meta.csv")

psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_1$observed
psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Go_2021_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Go_2021_filt_rarefied_meta.csv")

#for Ki
psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_1$observed
psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Ki_2020_filt_rarefied_meta.csv")

psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_1$observed
psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Ki_2021_filt_rarefied_meta.csv")

#Select variable of comparison
alfa_div_comparison_rotation_Go <- list(c("W1","W2"),c("W2","WM"),c("W1","WM"))
alfa_div_comparison_rotation_Ki <- list(c("W1","W3"))

#Plots SHANNON
#to compare rotations in Go
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size=22), strip.text.y = element_text(size=22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=22, colour = "black"))
  #stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, label = "p.format",  bracket.size = .3, size=4, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 15, height = 8, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size=22), strip.text.y = element_text(size=22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=22, colour = "black"))
  #stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, label = "p.format",  bracket.size = .3, size=4, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 15, height = 8, units = "cm",dpi = 300)




plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size=22), strip.text.y = element_text(size=22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=22, colour = "black"))
#stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, label = "p.format",  bracket.size = .3, size=4, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 8, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size=22), strip.text.y = element_text(size=22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size=22, colour = "black"))
  #stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, label = "p.format",  bracket.size = .3, size=4, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_rotation.eps", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 8, units = "cm",dpi = 300, device = "eps")



###########################################
#arrange multiple ggplots in one single figure ALL
alpha_rot_plots_all <- ggarrange(plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_rotation,
                                plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_rotation,
                                plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_rotation,
                                plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_rotation,
                                ncol =2, nrow =2, legend="right", widths = c(1, 1), common.legend = TRUE, label.y = 1)
alpha_rot_plots_all

ggsave("alpha_rot_plots_all.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 34, height = 20, units = "cm",dpi = 300, device = "svg")
###########################################




## I stopped here (19.02)

############################
#Plots CHAO1
#to compare rotations in Go
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_chao1_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)

plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_chao1_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)



#to compare rotations in Ki
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_chao1_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)

plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_chao1_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


####################################
#Plots PIELOU
#to compare rotations in Go
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_pielou_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_pielou_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)

#to compare rotations in Ki
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_pielou_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_pielou_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)



# The end!  Tchüss!!

