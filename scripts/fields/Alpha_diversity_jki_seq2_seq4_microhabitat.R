###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")
library("vegan")
library("ggpubr")
library("ggplot2")

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
alfa_div_comparison_microhabitat <- list(c("BS", "RH"), c("BS","RP"),c("RH","RP"))

#Plots SHANNON
#to compare microhabitats in Go
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "Microhabitat", "shannon", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 6.0, label.x = 1, size=5) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .2, size=3, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 16, height = 11, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "Microhabitat", "shannon", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 6.0, label.x = 1, size=5) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .2, size=3, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 16, height = 11, units = "cm",dpi = 300)


#to compare microhabitats in Ki
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "Microhabitat", "shannon", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 6.0, label.x = 1, size=5) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .2, size=3, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 11.5, height = 11, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "Microhabitat", "shannon", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(6.0, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 6.0, label.x = 1, size=5) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .2, size=3, label.y = c(7.7, 7.8, 7.9))
plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 11.5, height = 11, units = "cm",dpi = 300)



############################
#Plots CHAO1
#to compare microhabitats in Go
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_chao1_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "Microhabitat", "chao1", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_chao1_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_chao1_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)

plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_chao1_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "Microhabitat", "chao1", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_chao1_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_chao1_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)



#to compare microhabitats in Ki
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_chao1_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "Microhabitat", "chao1", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_chao1_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_chao1_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 12, units = "cm",dpi = 300)

plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_chao1_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "Microhabitat", "chao1", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(2500, 5000), breaks = seq(2500, 4500, by = 500), label = c("2500", "3000", "3500", "4000", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format",  bracket.size = .3, size=4, label.y = c(4500, 4650, 4800))
plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_chao1_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_chao1_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 12, units = "cm",dpi = 300)


####################################
#Plots PIELOU
#to compare microhabitats in Go
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_pielou_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Go_2020_filt_rarefied.meta, "Microhabitat", "pielou", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_pielou_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_pielou_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_pielou_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Go_2021_filt_rarefied.meta, "Microhabitat", "pielou", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_pielou_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_pielou_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)

#to compare microhabitats in Ki
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_pielou_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_filt_rarefied.meta, "Microhabitat", "pielou", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_pielou_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_pielou_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_pielou_microhabitat <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_filt_rarefied.meta, "Microhabitat", "pielou", fill = "Microhabitat", color = "Microhabitat", palette = c("#b05644", "#d9b967", "#57896a"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Rotation, scales = "free") +
  scale_y_continuous(limits = c(0.7, 0.95), breaks = seq(0.75, 0.95, by = 0.1), label = c("0.75", "0.85", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_microhabitat, label = "p.format", bracket.size = .3, size=4, label.y = c(0.92, 0.93, 0.94))
plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_pielou_microhabitat

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_pielou_microhabitat.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 12, units = "cm",dpi = 300)


###########################################
#arrange multiple ggplots in one single figure ALL
alpha_microh_plots_all <- ggarrange(plot_psO_jki_seq2_seq4_Go_2020_filt_rarefied_shannon_microhabitat,
                                    plot_psO_jki_seq2_seq4_Ki_2020_filt_rarefied_shannon_microhabitat,
                                    plot_psO_jki_seq2_seq4_Go_2021_filt_rarefied_shannon_microhabitat,
                                    plot_psO_jki_seq2_seq4_Ki_2021_filt_rarefied_shannon_microhabitat,
                                    ncol =2, nrow =2, legend="right", widths = c(1.5, 1), heights = c(1, 1), common.legend = TRUE, label.y = 1)
alpha_microh_plots_all

ggsave("alpha_microh_plots_all.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 40, height = 20, units = "cm",dpi = 300, device = "svg")
###########################################

# The end!  Tchüss!!

