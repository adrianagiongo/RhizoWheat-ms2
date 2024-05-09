###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")
library("vegan")
library("ggpubr")

###Select groups and subset data from the cleaned dataset

###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors
#(1)Location = Go or Ki
#(2)Year = T2 or T3 
#(3)Plant time = T1 or T2 or T3 
#(4)Rotation = W1, W2, WM (Go) or W1, W3 (Ki) 
#(5)Microhabitats = BS, RH, RP

#Go Stages
psO_jki_seq2_seq4_Go_2020_T2 <- subset_samples(psO_jki_seq2_seq4, Location == "Go" & Year == "Y2020" & Stage == "T2")
psO_jki_seq2_seq4_Go_2020_T2

psO_jki_seq2_seq4_Go_2021_T2 <- subset_samples(psO_jki_seq2_seq4, Location == "Go" & Year == "Y2021" & Stage == "T2")
psO_jki_seq2_seq4_Go_2021_T2

#Ki Stages
psO_jki_seq2_seq4_Ki_2020_T2 <- subset_samples(psO_jki_seq2_seq4, Location == "Ki" & Year == "Y2020" & Stage == "T2")
psO_jki_seq2_seq4_Ki_2020_T2

psO_jki_seq2_seq4_Ki_2021_T2 <- subset_samples(psO_jki_seq2_seq4, Location == "Ki" & Year == "Y2021" & Stage == "T2")
psO_jki_seq2_seq4_Ki_2021_T2

###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

#Go Year Stage T2
psO_jki_seq2_seq4_Go_2020_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_T2) > 0, psO_jki_seq2_seq4_Go_2020_T2)
psO_jki_seq2_seq4_Go_2020_T2_filt

psO_jki_seq2_seq4_Go_2021_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_T2) > 0, psO_jki_seq2_seq4_Go_2021_T2)
psO_jki_seq2_seq4_Go_2021_T2_filt

#Ki Year Stage T2
psO_jki_seq2_seq4_Ki_2020_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_T2) > 0, psO_jki_seq2_seq4_Ki_2020_T2)
psO_jki_seq2_seq4_Ki_2020_T2_filt

psO_jki_seq2_seq4_Ki_2021_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_T2) > 0, psO_jki_seq2_seq4_Ki_2021_T2)
psO_jki_seq2_seq4_Ki_2021_T2_filt


##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2020_T2_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2021_T2_filt)), step=20, cex=0.5, col = "red")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2020_T2_filt)), step=20, cex=0.5, col = "red")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2021_T2_filt)), step=20, cex=0.5, col = "green")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Go_2020_T2_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Go_2020_T2_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied)

psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Go_2021_T2_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Go_2021_T2_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied)

psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Ki_2020_T2_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Ki_2020_T2_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied)

psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Ki_2021_T2_filt, rngseed=2022,sample.size = 22000,trimOTUs=TRUE)  #Sample A291KIT2W1R3RH removed
psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied, index = 'shannon', zeroes = TRUE)

alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity
#Prepare file using meta function
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied)
psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied)

psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied)
psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for Go_T2
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_1$observed
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Go_2020_T2_T3_filt_rarefied_meta.csv")


psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_1$observed
psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Go_2021_T2_T3_filt_rarefied_meta.csv")


#for Ki_T2
psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_1$observed
psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Ki_2020_T2_T3_filt_rarefied_meta.csv")


psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_1$observed
psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Ki_2021_T2_T3_filt_rarefied_meta.csv")


#Select variable of comparison
alfa_div_comparison_rotation_Go <- list(c("W1","W2"),c("W2","WM"),c("W1","WM"))
alfa_div_comparison_rotation_Ki <- list(c("W1","W3"))


#Plots SHANNON
#to compare rotations in Go_T2
plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(5.5, 7.0), breaks = seq(5.5, 7.0, by = 0.5), label = c("5.5", "6.0", "6.5", "7.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 6.95, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(6.6, 6.7, 6.8))
plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_shannon_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)



plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(5.5, 7.0), breaks = seq(5.5, 7.0, by = 0.5), label = c("5.5", "6.0", "6.5", "7.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 6.95, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(6.6, 6.7, 6.8))
plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_shannon_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)




#to compare rotations in Ki_T2 
plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(5.5, 7.0), breaks = seq(5.5, 7.0, by = 0.5), label = c("5.5", "6.0", "6.5", "7.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(6.6, 6.7, 6.8))
plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_shannon_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)



plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(5.5, 7.0), breaks = seq(5.5, 7.0, by = 0.5), label = c("5.5", "6.0", "6.5", "7.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(6.6, 6.7, 6.8))
plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_shannon_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


############################
#Plots CHAO1
#to compare rotations in Go_T2
plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(800, 2000), breaks = seq(800, 2000, by = 300), label = c("800", "1100", "1400", "1700", "2000")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 2000, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(1800, 1850, 1900))
plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_chao1_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(800, 2000), breaks = seq(800, 2000, by = 300), label = c("800", "1100", "1400", "1700", "2000")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 2000, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(1800, 1850, 1900))
plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_chao1_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


#to compare rotations in Ki_T2
plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(800, 2000), breaks = seq(800, 2000, by = 300), label = c("800", "1100", "1400", "1700", "2000")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(1800, 1850, 1900))
plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_chao1_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(800, 2000), breaks = seq(800, 2000, by = 300), label = c("800", "1100", "1400", "1700", "2000")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(1800, 1850, 1900))
plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_chao1_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


####################################
#Plots PIELOU
#to compare rotations in Go_T2
plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.70, 1.00), breaks = seq(0.70, 1.00, by = 0.1), label = c("0.70", "0.80", "0.90", "1.00")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(0.91, 0.92, 0.93))
plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2020_T2_filt_rarefied_pielou_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.70, 1.00), breaks = seq(0.70, 1.00, by = 0.1), label = c("0.70", "0.80", "0.90", "1.00")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 0.98, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Go, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(0.91, 0.92, 0.93))
plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_2021_T2_filt_rarefied_pielou_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


#to compare rotations in Ki_T2
plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.70, 1.00), breaks = seq(0.70, 1.00, by = 0.1), label = c("0.70", "0.80", "0.90", "1.00")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(0.91, 0.92, 0.93))
plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2020_T2_filt_rarefied_pielou_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#cfceb7", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(~Microhabitat, scales = "free") +
  scale_y_continuous(limits = c(0.70, 1.00), breaks = seq(0.70, 1.00, by = 0.1), label = c("0.70", "0.80", "0.90", "1.00")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation_Ki, method = "wilcox.test", label = "p", bracket.size = .3, size=3, label.y = c(0.91, 0.92, 0.93))
plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_2021_T2_filt_rarefied_pielou_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 12, units = "cm",dpi = 300)


## I stopped here : )



###########################################
#arrange multiple ggplots in one single figure ALL
#alpha_plots_all_pielou_rotation <- ggarrange(plot_psO_jki_seq2_seq4_Go_T2_filt_rarefied_pielou_rotation, plot_psO_jki_seq2_seq4_Ki_T2_filt_rarefied_pielou_rotation, 
#                                             plot_psO_jki_seq2_seq4_Go_T3_filt_rarefied_pielou_rotation, plot_psO_jki_seq2_seq4_Ki_T3_filt_rarefied_pielou_rotation, 
#                                             ncol =2, nrow =2, legend="right", widths = c(1, 1), common.legend = TRUE, label.y = 1)
#alpha_plots_all_pielou_rotation

#ggsave("alpha_plots_all_pielou_rotation.png", path = "~/Documents/R_analysis_2/jki_seq2_seq4/output/Alpha_div/", width = 60, height = 50, units = "cm",dpi = 300)
###########################################

