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
#(2)Year = 2020 or 2021 
#(3)Plant time = T1 or T2 or T3 
#(4)Rotation = W1, W2, WM (Go) or W1, W3 (Ki) - / (MetBrewer:"VanGogh3") W1= "#e7e5cc", W2= "#9cc184", W3= "#447243", WM= "#1e3d14"
#(5)Microhabitats = BS, RH, RP


psO_jki_seq2_seq4_Go <- subset_samples(psO_jki_seq2_seq4, Location == "Go")
psO_jki_seq2_seq4_Go

psO_jki_seq2_seq4_Ki <- subset_samples(psO_jki_seq2_seq4, Location == "Ki")
psO_jki_seq2_seq4_Ki


###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

psO_jki_seq2_seq4_Go_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go) > 0, psO_jki_seq2_seq4_Go)
psO_jki_seq2_seq4_Go_filt

psO_jki_seq2_seq4_Ki_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki) > 0, psO_jki_seq2_seq4_Ki)
psO_jki_seq2_seq4_Ki_filt



##Creating rarefaction curve of non-rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_filt)), step=20, cex=0.5, col = "blue")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_filt)), step=20, cex=0.5, col = "green")

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq2_seq4_Go_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Go_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Go_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Go_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Go_filt_rarefied)

psO_jki_seq2_seq4_Ki_filt_rarefied<-rarefy_even_depth(psO_jki_seq2_seq4_Ki_filt, rngseed=2022,sample.size = min(sample_sums(psO_jki_seq2_seq4_Ki_filt)),trimOTUs=TRUE)
psO_jki_seq2_seq4_Ki_filt_rarefied
sample_sums(psO_jki_seq2_seq4_Ki_filt_rarefied)

##Creating rarefaction curve of rarefied samples
#Using rarecurve()
rarecurve(t(otu_table(psO_jki_seq2_seq4_Go_filt_rarefied)), step=20, cex=0.5, col = "black")
rarecurve(t(otu_table(psO_jki_seq2_seq4_Ki_filt_rarefied)), step=20, cex=0.5, col = "black")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Go_filt_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_1 <- richness(psO_jki_seq2_seq4_Ki_filt_rarefied, c("observed","chao1"), detection = 0)

alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Go_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_2 <- evenness(psO_jki_seq2_seq4_Ki_filt_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)

alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Go_filt_rarefied, index = 'shannon', zeroes = TRUE)
alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_Ki_filt_rarefied, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity
#Prepare file using meta function
psO_jki_seq2_seq4_Go_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Go_filt_rarefied)
psO_jki_seq2_seq4_Ki_filt_rarefied.meta <- meta(psO_jki_seq2_seq4_Ki_filt_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for Go
psO_jki_seq2_seq4_Go_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_1$observed
psO_jki_seq2_seq4_Go_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Go_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Go_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Go_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Go_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Go_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Go_filt_rarefied_meta.csv")

#for Ki
psO_jki_seq2_seq4_Ki_filt_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_1$observed
psO_jki_seq2_seq4_Ki_filt_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_1$chao1
psO_jki_seq2_seq4_Ki_filt_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_2$pielou
psO_jki_seq2_seq4_Ki_filt_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_Ki_filt_rarefied_3$shannon

head(psO_jki_seq2_seq4_Ki_filt_rarefied.meta)
write.csv(psO_jki_seq2_seq4_Ki_filt_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_Ki_filt_rarefied_meta.csv")


#Select variable of comparison
alfa_div_comparison_rotation <- list(c("W1","W2"),c("W2","WM"),c("W1","WM"))

#Plots SHANNON
#to compare rotations in Go
plot_psO_jki_seq2_seq4_Go_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(Microhabitat~Stage, scales = "free") +
  scale_y_continuous(limits = c(4.0, 8.0), breaks = seq(4, 8, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 7.7, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(7.0, 7.2, 7.4))
plot_psO_jki_seq2_seq4_Go_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_filt_rarefied_shannon_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 28, units = "cm",dpi = 300)


#to compare rotations in Ki 
plot_psO_jki_seq2_seq4_Ki_filt_rarefied_shannon_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_filt_rarefied.meta, "Rotation", "shannon", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(Microhabitat~Stage, scales = "free") +
  scale_y_continuous(limits = c(4.0, 8.0), breaks = seq(4, 8, by = 1.0), label = c("4.0", "5.0", "6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 7.7, label.x = 1, size=5)
plot_psO_jki_seq2_seq4_Ki_filt_rarefied_shannon_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_filt_rarefied_shannon_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 14, height = 28, units = "cm",dpi = 300)



############################
#Plots CHAO1
#to compare rotations in Go
plot_psO_jki_seq2_seq4_Go_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(Microhabitat~Stage, scales = "free") +
  scale_y_continuous(limits = c(600, 2100), breaks = seq(600, 2100, by = 300), label = c("600", "900", "1200", "1500", "1800", "2100")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "kruskal", label= "p", label.y = 1950, label.x = 1, size=5) +
  stat_compare_means(comparisons = alfa_div_comparison_rotation, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(1850, 1950, 2050))

plot_psO_jki_seq2_seq4_Go_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_filt_rarefied_chao1_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 18, height = 28, units = "cm",dpi = 300)

#to compare rotations in Ki
plot_psO_jki_seq2_seq4_Ki_filt_rarefied_chao1_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_filt_rarefied.meta, "Rotation", "chao1", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(Microhabitat~Stage, scales = "free") +
  scale_y_continuous(limits = c(600, 2100), breaks = seq(600, 2100, by = 300), label = c("600", "900", "1200", "1500", "1800", "2100")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black")) +
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 1950, label.x = 1, size=5)
  #stat_compare_means(comparisons = alfa_div_comparison_rotation, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(1850, 1950, 2050))

plot_psO_jki_seq2_seq4_Ki_filt_rarefied_chao1_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_filt_rarefied_chao1_rotation.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 14, height = 28, units = "cm",dpi = 300)

####################################
#Plots PIELOU
#to compare rotations in Go
plot_psO_jki_seq2_seq4_Go_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Go_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#9cc184", "#1e3d14"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(Microhabitat~Stage, scales = "free") +
  scale_y_continuous(limits = c(0.60, 1.0), breaks = seq(0.60, 1.0, by = 0.1), label = c("0.60", "0.70", "0.80", "0.90", "1.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black"))
plot_psO_jki_seq2_seq4_Go_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Go_filt_rarefied_pielou_rotation.png", path = "~/Documents/R_analysis_2/jki_seq2_seq4/output/Alpha_div/", width = 18, height = 28, units = "cm",dpi = 300)

#to compare rotations in Ki
plot_psO_jki_seq2_seq4_Ki_filt_rarefied_pielou_rotation <- ggboxplot(psO_jki_seq2_seq4_Ki_filt_rarefied.meta, "Rotation", "pielou", fill = "Rotation", color = "Rotation", palette = c("#e7e5cc", "#447243"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  facet_grid(Microhabitat~Stage, scales = "free") +
  scale_y_continuous(limits = c(0.60, 1.0), breaks = seq(0.60, 1.0, by = 0.1), label = c("0.60", "0.70", "0.80", "0.90", "1.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 30, colour = "black")) +
  theme(axis.text.x = element_text(size=30, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=30)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30), legend.text = element_text(size = 30, colour = "black"))
plot_psO_jki_seq2_seq4_Ki_filt_rarefied_pielou_rotation

ggsave("plot_psO_jki_seq2_seq4_Ki_filt_rarefied_pielou_rotation.png", path = "~/Documents/R_analysis_2/jki_seq2_seq4/output/Alpha_div/", width = 14, height = 28, units = "cm",dpi = 300)

## I stopped here : )

