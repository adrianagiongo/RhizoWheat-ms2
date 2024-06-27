##BoxPlot comparison taxa
#Palette for Rotation  (W1 = "#e7e5cc" ("#A5A492" for borders), W2 = "#9cc184", WM = "#1e3d14", W3 =  "#447243")
#Palette for Microhabitat  RA = "#b05644", RH = "#d9b967", RP = "#57896A"
#Palette for Depth 0-15 = "#dcc8ba", 15-30 = "#DEB3AD", 30-60 = "#DE847B", 60-90 = "#B95C50", 90-120 = "#3B0404"
#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")
library("microbiome")
library("gridExtra")
library("ggpubr")
library("RColorBrewer")
library("ggbreak")

###### MOST ABUNDANT TAXA
#####################
###Select groups and subset data from the cleaned dataset

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq2_seq4, taxonomic.rank = "Annotation"))

##Using subset_samples()
## Samples separated by site, year, microhabitat in "Ordination_plot_jki_seq2_seq4_Site

### Subset taxa by microhabitat
#ps_pgpr_Go_2020_BS= subset_taxa(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
#ps_pgpr_Go_2020_RH= subset_taxa(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
ps_pgpr_Go_2020_RP= subset_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel, Annotation == "Cutibacterium" | Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")

#ps_pgpr_Ki_2020_BS= subset_taxa(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
#ps_pgpr_Ki_2020_RH= subset_taxa(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
ps_pgpr_Ki_2020_RP= subset_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")

#ps_pgpr_Go_2021_BS= subset_taxa(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
#ps_pgpr_Go_2021_RH= subset_taxa(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
ps_pgpr_Go_2021_RP= subset_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")

#ps_pgpr_Ki_2021_BS= subset_taxa(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
#ps_pgpr_Ki_2021_RH= subset_taxa(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")
ps_pgpr_Ki_2021_RP= subset_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel, Annotation == "Bradyrhizobium" | Annotation == "Devosia" | Annotation == "Nocardioides" | Annotation == "Sphingomonas")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)

###### Go
plot_abundance_stat_Go_high = function(phyloseq, title = "",
                                           Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#A5A492", "#9cc184", "#1e3d14")) +
    scale_fill_manual(values =c("#e7e5cc", "#9cc184", "#1e3d14"))
  
}

### Ki
plot_abundance_stat_Ki_high = function(phyloseq, title = "",
                                           Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#A5A492", "#447243")) +
    scale_fill_manual(values =c("#e7e5cc", "#447243"))
  
}

#Select variable of comparison
Go_comparison <- list(c("W1","W2"),c("W2","WM"),c("W1","WM"))
Ki_comparison <- list(c("W1","W3"))

##### Go_2020
plot_jki_seq2_seq4_taxa_Go_2020_RP_rel_high <- plot_abundance_stat_Go_high(ps_pgpr_Go_2020_RP, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2.5), label = c("0", "2.5", "5.0", "7.5", "10.0")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Go_comparison, label= "p", bracket.size = .3, size=3, label.y = c(4.0, 4.25, 4.5))
plot_jki_seq2_seq4_taxa_Go_2020_RP_rel_high

ggsave("plot_jki_seq2_seq4_taxa_Go_2020_RP_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 28, height = 12, units = "cm",dpi = 300)

##### Go 2021
plot_jki_seq2_seq4_taxa_Go_2021_RP_rel_high <- plot_abundance_stat_Go_high(ps_pgpr_Go_2021_RP, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 90,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2.5), label = c("0", "2.5", "5.0", "7.5", "10.0")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Go_comparison, label= "p", bracket.size = .3, size=3, label.y = c(4.0, 4.25, 4.5))
plot_jki_seq2_seq4_taxa_Go_2021_RP_rel_high

ggsave("plot_jki_seq2_seq4_taxa_Go_2021_RP_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 28, height = 16, units = "cm",dpi = 300)


##### Ki_2020
plot_jki_seq2_seq4_taxa_Ki_2020_RP_rel_high <- plot_abundance_stat_Ki_high(ps_pgpr_Ki_2020_RP, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2.5), label = c("0", "2.5", "5.0", "7.5", "10.0")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Ki_comparison, label= "p", bracket.size = .3, size=3, label.y = c(4.0))
plot_jki_seq2_seq4_taxa_Ki_2020_RP_rel_high

ggsave("plot_jki_seq2_seq4_taxa_Ki_2020_RP_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 28, height = 16, units = "cm",dpi = 300)


##### Ki_2021
plot_jki_seq2_seq4_taxa_Ki_2021_RP_rel_high <- plot_abundance_stat_Ki_high(ps_pgpr_Ki_2021_RP, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 90,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2.5), label = c("0", "2.5", "5.0", "7.5", "10.0")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Ki_comparison, label= "p", bracket.size = .3, size=3, label.y = c(6.0, 6.5, 7.0))
plot_jki_seq2_seq4_taxa_Ki_2021_RP_rel_high

ggsave("plot_jki_seq2_seq4_taxa_Ki_2021_RP_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 28, height = 16, units = "cm",dpi = 300)


###### LOW ABUNDANT TAXA

### By microhabitat
#ps_pgpr_Go_2020_BS1= subset_taxa(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
#ps_pgpr_Go_2020_RH1= subset_taxa(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
ps_pgpr_Go_2020_RP1= subset_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Gaiella" | Annotation == "Terrimonas")

#ps_pgpr_Go_2021_BS1= subset_taxa(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
#ps_pgpr_Go_2021_RH1= subset_taxa(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
ps_pgpr_Go_2021_RP1= subset_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")

#ps_pgpr_Ki_2020_BS1= subset_taxa(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
#ps_pgpr_Ki_2020_RH1= subset_taxa(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
ps_pgpr_Ki_2020_RP1= subset_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")

#ps_pgpr_Ki_2021_BS1= subset_taxa(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
#ps_pgpr_Ki_2021_RH1= subset_taxa(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")
ps_pgpr_Ki_2021_RP1= subset_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel, Annotation == "Pseudomonas" | Annotation == "Bacillus" | Annotation == "Niabella" | Annotation == "Terrimonas")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)
#### 2020
plot_abundance_stat_Go_low = function(phyloseq, title = "",
                                       Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#A5A492", "#9cc184", "#1e3d14")) +
    scale_fill_manual(values =c("#e7e5cc", "#9cc184", "#1e3d14"))
  
}

plot_abundance_stat_Ki_low = function(phyloseq, title = "",
                                       Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#A5A492", "#447243")) +
    scale_fill_manual(values =c("#e7e5cc", "#447243"))
  
}


##Subset by taxonomy, set comparisons, and run statistics for a specific otu taxa
### 2020
p_abundance_stat_Go_2020_RP1 <- plot_abundance_stat_Go_low(ps_pgpr_Go_2020_RP1, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  #scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, by = 0.01), label = c("0", "0.01", "0.02", "0.03")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 0.04, size=3) +
  labs(x="", y= "Relative abundance (%)") +
  #scale_y_break(c(10, 15)) +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  stat_compare_means(method = "t.test", comparisons = Go_comparison, label= "p", bracket.size = .3, size=3, label.y = c(2.4, 2.5, 2.6))+
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abundance_stat_Go_2020_RP1

ggsave("p_abundance_stat_Go_2020_RP1.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 28, height = 12, units = "cm",dpi = 300)

p_abundance_stat_Ki_2020_RP1 <- plot2_abundance_stat_Ki1(ps_pgpr_Ki_2020_RP1, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, by = 0.01), label = c("0", "0.01", "0.02", "0.03")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 0.04, size=3) +
  #scale_y_break(c(20, 35)) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abundance_stat_Ki_2020_RP1

ggsave("p_abundance_stat_Ki_2020_RP1.png", path = "~/Documents/R_analysis_2/jki_seq2_seq4/output/Stats_abundance/", width = 24, height = 16, units = "cm",dpi = 300)


#### 2021
p_abundance_stat_Go_2021_RP1 <- plot2_abundance_stat_Go1(ps_pgpr_Go_2021_RP1, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, by = 0.01), label = c("0", "0.01", "0.02", "0.03")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 0.04, size=3) +
  labs(x="", y= "Relative abundance (%)") +
  #scale_y_break(c(10, 15)) +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abundance_stat_Go_2021_RP1

ggsave("p_abundance_stat_Go_2021_RP1.png", path = "~/Documents/R_analysis_2/jki_seq2_seq4/output/Stats_abundance/", width = 30, height = 16, units = "cm",dpi = 300)

p_abundance_stat_Ki_2021_RP1 <- plot2_abundance_stat_Ki1(ps_pgpr_Ki_2021_RP1, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=22, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size = 22), legend.title = element_text(size = 18)) +
  #facet_wrap(~facet_label, scales = "free_x", ncol = 4, nrow = 2) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, by = 0.01), label = c("0", "0.01", "0.02", "0.03")) +
  #stat_compare_means(method = "kruskal.test", label= "p", label.y = 0.04, size=3) +
  #scale_y_break(c(20, 35)) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black"))
p_abundance_stat_Ki_2021_RP1

ggsave("p_abundance_stat_Ki_2021_RP1.png", path = "~/Documents/R_analysis_2/jki_seq2_seq4/output/Stats_abundance/", width = 18, height = 16, units = "cm",dpi = 300)


## The end! Have fun!
