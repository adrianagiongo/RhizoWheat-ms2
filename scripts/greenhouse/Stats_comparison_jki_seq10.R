##BoxPlot comparison taxa

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
length(get_taxa_unique(psO_jki_seq10, taxonomic.rank = "Annotation"))

### Subset taxa by microhabitat
ps_pgpr_Go= subset_taxa(psO_jki_seq10_Go_filt_annotation_rel, Annotation == "Ensifer" | Annotation == "Pedobacter" | Annotation == "Microbacterium" | Annotation == "Pseudomonas" | Annotation == "Sphingomonas" | Annotation == "Flavobacterium" | Annotation == "Bacillus")

ps_pgpr_Ki= subset_taxa(psO_jki_seq10_Ki_filt_annotation_rel, Annotation == "Ensifer" | Annotation == "Pedobacter" | Annotation == "Microbacterium" | Annotation == "Pseudomonas" | Annotation == "Sphingomonas" | Annotation == "Flavobacterium" | Annotation == "Bacillus")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)

###### 2020
plot_abundance_stat_Go_high = function(phyloseq, title = "",
                                           Facet ="Annotation", Color = "Site_rot") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Site_rot_tre", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#6baea5", "#216974")) +
    scale_fill_manual(values =c("#6baea5", "#216974"))
  
}

plot_abundance_stat_Ki_high = function(phyloseq, title = "",
                                           Facet ="Annotation", Color = "Site_rot") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Site_rot_tre", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#ebaf7a", "#D1711F")) +
    scale_fill_manual(values =c("#ebaf7a", "#D1711F"))
  
}

#Select variable of comparison
treat_comp_Go <- list(c("Go_W1_C", "Go_W1_Ggt"), c("Go_WM_C", "Go_WM_Ggt"))
treat_comp_Ki <- list(c("Ki_W1_C", "Ki_W1_Ggt"), c("Ki_W3_C", "Ki_W3_Ggt"))

##Subset by taxonomy, set comparisons, and run statistics for a specific otu taxa

##### Go
plot_jki_seq10_taxa_Go_rel_high <- plot_abundance_stat_Go_high(ps_pgpr_Go, Facet = "Annotation", Color = "Site_rot") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=20, angle = 90,vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 2.5), label = c("0", "2.5", "5")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "wilcox.test", comparisons = treat_comp_Go, label= "p", bracket.size = .3, size=3, label.y = c(4.5,4.5))
plot_jki_seq10_taxa_Go_rel_high

ggsave("plot_jki_seq10_taxa_Go_rel_high.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Stats_comparison_jki_seq10/", width = 42, height = 14, units = "cm",dpi = 300)

##### Ki
plot_jki_seq10_taxa_Ki_rel_high <- plot_abundance_stat_Ki_high(ps_pgpr_Ki, Facet = "Annotation", Color = "Site_rot") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=20, angle = 90,vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 2.5), label = c("0", "2.5", "5")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "wilcox.test", comparisons = treat_comp_Ki, label= "p", bracket.size = .3, size=3, label.y = c(4.5,4.5))
plot_jki_seq10_taxa_Ki_rel_high

ggsave("plot_jki_seq10_taxa_Ki_rel_high.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Stats_comparison_jki_seq10/", width = 42, height = 14, units = "cm",dpi = 300)


# The end! Enjoy the sun!


