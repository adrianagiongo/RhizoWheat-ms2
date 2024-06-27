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
length(get_taxa_unique(psO_jki_seq10_rarefied, taxonomic.rank = "Phylum"))

### Subset taxa by microhabitat
ps_rarefied_pgpr_high_Go= subset_taxa(psO_jki_seq10_rarefied_Go_filt_phylum_rel, Phylum == "Actinobacteriota" | Phylum == "Proteobacteria")
ps_rarefied_pgpr_high_Ki= subset_taxa(psO_jki_seq10_rarefied_Ki_filt_phylum_rel, Phylum == "Actinobacteriota" | Phylum == "Proteobacteria")

ps_rarefied_pgpr_low_Go= subset_taxa(psO_jki_seq10_rarefied_Go_filt_phylum_rel, Phylum == "Firmicutes" | Phylum == "Acidobacteriota" | Phylum == "Bacteroidota")
ps_rarefied_pgpr_low_Ki= subset_taxa(psO_jki_seq10_rarefied_Ki_filt_phylum_rel, Phylum == "Firmicutes" | Phylum == "Acidobacteriota" | Phylum == "Bacteroidota")


##Create a function to plot abundances for a subset of samples (Subset Phylum and facet Phylum)

###### 2020
plot_abundance_stat_rarefied_Go = function(phyloseq, title = "",
                                           Facet ="Phylum", Color = "Site_rot") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Site_rot_tre", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "right") +
    scale_color_manual(values =c("#6baea5", "#216974")) +
    scale_fill_manual(values =c("#6baea5", "#216974"))
  
}

plot_abundance_stat_rarefied_Ki = function(phyloseq, title = "",
                                           Facet ="Phylum", Color = "Site_rot") {
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
plot_jki_seq10_rarefied_phylum_Go_rel_high <- plot_abundance_stat_rarefied_Go(ps_rarefied_pgpr_high_Go, Facet = "Phylum", Color = "Site_rot") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=20, angle = 90,vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 25), label = c("0", "25", "50")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = treat_comp_Go, label= "p", bracket.size = .3, size=3, label.y = c(45,45))
plot_jki_seq10_rarefied_phylum_Go_rel_high

ggsave("plot_jki_seq10_rarefied_phylum_Go_rel_high.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Stats_comparison_jki_seq10/", width = 20, height = 12, units = "cm",dpi = 300)

##### Ki
plot_jki_seq10_rarefied_phylum_Ki_rel_high <- plot_abundance_stat_rarefied_Ki(ps_rarefied_pgpr_high_Ki, Facet = "Phylum", Color = "Site_rot") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=20, angle = 90,vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 25), label = c("0", "25", "50")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = treat_comp_Ki, label= "p", bracket.size = .3, size=3, label.y = c(45, 45))
plot_jki_seq10_rarefied_phylum_Ki_rel_high

ggsave("plot_jki_seq10_rarefied_phylum_Ki_rel_high.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Stats_comparison_jki_seq10/", width = 20, height = 12, units = "cm",dpi = 300)




##Subset by taxonomy, set comparisons, and run statistics for a specific otu taxa

##### Go
plot_jki_seq10_rarefied_phylum_Go_rel_low <- plot_abundance_stat_rarefied_Go(ps_rarefied_pgpr_low_Go, Facet = "Phylum", Color = "Site_rot") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=20, angle = 90,vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 10), label = c("0", "10", "20")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = treat_comp_Go, label= "p", bracket.size = .3, size=3, label.y = c(18, 18))
plot_jki_seq10_rarefied_phylum_Go_rel_low

ggsave("plot_jki_seq10_rarefied_phylum_Go_rel_low.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Stats_comparison_jki_seq10/", width = 26, height = 12, units = "cm",dpi = 300)

##### Ki
plot_jki_seq10_rarefied_phylum_Ki_rel_low <- plot_abundance_stat_rarefied_Ki(ps_rarefied_pgpr_low_Ki, Facet = "Phylum", Color = "Site_rot") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=20, angle = 90,vjust = 0.5, hjust = 1), axis.text.y = element_text(size=18)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size = 18)) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 10), label = c("0", "10", "20")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 18, face="italic")) +
  theme(axis.title.x = element_text(size = 18, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = treat_comp_Ki, label= "p", bracket.size = .3, size=3, label.y = c(18, 18))
plot_jki_seq10_rarefied_phylum_Ki_rel_low

ggsave("plot_jki_seq10_rarefied_phylum_Ki_rel_low.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Stats_comparison_jki_seq10/", width = 26, height = 12, units = "cm",dpi = 300)


# The end! Enjoy the sun!


