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

#####################
###### MOST ABUNDANT TAXA
#####################
###Select groups and subset data from the cleaned dataset

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq2_seq4_rarefied, taxonomic.rank = "Annotation"))

### Subset taxa by microhabitat
ps_rarefied_Go_2020_BS_T2 = subset_taxa(psO_jki_seq2_seq4_rarefied_Go_2020_BS_T2_filt_annotation_rel, Annotation == "ANPR" | Annotation == "Bacillus" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas")
ps_rarefied_Go_2020_RH_T2 = subset_taxa(psO_jki_seq2_seq4_rarefied_Go_2020_RH_T2_filt_annotation_rel, Annotation == "ANPR" | Annotation == "Bacillus" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas")
ps_rarefied_Go_2020_RP_T2 = subset_taxa(psO_jki_seq2_seq4_rarefied_Go_2020_RP_T2_filt_annotation_rel, Annotation == "ANPR" | Annotation == "Bacillus" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas")

ps_rarefied_Ki_2020_BS_T2 = subset_taxa(psO_jki_seq2_seq4_rarefied_Ki_2020_BS_T2_filt_annotation_rel, Annotation == "ANPR" | Annotation == "Bacillus" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas")
ps_rarefied_Ki_2020_RH_T2 = subset_taxa(psO_jki_seq2_seq4_rarefied_Ki_2020_RH_T2_filt_annotation_rel, Annotation == "ANPR" | Annotation == "Bacillus" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas")
ps_rarefied_Ki_2020_RP_T2 = subset_taxa(psO_jki_seq2_seq4_rarefied_Ki_2020_RP_T2_filt_annotation_rel, Annotation == "ANPR" | Annotation == "Bacillus" | Annotation == "Sphingomonas" | Annotation == "Pseudomonas")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)

plot_abundance_stat_Go_high = function(phyloseq, title = "",
                                           Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "none") +
    scale_color_manual(values =c("#A5A492", "#9cc184", "#1e3d14")) +
    scale_fill_manual(values =c("#e7e5cc", "#9cc184", "#1e3d14"))
  
}

plot_abundance_stat_Ki_high = function(phyloseq, title = "",
                                           Facet ="Annotation", Color = "Rotation") {
  p1k_bac = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac = psmelt(p1k_bac)
  ggplot(data = mkk_bac, mapping = aes_string(x ="Rotation", y="Abundance", color = Color, fill = Color)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width=0.9) + 
    facet_wrap(facets = Facet, nrow = 1) +
    theme(legend.position = "none") +
    scale_color_manual(values =c("#A5A492", "#447243")) +
    scale_fill_manual(values =c("#e7e5cc", "#447243"))
  
}

#Select variable of comparison
Go_2020_comparison_rotation <- list(c("W1","W2"),c("W1","WM"),c("W2","WM"))
Ki_2020_comparison_rotation <- list(c("W1","W3"))

##Subset by taxonomy, set comparisons, and run statistics for a specific otu taxa

##### Go BS
plot_jki_seq2_seq4_rarefied_taxa_Go_2020_BS_T2_rel_high <- plot_abundance_stat_Go_high(ps_rarefied_Go_2020_BS_T2, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), label = c("0", "2", "4", "6", "8")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 20, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Go_2020_comparison_rotation, label= "p", bracket.size = .3, size=4, label.y = c(7, 7.6))
plot_jki_seq2_seq4_rarefied_taxa_Go_2020_BS_T2_rel_high

ggsave("plot_jki_seq2_seq4_rarefied_taxa_Go_2020_BS_T2_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 26, height = 10, units = "cm",dpi = 300)

##### Go RH
plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RH_T2_rel_high <- plot_abundance_stat_Go_high(ps_rarefied_Go_2020_RH_T2, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 20, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Go_2020_comparison_rotation, label= "p", bracket.size = .3, size=4, label.y = c(5, 5.6))
plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RH_T2_rel_high

ggsave("plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RH_T2_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/",  width = 26, height = 10, units = "cm",dpi = 300)

##### Go RP
plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RP_T2_rel_high <- plot_abundance_stat_Go_high(ps_rarefied_Go_2020_RP_T2, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 6.3), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 20, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Go_2020_comparison_rotation, label= "p", bracket.size = .3, size=4, label.y = c(5, 5.6))
plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RP_T2_rel_high

ggsave("plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RP_T2_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/",  width = 26, height = 10, units = "cm",dpi = 300)




##### Ki BS
plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_BS_T2_rel_high <- plot_abundance_stat_Ki_high(ps_rarefied_Ki_2020_BS_T2, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2), label = c("0", "2", "4", "6", "8")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 20, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Ki_2020_comparison_rotation, label= "p", bracket.size = .3, size=4, label.y = c(7.4))
plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_BS_T2_rel_high

ggsave("plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_BS_T2_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 26, height = 10, units = "cm",dpi = 300)

##### Ki RH
plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RH_T2_rel_high <- plot_abundance_stat_Ki_high(ps_rarefied_Ki_2020_RH_T2, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 20, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  stat_compare_means(method = "t.test", comparisons = Ki_2020_comparison_rotation, label= "p", bracket.size = .3, size=4, label.y = c(5.6))
plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RH_T2_rel_high

ggsave("plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RH_T2_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/",  width = 26, height = 10, units = "cm",dpi = 300)

##### Ki RP
plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RP_T2_rel_high <- plot_abundance_stat_Ki_high(ps_rarefied_Ki_2020_RP_T2, Facet = "Annotation", Color = "Rotation") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(0, 6.3), breaks = seq(0, 6, by = 2), label = c("0", "2", "4", "6")) +
  labs(x="", y= "Relative abundance (%)") +
  theme(strip.text.x = element_text(size = 20, face="italic")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  #stat_compare_means(method = "kruskal", label= "p", label.y = 8, label.x = 1, size=5) +
  stat_compare_means(method = "t.test", comparisons = Ki_2020_comparison_rotation, label= "p", bracket.size = .3, size=4, label.y = c(5.6))
#stat_compare_means(method = "t.test", comparisons = Ki_comparison_stage_rotation, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(6.0, 6.5, 7.0, 6.0, 6.5, 7.0))

plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RP_T2_rel_high

ggsave("plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RP_T2_rel_high.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/",  width = 26, height = 10, units = "cm",dpi = 300)




###########################################
#arrange multiple ggplots in one single figure ALL
stats_comp_plots_all <- ggarrange(plot_jki_seq2_seq4_rarefied_taxa_Go_2020_BS_T2_rel_high,
                                    plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_BS_T2_rel_high,
                                    plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RH_T2_rel_high,
                                    plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RH_T2_rel_high,
                                    plot_jki_seq2_seq4_rarefied_taxa_Go_2020_RP_T2_rel_high,
                                    plot_jki_seq2_seq4_rarefied_taxa_Ki_2020_RP_T2_rel_high,
                                    ncol =2, nrow =3, legend="right", widths = c(1.1, 1), heights = c(1, 1), common.legend = TRUE, label.y = 1)
stats_comp_plots_all

ggsave("stats_comp_plots_all.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 58, height = 36, units = "cm",dpi = 300, device = "svg")
###########################################

# The end!  Tchüss!!




