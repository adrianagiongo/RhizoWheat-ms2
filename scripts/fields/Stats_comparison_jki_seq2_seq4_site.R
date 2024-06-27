###########
##### Correspondent to Figure 1C and Suppl Figure S1C
###########
#Go = "#1d5b65"   Go 2020 = "#216974"  Go 2021 = "#6baea5"
#Ki = "#A34828"   Ki 2020 = "#D1711F"  Ki 2021 = "#ebaf7a"

##BoxPlot comparison 
#Palette for Phyla
#Palette for Taxa 

#load packages
library("phyloseq")
library("microbiome")
library("ggplot2")
library("dplyr")
library("gridExtra")
library("ggpubr")
library("RColorBrewer")

#######################
#Subset taxa
#phylum
ps_jki_seq2_seq4_phylum_rel_high= subset_taxa(psO_jki_seq2_seq4_phylum_rel, Phylum == "Actinobacteriota" | Phylum == "Proteobacteria" )

ps_jki_seq2_seq4_phylum_rel_low= subset_taxa(psO_jki_seq2_seq4_phylum_rel, Phylum == "Firmicutes" | Phylum == "Acidobacteriota" | Phylum == "Bacteroidota" )

ps_jki_seq2_seq4_phylum_rel_verylow= subset_taxa(psO_jki_seq2_seq4_phylum_rel,  Phylum == "Verrucomicrobiota" | Phylum == "Patescibacteria" | Phylum =="Desulfobacterota" | Phylum == "Myxococcota")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)
plot_abundance_compare_Site = function(phyloseq, title = "",
                                            Facet ="Phylum", Color = "Site") {
  p1k_bac_treatment = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac_treatment = psmelt(p1k_bac_treatment)
  ggplot(data = mkk_bac_treatment, mapping = aes_string(x ="Site", y="Abundance", color = Color, fill = Color)) +
    theme_bw() +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, color= "black") + 
    #geom_jitter(size = 0.3, color = "black", position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}

# Set color palette
Site_colors <- c("Go" = "#1d5b65", "Ki" = "#A34828")
Site_year_colors <- c("Go_2020" = "#216974", "Go_2021" = "#6baea5", "Ki_2020" = "#D1711F", "Ki_2021" = "#ebaf7a")

##Subset by taxonomy, set comparisons, and run statistics for a specific otu taxa
#high abundance
plot_jki_seq2_seq4_phylum_rel_high <- plot_abundance_compare_Site(ps_jki_seq2_seq4_phylum_rel_high, Facet = "Phylum", Color = "Site_year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(20, 50), breaks = seq(20, 50, by = 10), label = c("20", "30", "40", "50")) +
  scale_color_manual(values = Site_year_colors) +
  scale_fill_manual(values = Site_year_colors) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 49, label.x = 1, size=3)
plot_jki_seq2_seq4_phylum_rel_high

ggsave("plot_jki_seq2_seq4_phylum_rel_high.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 9, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_jki_seq2_seq4_phylum_rel_low <- plot_abundance_compare_Site(ps_jki_seq2_seq4_phylum_rel_low, Facet = "Phylum", Color = "Site_year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 3), label = c("0", "3", "6", "9", "12", "15", "18")) +
  scale_color_manual(values = Site_year_colors) +
  scale_fill_manual(values = Site_year_colors) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 17.5, label.x = 1, size=3)
plot_jki_seq2_seq4_phylum_rel_low

ggsave("plot_jki_seq2_seq4_phylum_rel_low.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 12.6, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_jki_seq2_seq4_phylum_rel_verylow <- plot_abundance_compare_Site(ps_jki_seq2_seq4_phylum_rel_verylow, Facet = "Phylum", Color = "Site_year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1), label = c("0", "1", "2", "3", "4", "5")) +
  scale_color_manual(values = Site_year_colors) +
  scale_fill_manual(values = Site_year_colors) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 4.9, label.x = 1, size=3)
plot_jki_seq2_seq4_phylum_rel_verylow

ggsave("plot_jki_seq2_seq4_phylum_rel_verylow.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 17, height = 10, units = "cm", dpi = 300, device = "tiff")




####### I stopped here for the nt version.

###############################
### Rename strip names (facet)
#pgpr2_name <- c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "Allo-Neo-Para-Rhi", 'Bosea' = "Bosea", 'Phyllobacterium' = "Phyllobacterium", 'Rhizobiales_Incertae_Sedis' = "Rhizobiales", 'Shinella' = "Shinella")

#ps_pgpr2= subset_taxa(psO_jki_seq1_annotation_rel, Annotation == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" | Annotation == "Bosea" | Annotation == "Phyllobacterium" | Annotation == "Rhizobiales_Incertae_Sedis" | Annotation == "Shinella")


#######################
##Agglomerating (using tax_glom()) taxa per sample at the appropriated taxonomic rank
#for psO_jki_seq2_seq4
colnames(tax_table(psO_jki_seq2_seq4_Go_filt))
psO_jki_seq2_seq4_Go_filt_phylum <- tax_glom(psO_jki_seq2_seq4_Go_filt, taxrank = "Phylum")
ntaxa(psO_jki_seq2_seq4_Go_filt); ntaxa(psO_jki_seq2_seq4_Go_filt_phylum)

colnames(tax_table(psO_jki_seq2_seq4_Ki_filt))
psO_jki_seq2_seq4_Ki_filt_phylum <- tax_glom(psO_jki_seq2_seq4_Ki_filt, taxrank = "Phylum")
ntaxa(psO_jki_seq2_seq4_Ki_filt); ntaxa(psO_jki_seq2_seq4_Ki_filt_phylum)

### Transform to relative abundance
psO_jki_seq2_seq4_Go_filt_phylum_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_filt_phylum, function(x) (x*100)/sum(x))
psO_jki_seq2_seq4_Go_filt_phylum_rel
head(otu_table(psO_jki_seq2_seq4_Go_filt_phylum_rel))

psO_jki_seq2_seq4_Ki_filt_phylum_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_filt_phylum, function(x) (x*100)/sum(x))
psO_jki_seq2_seq4_Ki_filt_phylum_rel
head(otu_table(psO_jki_seq2_seq4_Ki_filt_phylum_rel))

#Subset taxa
#phylum
ps_jki_seq2_seq4_Go_filt_phylum_rel_high= subset_taxa(psO_jki_seq2_seq4_Go_filt_phylum_rel, Phylum == "Actinobacteriota" | Phylum == "Proteobacteria")
ps_jki_seq2_seq4_Go_filt_phylum_rel_low= subset_taxa(psO_jki_seq2_seq4_Go_filt_phylum_rel, Phylum == "Chloroflexi" | Phylum == "Acidobacteriota" | Phylum == "Gemmatimonadota" | Phylum == "Firmicutes")
ps_jki_seq2_seq4_Go_filt_phylum_rel_verylow= subset_taxa(psO_jki_seq2_seq4_Go_filt_phylum_rel,  Phylum == "Nitrospirota" | Phylum == "Verrucomicrobiota" | Phylum == "Patescibacteria")

ps_jki_seq2_seq4_Ki_filt_phylum_rel_high= subset_taxa(psO_jki_seq2_seq4_Ki_filt_phylum_rel, Phylum == "Actinobacteriota" | Phylum == "Proteobacteria")
ps_jki_seq2_seq4_Ki_filt_phylum_rel_low= subset_taxa(psO_jki_seq2_seq4_Ki_filt_phylum_rel, Phylum == "Chloroflexi" | Phylum == "Acidobacteriota" | Phylum == "Gemmatimonadota" | Phylum == "Firmicutes")
ps_jki_seq2_seq4_Ki_filt_phylum_rel_verylow= subset_taxa(psO_jki_seq2_seq4_Ki_filt_phylum_rel,  Phylum == "Nitrospirota" | Phylum == "Verrucomicrobiota" | Phylum == "Patescibacteria")

##Create a function to plot abundances for a subset of samples (Subset Annotation and facet Annotation)

plot_abundance_compare_Year = function(phyloseq, title = "",
                                       Facet ="Phylum", Color = "Year") {
  p1k_bac_treatment = subset_taxa(phyloseq, Kingdom %in% c("Bacteria"))
  mkk_bac_treatment = psmelt(p1k_bac_treatment)
  ggplot(data = mkk_bac_treatment, mapping = aes_string(x ="Year", y="Abundance", color = Color, fill = Color)) +
    theme_bw() +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color= "black") + geom_point(size = 0.3, color = "black", position = position_jitter(width = 0.2)) +
    facet_wrap(facets = Facet, nrow = 1, strip.position = "top")
}


##Subset by taxonomy, set comparisons, and run statistics for a specific otu taxa
#high abundance
plot_jki_seq2_seq4_Go_filt_phylum_rel_high <- plot_abundance_compare_Year(ps_jki_seq2_seq4_Go_filt_phylum_rel_high, Facet = "Phylum", Color = "Year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20), label = c("0", "20", "40", "60", "80")) +
  scale_color_manual(values = c("#B47474", "#831818")) +
  scale_fill_manual(values = c("#B47474", "#831818")) +
  scale_x_discrete("Year", labels = c("Y2020" = "2020", "Y2021" = "2021")) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 70, label.x = 1, size=4)
plot_jki_seq2_seq4_Go_filt_phylum_rel_high

ggsave("plot_jki_seq2_seq4_Go_filt_phylum_rel_high.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 9, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_jki_seq2_seq4_Go_filt_phylum_rel_low <- plot_abundance_compare_Year(ps_jki_seq2_seq4_Go_filt_phylum_rel_low, Facet = "Phylum", Color = "Year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), label = c("0", "5", "10", "15", "20")) +
  scale_color_manual(values = c("#B47474", "#831818")) +
  scale_fill_manual(values = c("#B47474", "#831818")) +
  scale_x_discrete("Year", labels = c("Y2020" = "2020", "Y2021" = "2021")) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 18, label.x = 1, size=4)
plot_jki_seq2_seq4_Go_filt_phylum_rel_low

ggsave("plot_jki_seq2_seq4_Go_filt_phylum_rel_low.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 17, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_jki_seq2_seq4_Go_filt_phylum_rel_verylow <- plot_abundance_compare_Year(ps_jki_seq2_seq4_Go_filt_phylum_rel_verylow, Facet = "Phylum", Color = "Year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1), label = c("0", "1", "2", "3")) +
  scale_color_manual(values = c("#B47474", "#831818")) +
  scale_fill_manual(values = c("#B47474", "#831818")) +
  scale_x_discrete("Year", labels = c("Y2020" = "2020", "Y2021" = "2021")) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 2.8, label.x = 1, size=4)
plot_jki_seq2_seq4_Go_filt_phylum_rel_verylow

ggsave("plot_jki_seq2_seq4_Go_filt_phylum_rel_verylow.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 14, height = 10, units = "cm", dpi = 300, device = "tiff")




#high abundance
plot_jki_seq2_seq4_Ki_filt_phylum_rel_high <- plot_abundance_compare_Year(ps_jki_seq2_seq4_Ki_filt_phylum_rel_high, Facet = "Phylum", Color = "Year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20), label = c("0", "20", "40", "60", "80")) +
  scale_color_manual(values = c("#f69c8e", "#f05b43")) +
  scale_fill_manual(values = c("#f69c8e", "#f05b43")) +
  scale_x_discrete("Year", labels = c("Y2020" = "2020", "Y2021" = "2021")) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 70, label.x = 1, size=4)
plot_jki_seq2_seq4_Ki_filt_phylum_rel_high

ggsave("plot_jki_seq2_seq4_Ki_filt_phylum_rel_high.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 9, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_jki_seq2_seq4_Ki_filt_phylum_rel_low <- plot_abundance_compare_Year(ps_jki_seq2_seq4_Ki_filt_phylum_rel_low, Facet = "Phylum", Color = "Year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5), label = c("0", "5", "10", "15", "20")) +
  scale_color_manual(values = c("#f69c8e", "#f05b43")) +
  scale_fill_manual(values = c("#f69c8e", "#f05b43")) +
  scale_x_discrete("Year", labels = c("Y2020" = "2020", "Y2021" = "2021")) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 18, label.x = 1, size=4)
plot_jki_seq2_seq4_Ki_filt_phylum_rel_low

ggsave("plot_jki_seq2_seq4_Ki_filt_phylum_rel_low.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 17, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_jki_seq2_seq4_Ki_filt_phylum_rel_verylow <- plot_abundance_compare_Year(ps_jki_seq2_seq4_Ki_filt_phylum_rel_verylow, Facet = "Phylum", Color = "Year") +
  ylab("Relative abundance (%)") +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_text(size=14, color = "black")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 10, face="italic", color = "black")) +
  theme(axis.text.x = element_text(size = 14, colour = "black"), axis.title.y = element_text(size = 14, colour = "black")) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1), label = c("0", "1", "2", "3")) +
  scale_color_manual(values = c("#f69c8e", "#f05b43")) +
  scale_fill_manual(values = c("#f69c8e", "#f05b43")) +
  scale_x_discrete("Year", labels = c("Y2020" = "2020", "Y2021" = "2021")) +
  #stat_compare_means(comparisons = Comparison_treatment_exp3, method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), bracket.size = .3, size=3, label.y = c(0.04, 0.045, 0.05))+
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 2.8, label.x = 1, size=4)
plot_jki_seq2_seq4_Ki_filt_phylum_rel_verylow

ggsave("plot_jki_seq2_seq4_Ki_filt_phylum_rel_verylow.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Stats_comparison_jki_seq2_seq4/", width = 14, height = 10, units = "cm", dpi = 300, device = "tiff")

## I stopped here : )
