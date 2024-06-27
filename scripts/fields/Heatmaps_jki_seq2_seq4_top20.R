##Creating heatmap for selected dataset
rmarkdown::render("~/Documents/R_analysis/jki_seq2_seq4/scripts_jki_seq2_seq4/Heatmaps_jki_seq2_seq4.R")
Site_colors <- c("Go" = "#1d5b65", "Ki" = "#A34828")
##load libraries
library("phyloseq")
library("microbiome")
library("ggplot2")
library("RColorBrewer")
library("gridExtra")

## Heatmaps with annotations with TOP 20
##for psO_jki_seq2_Go_filt_annotation_rel
#List top abundant taxa
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10 <- top_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel, n = 10)
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10

#Subset top abundant taxa
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa <- prune_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10, psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel)

#Create tables
df_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa),otu_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa))
write.csv(df_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa.csv")

##Using plot_Heatmap()
###Heatmap for Go
heatmap_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa = plot_heatmap(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white", high="#1d5b65", na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 22), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 22)) +
  scale_fill_gradientn (limits = c(0, 7), breaks = seq(0, 7, by = 3), label = c("0", "3", "6"), colours = c("white", "#1d5b65")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa

ggsave("heatmap_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_top10_taxa.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 18, height = 10, units = "cm",dpi = 300)


## Heatmaps with annotations with TOP 20
##for psO_jki_seq2_Ki_filt_annotation_rel
#List top abundant taxa
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10 <- top_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel, n = 10)
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10

#Subset top abundant taxa
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa <- prune_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10, psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel)

#Create tables
df_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa),otu_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa))
write.csv(df_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa.csv")

##Using plot_Heatmap()
###Heatmap for Ki
heatmap_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa = plot_heatmap(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#A34828",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=22, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 22), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 22)) +
  scale_fill_gradientn (limits = c(0, 7), breaks = seq(0, 7, by = 3), label = c("0", "3", "6"), colours = c("white", "#A34828")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa

ggsave("heatmap_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_top10_taxa.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 16, height = 10, units = "cm",dpi = 300)


## Enjoy the day!  : )


#List top abundant taxa
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10 <- top_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel, n = 10)
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10

#Subset top abundant taxa
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa <- prune_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10, psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel)

#Create tables
df_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa),otu_table(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa))
write.csv(df_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa.csv")

##Using plot_Heatmap()
###Heatmap for Go
heatmap_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa = plot_heatmap(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#1d5b65",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 22), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 22, face="bold")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa

ggsave("heatmap_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_top10_taxa.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 16, height = 10, units = "cm",dpi = 300)


## Heatmaps with annotations with TOP 20
##for psO_jki_seq2_Ki_filt_annotation_rel
#List top abundant taxa
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10 <- top_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel, n = 10)
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10

#Subset top abundant taxa
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa <- prune_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10, psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel)

#Create tables
df_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa),otu_table(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa))
write.csv(df_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa.csv")

##Using plot_Heatmap()
###Heatmap for Ki
heatmap_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa = plot_heatmap(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#A34828",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=22, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 22), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 22, face="bold")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa

ggsave("heatmap_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_top10_taxa.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 16, height = 10, units = "cm",dpi = 300)

