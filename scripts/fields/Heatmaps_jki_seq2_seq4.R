##Creating heatmap for selected dataset
rmarkdown::render("~/Documents/R_analysis/jki_seq2_seq4/scripts_jki_seq2_seq4/Heatmaps_jki_seq2_seq4.R")
Site_colors <- c("Go" = "#1d5b65", "Ki" = "#A34828")
##load libraries
library("phyloseq")
library("ggplot2")
library("RColorBrewer")

## Heatmaps with annotations with + 1% 
##for psO_jki_seq2_Go_filt_annotation_rel
#subset annotation with >1% abundance
totalGo_2020_RP_T2 = sample_sums(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel)
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_abund1 <- filter_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel, function(x) sum(x > totalGo_2020_RP_T2*0.01) > 0, TRUE)
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_abund1

totalGo_2021_RP_T2 = sample_sums(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel)
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_abund1 <- filter_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel, function(x) sum(x > totalGo_2021_RP_T2*0.01) > 0, TRUE)
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_abund1


totalKi_2020_RP_T2 = sample_sums(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel)
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_abund1 <- filter_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel, function(x) sum(x > totalKi_2020_RP_T2*0.01) > 0, TRUE)
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_abund1

totalKi_2021_RP_T2 = sample_sums(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel)
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_abund1 <- filter_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel, function(x) sum(x > totalKi_2021_RP_T2*0.01) > 0, TRUE)
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_abund1


##Using plot_Heatmap()
###Heatmap for Go
heatmap_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_abund1 = plot_heatmap(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_abund1, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#1d5b65",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 18), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 20, face="bold")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_abund1

ggsave("heatmap_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel_abund1.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 19, height = 25, units = "cm",dpi = 300)


heatmap_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_abund1 = plot_heatmap(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_abund1, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#57896a",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 20, face="bold")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_abund1

ggsave("heatmap_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel_abund1.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 22, height = 30, units = "cm",dpi = 300)


#Heatmap for Ki
heatmap_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_abund1 = plot_heatmap(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_abund1, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#A34828",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 20, face="bold")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_abund1

ggsave("heatmap_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel_abund1.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 20, height = 30, units = "cm",dpi = 300)


heatmap_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_abund1 = plot_heatmap(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_abund1, sample.order = "Stage",sample.label = "Sample_name", taxa.label = "Annotation",taxa.order="Annotation",low="white",high="#57896a",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free") +
  labs(fill="Relative Abundance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face = "italic")) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
  theme(strip.text.x = element_text(size = 20, face="bold")) +
  #scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_abund1

ggsave("heatmap_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel_abund1.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4/", width = 22, height = 30, units = "cm",dpi = 300)


## Enjoy the day!  : )
