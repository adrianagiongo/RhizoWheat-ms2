##Creating heatmap for selected dataset

##load libraries
library("phyloseq")
library("microbiome")
library("ggplot2")
library("RColorBrewer")
library("gridExtra")

#### All soils together
####
#Aggregate taxa
psO_jki_seq10_rarefied_agreg_annotation <- aggregate_taxa(psO_jki_seq10_rarefied, level = "Annotation")

#Transform abundance to relative abundance
psO_jki_seq10_rarefied_agreg_annotation_rel <- transform_sample_counts(psO_jki_seq10_rarefied_agreg_annotation, function(x) (x*100)/sum(x))
psO_jki_seq10_rarefied_agreg_annotation_rel

#List top abundant taxa
psO_jki_seq10_rarefied_agreg_annotation_rel_top20 <- top_taxa(psO_jki_seq10_rarefied_agreg_annotation_rel, n = 20)
psO_jki_seq10_rarefied_agreg_annotation_rel_top20

#Subset top abundant taxa
psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa <- prune_taxa(psO_jki_seq10_rarefied_agreg_annotation_rel_top20, psO_jki_seq10_rarefied_agreg_annotation_rel)

#Create tables
df_psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa <- data.frame(tax_table(psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa),otu_table(psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa))
write.csv(df_psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa.csv")

###Heatmap plot
##Using plot_Heatmap()
heatmap_psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa = plot_heatmap(psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa, sample.order = "Site_rot_tre", sample.label = "Site_rot_tre", taxa.order = "Annotation", taxa.label = "Annotation", verbose = FALSE, low="white",high="black",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Site, scales = "free", space="free_x") +
  labs(fill="Relative Abundance (%)") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 24, angle = 90,vjust = 0.5, hjust = 1, colour = "black"), axis.text.y = element_text(size=24, face = "italic", colour = "black")) +
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 24), legend.text.align = 1) +
  theme(strip.text.x = element_text(size = 24, face="bold")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa

ggsave("heatmap_psO_jki_seq10_rarefied_agreg_annotation_rel_top20_taxa.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Heatmap_jki_seq10/", width = 60, height = 30, units = "cm",dpi = 300)

## distance = "bray"


#### Each soil individually

## Go
#Aggregate taxa
psO_jki_seq10_rarefied_Go_filt_agreg_annotation <- aggregate_taxa(psO_jki_seq10_rarefied_Go_filt, level = "Annotation")

#Transform abundance to relative abundance
psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel <- transform_sample_counts(psO_jki_seq10_rarefied_Go_filt_agreg_annotation, function(x) (x*100)/sum(x))
psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel

#List top abundant taxa
psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20 <- top_taxa(psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel, n = 20)
psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20

#Subset top abundant taxa
psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa <- prune_taxa(psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20, psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel)

#Create tables
df_psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa <- data.frame(tax_table(psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa),otu_table(psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa))
write.csv(df_psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa.csv")

###Heatmap plot
##Using plot_Heatmap()
heatmap_psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa = plot_heatmap(psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa, sample.order = "Site_rot_tre", sample.label = "Site_rot_tre", taxa.order = "Annotation", taxa.label = "Annotation", verbose = FALSE, low="white",high="#486c71",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free", space="free_x") +
  labs(fill="(%)") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 24, angle = 90,vjust = 0.5, hjust = 1, colour = "black"), axis.text.y = element_text(size=24, face = "italic", colour = "black")) +
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 24), legend.text.align = 1) +
  theme(strip.text.x = element_text(size = 24, face="bold")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa

ggsave("heatmap_psO_jki_seq10_rarefied_Go_filt_agreg_annotation_rel_top20_taxa.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Heatmap_jki_seq10/", width = 32, height = 35, units = "cm",dpi = 300)



## Ki
#Aggregate taxa
psO_jki_seq10_rarefied_Ki_filt_agreg_annotation <- aggregate_taxa(psO_jki_seq10_rarefied_Ki_filt, level = "Annotation")

#Transform abundance to relative abundance
psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel <- transform_sample_counts(psO_jki_seq10_rarefied_Ki_filt_agreg_annotation, function(x) (x*100)/sum(x))
psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel

#List top abundant taxa
psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20 <- top_taxa(psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel, n = 20)
psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20

#Subset top abundant taxa
psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa <- prune_taxa(psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20, psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel)

#Create tables
df_psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa <- data.frame(tax_table(psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa),otu_table(psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa))
write.csv(df_psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa.csv")

###Heatmap plot
##Using plot_Heatmap()
heatmap_psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa = plot_heatmap(psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa, sample.order = "Site_rot_tre", sample.label = "Site_rot_tre", taxa.order = "Annotation", taxa.label = "Annotation", verbose = FALSE, low="white",high="#c33935",na.value="white", trans = NULL) +
  theme_bw() +
  facet_grid(~Rotation, scales = "free", space="free_x") +
  labs(fill="(%)") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 24, angle = 90,vjust = 0.5, hjust = 1, colour = "black"), axis.text.y = element_text(size=24, face = "italic", colour = "black")) +
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 24), legend.text.align = 1) +
  theme(strip.text.x = element_text(size = 24, face="bold")) +
  geom_tile(colour="gray",size=0.15)
heatmap_psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa

ggsave("heatmap_psO_jki_seq10_rarefied_Ki_filt_agreg_annotation_rel_top20_taxa.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Heatmap_jki_seq10/", width = 32.5, height = 35, units = "cm",dpi = 300)


#The end! Enjoy your day! : )
