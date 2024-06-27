##Creating MDS plot for selected dataset
#load package
library("phyloseq")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")
library("ggpubr")

################################
set.seed(2022)

# Set color palette
Location_colors <- c("Go" = "#831818", "Ki" = "#f05b43")


###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors
#(1)Location = Go or Ki
#(2)Year = 2020 or 2021 
#(3)Rotation = W1, W2, WM (Go) or W1, W3 (Ki)

#### Separated by location, stage and microhabitat 
psO_jki_seq2_seq4_Go_T2_BS <- subset_samples(psO_jki_seq2_seq4_Go_T2_filt, Microhabitat == "BS")
psO_jki_seq2_seq4_Go_T2_BS

psO_jki_seq2_seq4_Go_T2_RH <- subset_samples(psO_jki_seq2_seq4_Go_T2_filt, Microhabitat == "RH")
psO_jki_seq2_seq4_Go_T2_RH

psO_jki_seq2_seq4_Go_T2_RP <- subset_samples(psO_jki_seq2_seq4_Go_T2_filt, Microhabitat == "RP")
psO_jki_seq2_seq4_Go_T2_RP

psO_jki_seq2_seq4_Go_T3_BS <- subset_samples(psO_jki_seq2_seq4_Go_T3_filt, Microhabitat == "BS")
psO_jki_seq2_seq4_Go_T3_BS

psO_jki_seq2_seq4_Go_T3_RH <- subset_samples(psO_jki_seq2_seq4_Go_T3_filt, Microhabitat == "RH")
psO_jki_seq2_seq4_Go_T3_RH

psO_jki_seq2_seq4_Go_T3_RP <- subset_samples(psO_jki_seq2_seq4_Go_T3_filt, Microhabitat == "RP")
psO_jki_seq2_seq4_Go_T3_RP


psO_jki_seq2_seq4_Ki_T2_BS <- subset_samples(psO_jki_seq2_seq4_Ki_T2_filt, Microhabitat == "BS")
psO_jki_seq2_seq4_Ki_T2_BS

psO_jki_seq2_seq4_Ki_T2_RH <- subset_samples(psO_jki_seq2_seq4_Ki_T2_filt, Microhabitat == "RH")
psO_jki_seq2_seq4_Ki_T2_RH

psO_jki_seq2_seq4_Ki_T2_RP <- subset_samples(psO_jki_seq2_seq4_Ki_T2_filt, Microhabitat == "RP")
psO_jki_seq2_seq4_Ki_T2_RP

psO_jki_seq2_seq4_Ki_T3_BS <- subset_samples(psO_jki_seq2_seq4_Ki_T3_filt, Microhabitat == "BS")
psO_jki_seq2_seq4_Ki_T3_BS

psO_jki_seq2_seq4_Ki_T3_RH <- subset_samples(psO_jki_seq2_seq4_Ki_T3_filt, Microhabitat == "RH")
psO_jki_seq2_seq4_Ki_T3_RH

psO_jki_seq2_seq4_Ki_T3_RP <- subset_samples(psO_jki_seq2_seq4_Ki_T3_filt, Microhabitat == "RP")
psO_jki_seq2_seq4_Ki_T3_RP


#### Separated by location, stage T2-T3 and microhabitat 
psO_jki_seq2_seq4_Go_T2_T3_BS <- subset_samples(psO_jki_seq2_seq4_Go_filt, Microhabitat == "BS" & Stage %in%c("T2", "T3"))
psO_jki_seq2_seq4_Go_T2_T3_BS

psO_jki_seq2_seq4_Go_T2_T3_RH <- subset_samples(psO_jki_seq2_seq4_Go_filt, Microhabitat == "RH" & Stage %in%c("T2", "T3"))
psO_jki_seq2_seq4_Go_T2_T3_RH

psO_jki_seq2_seq4_Go_T2_T3_RP <- subset_samples(psO_jki_seq2_seq4_Go_filt, Microhabitat == "RP" & Stage %in%c("T2", "T3"))
psO_jki_seq2_seq4_Go_T2_T3_RP

psO_jki_seq2_seq4_Ki_T2_T3_BS <- subset_samples(psO_jki_seq2_seq4_Ki_filt, Microhabitat == "BS" & Stage %in%c("T2", "T3"))
psO_jki_seq2_seq4_Ki_T2_T3_BS

psO_jki_seq2_seq4_Ki_T2_T3_RH <- subset_samples(psO_jki_seq2_seq4_Ki_filt, Microhabitat == "RH" & Stage %in%c("T2", "T3"))
psO_jki_seq2_seq4_Ki_T2_T3_RH

psO_jki_seq2_seq4_Ki_T2_T3_RP <- subset_samples(psO_jki_seq2_seq4_Ki_filt, Microhabitat == "RP" & Stage %in%c("T2", "T3"))
psO_jki_seq2_seq4_Ki_T2_T3_RP


###Filtering table removing ASV that do not have counts in any sample on the group subset
##Using prune_taxa()

psO_jki_seq2_seq4_Go_T2_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T2_BS) > 0, psO_jki_seq2_seq4_Go_T2_BS)
psO_jki_seq2_seq4_Go_T2_BS_filt

psO_jki_seq2_seq4_Go_T2_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T2_RH) > 0, psO_jki_seq2_seq4_Go_T2_RH)
psO_jki_seq2_seq4_Go_T2_RH_filt

psO_jki_seq2_seq4_Go_T2_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T2_RP) > 0, psO_jki_seq2_seq4_Go_T2_RP)
psO_jki_seq2_seq4_Go_T2_RP_filt

psO_jki_seq2_seq4_Go_T3_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T3_BS) > 0, psO_jki_seq2_seq4_Go_T3_BS)
psO_jki_seq2_seq4_Go_T3_BS_filt

psO_jki_seq2_seq4_Go_T3_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T3_RH) > 0, psO_jki_seq2_seq4_Go_T3_RH)
psO_jki_seq2_seq4_Go_T3_RH_filt

psO_jki_seq2_seq4_Go_T3_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T3_RP) > 0, psO_jki_seq2_seq4_Go_T3_RP)
psO_jki_seq2_seq4_Go_T3_RP_filt


psO_jki_seq2_seq4_Ki_T2_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T2_BS) > 0, psO_jki_seq2_seq4_Ki_T2_BS)
psO_jki_seq2_seq4_Ki_T2_BS_filt

psO_jki_seq2_seq4_Ki_T2_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T2_RH) > 0, psO_jki_seq2_seq4_Ki_T2_RH)
psO_jki_seq2_seq4_Ki_T2_RH_filt

psO_jki_seq2_seq4_Ki_T2_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T2_RP) > 0, psO_jki_seq2_seq4_Ki_T2_RP)
psO_jki_seq2_seq4_Ki_T2_RP_filt

psO_jki_seq2_seq4_Ki_T3_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T3_BS) > 0, psO_jki_seq2_seq4_Ki_T3_BS)
psO_jki_seq2_seq4_Ki_T3_BS_filt

psO_jki_seq2_seq4_Ki_T3_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T3_RH) > 0, psO_jki_seq2_seq4_Ki_T3_RH)
psO_jki_seq2_seq4_Ki_T3_RH_filt

psO_jki_seq2_seq4_Ki_T3_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T3_RP) > 0, psO_jki_seq2_seq4_Ki_T3_RP)
psO_jki_seq2_seq4_Ki_T3_RP_filt


psO_jki_seq2_seq4_Go_T2_T3_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T2_T3_BS) > 0, psO_jki_seq2_seq4_Go_T2_T3_BS)
psO_jki_seq2_seq4_Go_T2_T3_BS_filt

psO_jki_seq2_seq4_Go_T2_T3_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T2_T3_RH) > 0, psO_jki_seq2_seq4_Go_T2_T3_RH)
psO_jki_seq2_seq4_Go_T2_T3_RH_filt

psO_jki_seq2_seq4_Go_T2_T3_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_T2_T3_RP) > 0, psO_jki_seq2_seq4_Go_T2_T3_RP)
psO_jki_seq2_seq4_Go_T2_T3_RP_filt

psO_jki_seq2_seq4_Ki_T2_T3_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T2_T3_BS) > 0, psO_jki_seq2_seq4_Ki_T2_T3_BS)
psO_jki_seq2_seq4_Ki_T2_T3_BS_filt

psO_jki_seq2_seq4_Ki_T2_T3_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T2_T3_RH) > 0, psO_jki_seq2_seq4_Ki_T2_T3_RH)
psO_jki_seq2_seq4_Ki_T2_T3_RH_filt

psO_jki_seq2_seq4_Ki_T2_T3_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_T2_T3_RP) > 0, psO_jki_seq2_seq4_Ki_T2_T3_RP)
psO_jki_seq2_seq4_Ki_T2_T3_RP_filt

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
### Different locations
MDS_bray_psO_jki_seq2_seq4_Go_T2_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T2_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_T2_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T2_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_T2_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T2_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Go_T3_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T3_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_T3_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T3_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_T3_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T3_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T2_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T2_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T2_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T3_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T3_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T3_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T2_T3_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T2_T3_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_T2_T3_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T2_T3_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T2_T3_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_T2_T3_RP_filt, "MDS","bray", autotransform=TRUE)

head(MDS_bray_psO_jki_seq2_seq4_Go_T2_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_T2_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_T2_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Go_T3_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_T3_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_T3_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr)


#Create a MDS plot 
#### ASVs
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T2_BS_filt, MDS_bray_psO_jki_seq2_seq4_Go_T2_BS_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T2_RH_filt, MDS_bray_psO_jki_seq2_seq4_Go_T2_RH_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T2_RP_filt, MDS_bray_psO_jki_seq2_seq4_Go_T2_RP_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T3_BS_filt, MDS_bray_psO_jki_seq2_seq4_Go_T3_BS_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T3_RH_filt, MDS_bray_psO_jki_seq2_seq4_Go_T3_RH_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T3_RP_filt, MDS_bray_psO_jki_seq2_seq4_Go_T3_RP_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T3_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T2_BS_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T2_RH_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T2_RP_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T3_BS_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T3_RH_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T3_RP_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T2_T3_BS_filt, MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T2_T3_RH_filt, MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_T2_T3_RP_filt, MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T2_T3_BS_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T2_T3_RH_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_T2_T3_RP_filt, MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr, type="sample",color="Rotation") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")




###########################################
#arrange multiple ggplots in one single figure 
#alpha_plots_rotation_depth <- ggarrange(plot_MDS_bray_psO_jki_seq1_W1_filt_sqr, plot_MDS_bray_psO_jki_seq1_W2_filt_sqr, plot_MDS_bray_psO_jki_seq1_WM_filt_sqr, ncol =3, nrow =1, common.legend = TRUE, legend="right")
#ggsave("alpha_plots_rotation_depth.png", path = "~/Documents/R_analysis_2/jki_seq1/output/Ordination/", width = 70, height = 20, units = "cm",dpi = 300)
###########################################


####################################### Statistics
set.seed(2022)
##Multivariate Anova for ordinate files

#Transform abundances to square root
psO_jki_seq2_seq4_Go_T2_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T2_BS_filt, "hellinger")
psO_jki_seq2_seq4_Go_T2_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T2_RH_filt, "hellinger")
psO_jki_seq2_seq4_Go_T2_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T2_RP_filt, "hellinger")

psO_jki_seq2_seq4_Go_T3_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T3_BS_filt, "hellinger")
psO_jki_seq2_seq4_Go_T3_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T3_RH_filt, "hellinger")
psO_jki_seq2_seq4_Go_T3_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T3_RP_filt, "hellinger")

psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T2_BS_filt, "hellinger")
psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T2_RH_filt, "hellinger")
psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T2_RP_filt, "hellinger")

psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T3_BS_filt, "hellinger")
psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T3_RH_filt, "hellinger")
psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T3_RP_filt, "hellinger")

psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T2_T3_BS_filt, "hellinger")
psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T2_T3_RH_filt, "hellinger")
psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_T2_T3_RP_filt, "hellinger")

psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T2_T3_BS_filt, "hellinger")
psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T2_T3_RH_filt, "hellinger")
psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_T2_T3_RP_filt, "hellinger")



#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq2_seq4_Go_T2_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T2_BS_filt_sqr)
psO_jki_seq2_seq4_Go_T2_BS_filt_meta <- meta(psO_jki_seq2_seq4_Go_T2_BS_filt)

psO_jki_seq2_seq4_Go_T2_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T2_RH_filt_sqr)
psO_jki_seq2_seq4_Go_T2_RH_filt_meta <- meta(psO_jki_seq2_seq4_Go_T2_RH_filt)

psO_jki_seq2_seq4_Go_T2_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T2_RP_filt_sqr)
psO_jki_seq2_seq4_Go_T2_RP_filt_meta <- meta(psO_jki_seq2_seq4_Go_T2_RP_filt)

psO_jki_seq2_seq4_Go_T3_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T3_BS_filt_sqr)
psO_jki_seq2_seq4_Go_T3_BS_filt_meta <- meta(psO_jki_seq2_seq4_Go_T3_BS_filt)

psO_jki_seq2_seq4_Go_T3_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T3_RH_filt_sqr)
psO_jki_seq2_seq4_Go_T3_RH_filt_meta <- meta(psO_jki_seq2_seq4_Go_T3_RH_filt)

psO_jki_seq2_seq4_Go_T3_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T3_RP_filt_sqr)
psO_jki_seq2_seq4_Go_T3_RP_filt_meta <- meta(psO_jki_seq2_seq4_Go_T3_RP_filt)


psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr)
psO_jki_seq2_seq4_Ki_T2_BS_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T2_BS_filt)

psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr)
psO_jki_seq2_seq4_Ki_T2_RH_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T2_RH_filt)

psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr)
psO_jki_seq2_seq4_Ki_T2_RP_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T2_RP_filt)

psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T3_BS_filt_sqr)
psO_jki_seq2_seq4_Ki_T3_BS_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T3_BS_filt)

psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T3_RH_filt_sqr)
psO_jki_seq2_seq4_Ki_T3_RH_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T3_RH_filt)

psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T3_RP_filt_sqr)
psO_jki_seq2_seq4_Ki_T3_RP_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T3_RP_filt)


psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr)
psO_jki_seq2_seq4_Go_T2_T3_BS_filt_meta <- meta(psO_jki_seq2_seq4_Go_T2_T3_BS_filt)

psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr)
psO_jki_seq2_seq4_Go_T2_T3_RH_filt_meta <- meta(psO_jki_seq2_seq4_Go_T2_T3_RH_filt)

psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr)
psO_jki_seq2_seq4_Go_T2_T3_RP_filt_meta <- meta(psO_jki_seq2_seq4_Go_T2_T3_RP_filt)

psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr)
psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T2_T3_BS_filt)

psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr)
psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T2_T3_RH_filt)

psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr)
psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_meta <- meta(psO_jki_seq2_seq4_Ki_T2_T3_RP_filt)


#######Permanova using adonis function from vegan and print p value
#ASVs
set.seed(2022)

Permanova_terms_psO_jki_seq2_seq4_Go_T2_BS_year <- adonis2(t(psO_jki_seq2_seq4_Go_T2_BS_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Go_T2_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_T2_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_T2_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_T2_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_T2_RH_year <- adonis2(t(psO_jki_seq2_seq4_Go_T2_RH_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Go_T2_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_T2_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_T2_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_T2_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_T2_RP_year <- adonis2(t(psO_jki_seq2_seq4_Go_T2_RP_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Go_T2_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_T2_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_T2_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_T2_RP_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_T2_BS_year <- adonis2(t(psO_jki_seq2_seq4_Ki_T2_BS_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Ki_T2_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_T2_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_T2_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_T2_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RH_year <- adonis2(t(psO_jki_seq2_seq4_Ki_T2_RH_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Ki_T2_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RP_year <- adonis2(t(psO_jki_seq2_seq4_Ki_T2_RP_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Ki_T2_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_T2_RP_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_BS_year <- adonis2(t(psO_jki_seq2_seq4_Go_T2_T3_BS_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Go_T2_T3_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RH_year <- adonis2(t(psO_jki_seq2_seq4_Go_T2_T3_RH_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Go_T2_T3_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RP_year <- adonis2(t(psO_jki_seq2_seq4_Go_T2_T3_RP_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Go_T2_T3_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_T2_T3_RP_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_BS_year <- adonis2(t(psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Ki_T2_T3_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RH_year <- adonis2(t(psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Ki_T2_T3_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RP_year <- adonis2(t(psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_sqr_abundances) ~Year*Rotation, data = psO_jki_seq2_seq4_Ki_T2_T3_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_T2_T3_RP_year.csv")



## The end : )



