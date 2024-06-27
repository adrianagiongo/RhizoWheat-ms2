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

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
### Different locations
MDS_bray_psO_jki_seq2_seq4_Go_2020_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_2020_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_2020_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_2020_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_2020_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_2020_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Go_2021_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_2021_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_2021_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_2021_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_2021_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Go_2021_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_2020_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_2020_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_2020_RP_filt, "MDS","bray", autotransform=TRUE)

MDS_bray_psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_2021_BS_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_2021_RH_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr<-ordinate(psO_jki_seq2_seq4_Ki_2021_RP_filt, "MDS","bray", autotransform=TRUE)


head(MDS_bray_psO_jki_seq2_seq4_Go_2020_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_2020_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_2020_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Go_2021_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_2021_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_2021_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr)

head(MDS_bray_psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr)


#Create a MDS plot 
#### ASVs  Go
plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_2020_BS_filt, MDS_bray_psO_jki_seq2_seq4_Go_2020_BS_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_2020_RH_filt, MDS_bray_psO_jki_seq2_seq4_Go_2020_RH_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_2020_RP_filt, MDS_bray_psO_jki_seq2_seq4_Go_2020_RP_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_2020_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")




plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_2021_BS_filt, MDS_bray_psO_jki_seq2_seq4_Go_2021_BS_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_2021_RH_filt, MDS_bray_psO_jki_seq2_seq4_Go_2021_RH_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_2021_RP_filt, MDS_bray_psO_jki_seq2_seq4_Go_2021_RP_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_2021_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


#### ASVs  Ki
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_2020_BS_filt, MDS_bray_psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_2020_RH_filt, MDS_bray_psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_2020_RP_filt, MDS_bray_psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")




plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_2021_BS_filt, MDS_bray_psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_2021_RH_filt, MDS_bray_psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_2021_RP_filt, MDS_bray_psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr, type="sample",color="Rotation", shape = "Stage") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c("#cfceb7", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


####################################### Statistics
set.seed(2022)
##Multivariate Anova for ordinate files

#Transform abundances to square root
psO_jki_seq2_seq4_Go_2020_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_2020_BS_filt, "hellinger")
psO_jki_seq2_seq4_Go_2020_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_2020_RH_filt, "hellinger")
psO_jki_seq2_seq4_Go_2020_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_2020_RP_filt, "hellinger")

psO_jki_seq2_seq4_Go_2021_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_2021_BS_filt, "hellinger")
psO_jki_seq2_seq4_Go_2021_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_2021_RH_filt, "hellinger")
psO_jki_seq2_seq4_Go_2021_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_2021_RP_filt, "hellinger")

psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_2020_BS_filt, "hellinger")
psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_2020_RH_filt, "hellinger")
psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_2020_RP_filt, "hellinger")

psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_2021_BS_filt, "hellinger")
psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_2021_RH_filt, "hellinger")
psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_2021_RP_filt, "hellinger")


#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq2_seq4_Go_2020_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_2020_BS_filt_sqr)
psO_jki_seq2_seq4_Go_2020_BS_filt_meta <- meta(psO_jki_seq2_seq4_Go_2020_BS_filt)

psO_jki_seq2_seq4_Go_2020_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_2020_RH_filt_sqr)
psO_jki_seq2_seq4_Go_2020_RH_filt_meta <- meta(psO_jki_seq2_seq4_Go_2020_RH_filt)

psO_jki_seq2_seq4_Go_2020_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_2020_RP_filt_sqr)
psO_jki_seq2_seq4_Go_2020_RP_filt_meta <- meta(psO_jki_seq2_seq4_Go_2020_RP_filt)


psO_jki_seq2_seq4_Go_2021_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_2021_BS_filt_sqr)
psO_jki_seq2_seq4_Go_2021_BS_filt_meta <- meta(psO_jki_seq2_seq4_Go_2021_BS_filt)

psO_jki_seq2_seq4_Go_2021_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_2021_RH_filt_sqr)
psO_jki_seq2_seq4_Go_2021_RH_filt_meta <- meta(psO_jki_seq2_seq4_Go_2021_RH_filt)

psO_jki_seq2_seq4_Go_2021_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_2021_RP_filt_sqr)
psO_jki_seq2_seq4_Go_2021_RP_filt_meta <- meta(psO_jki_seq2_seq4_Go_2021_RP_filt)


psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr)
psO_jki_seq2_seq4_Ki_2020_BS_filt_meta <- meta(psO_jki_seq2_seq4_Ki_2020_BS_filt)

psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr)
psO_jki_seq2_seq4_Ki_2020_RH_filt_meta <- meta(psO_jki_seq2_seq4_Ki_2020_RH_filt)

psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr)
psO_jki_seq2_seq4_Ki_2020_RP_filt_meta <- meta(psO_jki_seq2_seq4_Ki_2020_RP_filt)


psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr)
psO_jki_seq2_seq4_Ki_2021_BS_filt_meta <- meta(psO_jki_seq2_seq4_Ki_2021_BS_filt)

psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr)
psO_jki_seq2_seq4_Ki_2021_RH_filt_meta <- meta(psO_jki_seq2_seq4_Ki_2021_RH_filt)

psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr)
psO_jki_seq2_seq4_Ki_2021_RP_filt_meta <- meta(psO_jki_seq2_seq4_Ki_2021_RP_filt)

#######Permanova using adonis function from vegan and print p value
#ASVs
set.seed(2022)

Permanova_terms_psO_jki_seq2_seq4_Go_2020_BS_year <- adonis2(t(psO_jki_seq2_seq4_Go_2020_BS_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Go_2020_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_2020_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_2020_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_2020_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_2020_RH_year <- adonis2(t(psO_jki_seq2_seq4_Go_2020_RH_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Go_2020_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_2020_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_2020_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_2020_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_2020_RP_year <- adonis2(t(psO_jki_seq2_seq4_Go_2020_RP_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Go_2020_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_2020_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_2020_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_2020_RP_year.csv")


Permanova_terms_psO_jki_seq2_seq4_Go_2021_BS_year <- adonis2(t(psO_jki_seq2_seq4_Go_2021_BS_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Go_2021_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_2021_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_2021_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_2021_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_2021_RH_year <- adonis2(t(psO_jki_seq2_seq4_Go_2021_RH_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Go_2021_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_2021_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_2021_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_2021_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go_2021_RP_year <- adonis2(t(psO_jki_seq2_seq4_Go_2021_RP_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Go_2021_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go_2021_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go_2021_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go_2021_RP_year.csv")


Permanova_terms_psO_jki_seq2_seq4_Ki_2020_BS_year <- adonis2(t(psO_jki_seq2_seq4_Ki_2020_BS_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Ki_2020_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_2020_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_2020_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_2020_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RH_year <- adonis2(t(psO_jki_seq2_seq4_Ki_2020_RH_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Ki_2020_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RP_year <- adonis2(t(psO_jki_seq2_seq4_Ki_2020_RP_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Ki_2020_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_2020_RP_year.csv")


Permanova_terms_psO_jki_seq2_seq4_Ki_2021_BS_year <- adonis2(t(psO_jki_seq2_seq4_Ki_2021_BS_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Ki_2021_BS_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_2021_BS_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_2021_BS_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_2021_BS_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RH_year <- adonis2(t(psO_jki_seq2_seq4_Ki_2021_RH_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Ki_2021_RH_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RH_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RH_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RH_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RP_year <- adonis2(t(psO_jki_seq2_seq4_Ki_2021_RP_filt_sqr_abundances) ~Rotation*Stage, data = psO_jki_seq2_seq4_Ki_2021_RP_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RP_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RP_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki_2021_RP_year.csv")


## The end : )



