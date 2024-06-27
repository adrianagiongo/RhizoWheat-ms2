##Creating MDS plot for selected dataset
#load package
library("phyloseq")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")
library("ggpubr")

###########
#Go = "#1d5b65"   Go 2020 = "#216974"  Go 2021 = "#6baea5"
#Ki = "#A34828"   Ki 2020 = "#D1711F"  Ki 2021 = "#ebaf7a"

################################
set.seed(2022)

# Set color palette
Location_colors <- c("Go" = "#1d5b65", "Ki" = "#A34828")


###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors
#(1)Location = Go or Ki
#(2)Year = 2020 or 2021 
#(3)Microhabitat = W1, W2, WM (Go) or W1, W3 (Ki)

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
### Different locations
MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr<-ordinate(psO_jki_seq2_seq4_rarefied_Go_2020_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr<-ordinate(psO_jki_seq2_seq4_rarefied_Go_2021_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr<-ordinate(psO_jki_seq2_seq4_rarefied_Ki_2020_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr<-ordinate(psO_jki_seq2_seq4_rarefied_Ki_2021_filt, "MDS","bray", autotransform=TRUE)

head(MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr)


#Create a MDS plot 
#### ASVs  Go
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_rotation<-plot_ordination(psO_jki_seq2_seq4_rarefied_Go_2020_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr, type="sample", color="Rotation", shape="Microhabitat") + 
  geom_point(size=5) +
  theme_bw() +
  theme(axis.text = element_text(size=22)) +
  theme(axis.title = element_text(size=22)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=22)) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_x_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_color_manual(values = c("#e7e5cc", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 14, height = 10, units = "cm", dpi = 300, device = "svg")


plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_rotation<-plot_ordination(psO_jki_seq2_seq4_rarefied_Go_2021_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr, type="sample",color="Rotation", shape="Microhabitat") + 
  geom_point(size=5) +
  theme_bw() +
  theme(axis.text = element_text(size=22)) +
  theme(axis.title = element_text(size=22)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=22)) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_x_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_color_manual(values = c("#e7e5cc", "#9cc184", "#1e3d14"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 14, height = 10, units = "cm", dpi = 300, device = "svg")

#### ASVs  Ki
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_rotation<-plot_ordination(psO_jki_seq2_seq4_rarefied_Ki_2020_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr, type="sample",color="Rotation", shape="Microhabitat") + 
  geom_point(size=5) +
  theme_bw() +
  theme(axis.text = element_text(size=22)) +
  theme(axis.title = element_text(size=22)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=22)) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_x_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_color_manual(values = c("#e7e5cc", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 14, height = 10, units = "cm", dpi = 300, device = "svg")

plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_rotation<-plot_ordination(psO_jki_seq2_seq4_rarefied_Ki_2021_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr, type="sample",color="Rotation", shape="Microhabitat") + 
  geom_point(size=5) +
  theme_bw() +
  theme(axis.text = element_text(size=22)) +
  theme(axis.title = element_text(size=22)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=22)) +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_x_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, by = 0.3), label = c("-0.3", "0", "0.3")) +
  scale_color_manual(values = c("#e7e5cc", "#447243"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_rotation

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_rotation.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 14, height = 10, units = "cm", dpi = 300, device = "svg")



###########################################
#arrange multiple ggplots in one single figure ALL
beta_rotations_plots_all <- ggarrange(plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_rotation,
                                      plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_rotation,
                                      plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_rotation,
                                      plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_rotation, 
                             ncol =2, nrow =2, legend="right", widths = c(1, 1), common.legend = TRUE, label.y = 1)
beta_rotations_plots_all

ggsave("beta_rotations_plots_all.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 25, height = 20, units = "cm",dpi = 300, device = "svg")
###########################################

# ####################################### Statistics   ---> done already in the Ordination....microhabitat

## The end : )



