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
head(MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr)
head(MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr)


#Create a MDS plot 
#### ASVs  Go
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_microhabitat<-plot_ordination(psO_jki_seq2_seq4_rarefied_Go_2020_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr, type="sample", color="Microhabitat", shape="Rotation") + 
  geom_point(size=4) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_shape_manual(values=c(16, 17, 15))+
  scale_color_manual(values = c("#b05644", "#d9b967", "#57896a"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_microhabitat

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_microhabitat.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_microhabitat<-plot_ordination(psO_jki_seq2_seq4_rarefied_Go_2021_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr, type="sample",color="Microhabitat", shape="Rotation") + 
  geom_point(size=4) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_shape_manual(values=c(16, 17, 15))+
  scale_color_manual(values = c("#b05644", "#d9b967", "#57896a"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_microhabitat

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_microhabitat.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")

#### ASVs  Ki
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_microhabitat<-plot_ordination(psO_jki_seq2_seq4_rarefied_Ki_2020_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr, type="sample",color="Microhabitat", shape="Rotation") + 
  geom_point(size=4) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_shape_manual(values=c(16, 18))+
  scale_color_manual(values = c("#b05644", "#d9b967", "#57896a"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_microhabitat

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_microhabitat.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")

plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_microhabitat<-plot_ordination(psO_jki_seq2_seq4_rarefied_Ki_2021_filt, MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr, type="sample",color="Microhabitat", shape="Rotation") + 
  geom_point(size=4) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_shape_manual(values=c(16, 18))+
  scale_color_manual(values = c("#b05644", "#d9b967", "#57896a"))
plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_microhabitat

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_microhabitat.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


####################################### Statistics
set.seed(2022)
##Multivariate Anova for ordinate files

#Transform abundances to square root
psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_rarefied_Go_2020_filt, "hellinger")
psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_rarefied_Go_2021_filt, "hellinger")
psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_rarefied_Ki_2020_filt, "hellinger")
psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr <- microbiome::transform(psO_jki_seq2_seq4_rarefied_Ki_2021_filt, "hellinger")


#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr)
psO_jki_seq2_seq4_rarefied_Go_2020_filt_meta <- meta(psO_jki_seq2_seq4_rarefied_Go_2020_filt)

psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr)
psO_jki_seq2_seq4_rarefied_Go_2021_filt_meta <- meta(psO_jki_seq2_seq4_rarefied_Go_2021_filt)

psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr)
psO_jki_seq2_seq4_rarefied_Ki_2020_filt_meta <- meta(psO_jki_seq2_seq4_rarefied_Ki_2020_filt)

psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_abundances <- abundances(psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr)
psO_jki_seq2_seq4_rarefied_Ki_2021_filt_meta <- meta(psO_jki_seq2_seq4_rarefied_Ki_2021_filt)

#######Permanova using adonis function from vegan and print p value
#ASVs
set.seed(2022)

Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2020_year <- adonis2(t(psO_jki_seq2_seq4_rarefied_Go_2020_filt_sqr_abundances) ~Microhabitat*Rotation*Stage, data = psO_jki_seq2_seq4_rarefied_Go_2020_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2020_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2020_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2020_year.csv")

Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2021_year <- adonis2(t(psO_jki_seq2_seq4_rarefied_Go_2021_filt_sqr_abundances) ~Microhabitat*Rotation*Stage, data = psO_jki_seq2_seq4_rarefied_Go_2021_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2021_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2021_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_rarefied_Go_2021_year.csv")

Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2020_year <- adonis2(t(psO_jki_seq2_seq4_rarefied_Ki_2020_filt_sqr_abundances) ~Microhabitat*Rotation*Stage, data = psO_jki_seq2_seq4_rarefied_Ki_2020_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2020_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2020_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2020_year.csv")

Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2021_year <- adonis2(t(psO_jki_seq2_seq4_rarefied_Ki_2021_filt_sqr_abundances) ~Microhabitat*Rotation*Stage, data = psO_jki_seq2_seq4_rarefied_Ki_2021_filt_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2021_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2021_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_rarefied_Ki_2021_year.csv")

## The end : )






