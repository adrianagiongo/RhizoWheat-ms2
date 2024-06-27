###########
##### Correspondent to Figure 1B and Suppl Figure S1B
###########

##Creating MDS plot for selected dataset
#load package
library("phyloseq")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")
library("ggpubr")

################################   All samples    ---------> ASV  (Figure 1B)
#set seed
set.seed(2022)

# Set color palette
Site_colors <- c("Go" = "#1d5b65", "Ki" = "#A34828")
#Go = "#1d5b65"   Go 2020 = "#216974"  Go 2021 = "#6baea5"
#Ki = "#A34828"   Ki 2020 = "#D1711F"  Ki 2021 = "#ebaf7a"

### All samples ----> MDS
MDS_bray_psO_jki_seq2_seq4_sqr<-ordinate(psO_jki_seq2_seq4, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Go_sqr<-ordinate(psO_jki_seq2_seq4_Go_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq2_seq4_Ki_sqr<-ordinate(psO_jki_seq2_seq4_Ki_filt, "MDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_jki_seq2_seq4_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Go_sqr)
head(MDS_bray_psO_jki_seq2_seq4_Ki_sqr)

#Create a MDS plot 
plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site<-plot_ordination(psO_jki_seq2_seq4, MDS_bray_psO_jki_seq2_seq4_sqr, type="sample",color="Site", shape = "Microhabitat") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_colors <- c("#1d5b65", "#A34828")))
plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site 

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site_year_microh<-plot_ordination(psO_jki_seq2_seq4, MDS_bray_psO_jki_seq2_seq4_sqr, type="sample",color="Site_year", shape = "Microhabitat") + 
  geom_point(size=4) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_text(size=18)) +
  theme(legend.text = element_text(size=18)) +
  scale_color_manual(values = c(Site_year_colors <- c("#377881", "#6baea5", "#D1711F", "#ebaf7a")))
plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site_year_microh 

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site_year_microh.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


###### To show Stage 
plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site_year_stage<-plot_ordination(psO_jki_seq2_seq4, MDS_bray_psO_jki_seq2_seq4_sqr, type="sample",color="Site_year", shape = "Stage") + 
  geom_point(size=4) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_text(size=18)) +
  theme(legend.text = element_text(size=18)) +
  scale_color_manual(values = c(Site_year_colors <- c("#377881", "#6baea5", "#D1711F", "#ebaf7a")))
plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site_year_stage 

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_sqr_Site_year_stage.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")




plot_MDS_bray_psO_jki_seq2_seq4_Go_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_filt, MDS_bray_psO_jki_seq2_seq4_Go_sqr, type="sample",color="Site_year") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_year_colors <- c("#216974", "#6baea5")))
plot_MDS_bray_psO_jki_seq2_seq4_Go_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Go_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_MDS_bray_psO_jki_seq2_seq4_Ki_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_filt, MDS_bray_psO_jki_seq2_seq4_Ki_sqr, type="sample",color="Site_year") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_year_colors <- c("#D1711F", "#ebaf7a")))
plot_MDS_bray_psO_jki_seq2_seq4_Ki_sqr

ggsave("plot_MDS_bray_psO_jki_seq2_seq4_Ki_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



### All samples   -----> NMDS
NMDS_bray_psO_jki_seq2_seq4_sqr<-ordinate(psO_jki_seq2_seq4, "NMDS","bray", autotransform=TRUE)
NMDS_bray_psO_jki_seq2_seq4_Go_sqr<-ordinate(psO_jki_seq2_seq4_Go_filt, "NMDS","bray", autotransform=TRUE)
NMDS_bray_psO_jki_seq2_seq4_Ki_sqr<-ordinate(psO_jki_seq2_seq4_Ki_filt, "NMDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(NMDS_bray_psO_jki_seq2_seq4_sqr)
head(NMDS_bray_psO_jki_seq2_seq4_Go_sqr)
head(NMDS_bray_psO_jki_seq2_seq4_Ki_sqr)

#Create a NMDS plot 
plot_NMDS_bray_psO_jki_seq2_seq4_sqr_Site<-plot_ordination(psO_jki_seq2_seq4, NMDS_bray_psO_jki_seq2_seq4_sqr, type="sample",color="Site") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_colors <- c("#1d5b65", "#A34828")))
plot_NMDS_bray_psO_jki_seq2_seq4_sqr_Site 

ggsave("plot_NMDS_bray_psO_jki_seq2_seq4_sqr_Site.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")


plot_NMDS_bray_psO_jki_seq2_seq4_sqr_Site_year<-plot_ordination(psO_jki_seq2_seq4, NMDS_bray_psO_jki_seq2_seq4_sqr, type="sample",color="Site_year") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_year_colors <- c("#216974", "#6baea5", "#D1711F", "#ebaf7a")))
plot_NMDS_bray_psO_jki_seq2_seq4_sqr_Site_year 

ggsave("plot_NMDS_bray_psO_jki_seq2_seq4_sqr_Site_year.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Ordination_jki_seq2_seq4/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")



plot_NMDS_bray_psO_jki_seq2_seq4_Go_sqr<-plot_ordination(psO_jki_seq2_seq4_Go_filt, NMDS_bray_psO_jki_seq2_seq4_Go_sqr, type="sample",color="Site_year") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_year_colors <- c("#216974", "#6baea5")))
plot_NMDS_bray_psO_jki_seq2_seq4_Go_sqr 

ggsave("plot_NMDS_bray_psO_jki_seq2_seq4_Go_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4_Go/output_jki_seq2_seq4_Go/Ordination_jki_seq2_seq4_Go/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")

plot_NMDS_bray_psO_jki_seq2_seq4_Ki_sqr<-plot_ordination(psO_jki_seq2_seq4_Ki_filt, NMDS_bray_psO_jki_seq2_seq4_Ki_sqr, type="sample",color="Site_year") + 
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.text = element_text(size=20)) +
  scale_color_manual(values = c(Site_year_colors <- c("#D1711F", "#ebaf7a")))
plot_NMDS_bray_psO_jki_seq2_seq4_Ki_sqr 

ggsave("plot_NMDS_bray_psO_jki_seq2_seq4_Ki_sqr.tiff", path = "~/Documents/R_analysis/jki_seq2_seq4_Ki/output_jki_seq2_seq4_Ki/Ordination_jki_seq2_seq4_Ki/", width = 15, height = 10, units = "cm", dpi = 300, device = "tiff")




################################   All samples    ---------> ASV    ----------> Statistics
set.seed(2022)

##Multivariate Anova for ordinate files
#Transform abundances to square root
psO_jki_seq2_seq4_sqr <- microbiome::transform(psO_jki_seq2_seq4, "hellinger")
psO_jki_seq2_seq4_Go_sqr <- microbiome::transform(psO_jki_seq2_seq4_Go_filt, "hellinger")
psO_jki_seq2_seq4_Ki_sqr <- microbiome::transform(psO_jki_seq2_seq4_Ki_filt, "hellinger")

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq2_seq4_sqr_abundances <- abundances(psO_jki_seq2_seq4_sqr)
psO_jki_seq2_seq4_meta <- meta(psO_jki_seq2_seq4)

psO_jki_seq2_seq4_Go_sqr_abundances <- abundances(psO_jki_seq2_seq4_Go_sqr)
psO_jki_seq2_seq4_Go_meta <- meta(psO_jki_seq2_seq4_Go_filt)

psO_jki_seq2_seq4_Ki_sqr_abundances <- abundances(psO_jki_seq2_seq4_Ki_sqr)
psO_jki_seq2_seq4_Ki_meta <- meta(psO_jki_seq2_seq4_Ki_filt)

#######Permanova using adonis function from vegan and print p value
#ASVs
#by terms (it is default)
Permanova_terms_psO_jki_seq2_seq4_Site_year <- adonis2(t(psO_jki_seq2_seq4_sqr_abundances) ~Site*Year*Microhabitat*Stage*Rotation, data = psO_jki_seq2_seq4_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Site_year
write.csv(Permanova_terms_psO_jki_seq2_seq4_Site_year, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Site_year.csv")

Permanova_terms_psO_jki_seq2_seq4_Go <- adonis2(t(psO_jki_seq2_seq4_Go_sqr_abundances) ~Year*Microhabitat*Stage*Rotation, data = psO_jki_seq2_seq4_Go_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Go
write.csv(Permanova_terms_psO_jki_seq2_seq4_Go, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Go.csv")

Permanova_terms_psO_jki_seq2_seq4_Ki <- adonis2(t(psO_jki_seq2_seq4_Ki_sqr_abundances) ~Year*Microhabitat*Stage*Rotation, data = psO_jki_seq2_seq4_Ki_meta, permutations = 10000, method = "bray", by= "terms")
Permanova_terms_psO_jki_seq2_seq4_Ki
write.csv(Permanova_terms_psO_jki_seq2_seq4_Ki, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/Permanova_terms_psO_jki_seq2_seq4_Ki.csv")





#### ANOSIM    -------> MICROHABITATS
## testing of significance for the bray-curtis dissimilarity using ANOSIM
#calculate distance values between samples
dist_psO_jki_seq2_seq4 = phyloseq::distance(psO_jki_seq2_seq4, method = "bray")

dist_psO_jki_seq2_seq4_Go = phyloseq::distance(psO_jki_seq2_seq4_Go_filt, method = "bray")
dist_psO_jki_seq2_seq4_Ki = phyloseq::distance(psO_jki_seq2_seq4_Ki_filt, method = "bray")

dist_psO_jki_seq2_seq4_Go_2020 = phyloseq::distance(psO_jki_seq2_seq4_Go_2020_filt, method = "bray")
dist_psO_jki_seq2_seq4_Go_2021 = phyloseq::distance(psO_jki_seq2_seq4_Go_2021_filt, method = "bray")

dist_psO_jki_seq2_seq4_Ki_2020 = phyloseq::distance(psO_jki_seq2_seq4_Ki_2020_filt, method = "bray")
dist_psO_jki_seq2_seq4_Ki_2021 = phyloseq::distance(psO_jki_seq2_seq4_Ki_2021_filt, method = "bray")

#create a dataframe of the metadata
metadata_psO_jki_seq2_seq4 <- data.frame(sample_data(psO_jki_seq2_seq4))
metadata_psO_jki_seq2_seq4_Go <- data.frame(sample_data(psO_jki_seq2_seq4_Go_filt))
metadata_psO_jki_seq2_seq4_Ki <- data.frame(sample_data(psO_jki_seq2_seq4_Ki_filt))

metadata_psO_jki_seq2_seq4_Go_2020 <- data.frame(sample_data(psO_jki_seq2_seq4_Go_2020_filt))
metadata_psO_jki_seq2_seq4_Go_2021 <- data.frame(sample_data(psO_jki_seq2_seq4_Go_2021_filt))

metadata_psO_jki_seq2_seq4_Ki_2020 <- data.frame(sample_data(psO_jki_seq2_seq4_Ki_2020_filt))
metadata_psO_jki_seq2_seq4_Ki_2021 <- data.frame(sample_data(psO_jki_seq2_seq4_Ki_2021_filt))


#run Anosim using 10000 permutations
anosim_dist_psO_jki_seq2_seq4_site_year <- anosim(dist_psO_jki_seq2_seq4, metadata_psO_jki_seq2_seq4$Site_year, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Go_Microhabitat <- anosim(dist_psO_jki_seq2_seq4_Go, metadata_psO_jki_seq2_seq4_Go$Microhabitat, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Ki_Microhabitat <- anosim(dist_psO_jki_seq2_seq4_Ki, metadata_psO_jki_seq2_seq4_Ki$Microhabitat, permutations = 10000)

anosim_dist_psO_jki_seq2_seq4_Go_2020_Microhabitat <- anosim(dist_psO_jki_seq2_seq4_Go_2020, metadata_psO_jki_seq2_seq4_Go_2020$Microhabitat, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Go_2021_Microhabitat <- anosim(dist_psO_jki_seq2_seq4_Go_2021, metadata_psO_jki_seq2_seq4_Go_2021$Microhabitat, permutations = 10000)

anosim_dist_psO_jki_seq2_seq4_Ki_2020_Microhabitat <- anosim(dist_psO_jki_seq2_seq4_Ki_2020, metadata_psO_jki_seq2_seq4_Ki_2020$Microhabitat, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Ki_2021_Microhabitat <- anosim(dist_psO_jki_seq2_seq4_Ki_2021, metadata_psO_jki_seq2_seq4_Ki_2021$Microhabitat, permutations = 10000)

#print results
print(anosim_dist_psO_jki_seq2_seq4_site_year)
print(anosim_dist_psO_jki_seq2_seq4_Go_Microhabitat)
print(anosim_dist_psO_jki_seq2_seq4_Ki_Microhabitat)

print(anosim_dist_psO_jki_seq2_seq4_Go_2020_Microhabitat)
print(anosim_dist_psO_jki_seq2_seq4_Go_2021_Microhabitat)

print(anosim_dist_psO_jki_seq2_seq4_Ki_2020_Microhabitat)
print(anosim_dist_psO_jki_seq2_seq4_Ki_2021_Microhabitat)


#### PAIRWISE ANOSIM
## testing of pairwise significance for the bray-curtis dissimilarity using ANOSIM
#Create metadata with an unique variable and run Anosim using 10000 permutations
metadata_site_year_psO_jki_seq2_seq4 <- combn(x=unique(metadata_psO_jki_seq2_seq4$Site_year), m=2)
p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- c()

metadata_microhabitat_psO_jki_seq2_seq4_Go <- combn(x=unique(metadata_psO_jki_seq2_seq4_Go$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go <- c()

metadata_microhabitat_psO_jki_seq2_seq4_Ki <- combn(x=unique(metadata_psO_jki_seq2_seq4_Ki$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki <- c()

metadata_microhabitat_psO_jki_seq2_seq4_Go_2020 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Go_2020$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020 <- c()
metadata_microhabitat_psO_jki_seq2_seq4_Go_2021 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Go_2021$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021 <- c()

metadata_microhabitat_psO_jki_seq2_seq4_Ki_2020 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Ki_2020$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020 <- c()
metadata_microhabitat_psO_jki_seq2_seq4_Ki_2021 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Ki_2021$Microhabitat), m=2)
p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021 <- c()



for (i in 1:ncol(metadata_site_year_psO_jki_seq2_seq4)){
  ps_subs_psO_jki_seq2_seq4 <- subset_samples(psO_jki_seq2_seq4, Site_year %in% metadata_site_year_psO_jki_seq2_seq4 [,i])
  metadata_subs_psO_jki_seq2_seq4 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4, method= "bray"), metadata_subs_psO_jki_seq2_seq4$Site_year)
  p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- c(p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4, pairwise_anosim_bray_dist_psO_jki_seq2_seq4$signif[1])
}

for (i in 1:ncol(metadata_microhabitat_psO_jki_seq2_seq4_Go)){
  ps_subs_psO_jki_seq2_seq4_Go <- subset_samples(psO_jki_seq2_seq4_Go, Microhabitat %in% metadata_microhabitat_psO_jki_seq2_seq4_Go [,i])
  metadata_subs_psO_jki_seq2_seq4_Go <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Go))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Go, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Go$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go$signif[1])
}

for (i in 1:ncol(metadata_microhabitat_psO_jki_seq2_seq4_Ki)){
  ps_subs_psO_jki_seq2_seq4_Ki <- subset_samples(psO_jki_seq2_seq4_Ki, Microhabitat %in% metadata_microhabitat_psO_jki_seq2_seq4_Ki [,i])
  metadata_subs_psO_jki_seq2_seq4_Ki <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Ki))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Ki, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Ki$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki$signif[1])
}


for (i in 1:ncol(metadata_microhabitat_psO_jki_seq2_seq4_Go_2020)){
  ps_subs_psO_jki_seq2_seq4_Go_2020 <- subset_samples(psO_jki_seq2_seq4_Go_2020, Microhabitat %in% metadata_microhabitat_psO_jki_seq2_seq4_Go_2020 [,i])
  metadata_subs_psO_jki_seq2_seq4_Go_2020 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Go_2020))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2020 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Go_2020, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Go_2020$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020 <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2020$signif[1])
}

for (i in 1:ncol(metadata_microhabitat_psO_jki_seq2_seq4_Go_2021)){
  ps_subs_psO_jki_seq2_seq4_Go_2021 <- subset_samples(psO_jki_seq2_seq4_Go_2021, Microhabitat %in% metadata_microhabitat_psO_jki_seq2_seq4_Go_2021 [,i])
  metadata_subs_psO_jki_seq2_seq4_Go_2021 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Go_2021))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2021 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Go_2021, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Go_2021$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021 <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2021$signif[1])
}

for (i in 1:ncol(metadata_microhabitat_psO_jki_seq2_seq4_Ki_2020)){
  ps_subs_psO_jki_seq2_seq4_Ki_2020 <- subset_samples(psO_jki_seq2_seq4_Ki_2020, Microhabitat %in% metadata_microhabitat_psO_jki_seq2_seq4_Ki_2020 [,i])
  metadata_subs_psO_jki_seq2_seq4_Ki_2020 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Ki_2020))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2020 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Ki_2020, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Ki_2020$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020 <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2020$signif[1])
}

for (i in 1:ncol(metadata_microhabitat_psO_jki_seq2_seq4_Ki_2021)){
  ps_subs_psO_jki_seq2_seq4_Ki_2021 <- subset_samples(psO_jki_seq2_seq4_Ki_2021, Microhabitat %in% metadata_microhabitat_psO_jki_seq2_seq4_Ki_2021 [,i])
  metadata_subs_psO_jki_seq2_seq4_Ki_2021 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Ki_2021))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2021 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Ki_2021, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Ki_2021$Microhabitat)
  p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021 <- c(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2021$signif[1])
}


#Adjust statistics values using BH
p_adj_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- p.adjust(p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4, method = "BH")
p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go, method = "BH")
p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki, method = "BH")

p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020 <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020, method = "BH")
p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021 <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021, method = "BH")

p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020 <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020, method = "BH")
p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021 <- p.adjust(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021, method = "BH")


#Create a table with statistics results
p_table_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- cbind.data.frame(t(metadata_site_year_psO_jki_seq2_seq4), p=p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4, p.adj=p_adj_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4)
p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq2_seq4_Go), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go)
p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq2_seq4_Ki), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki)

p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020 <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq2_seq4_Go_2020), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020)
p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021 <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq2_seq4_Go_2021), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021)

p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020 <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq2_seq4_Ki_2020), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020)
p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021 <- cbind.data.frame(t(metadata_microhabitat_psO_jki_seq2_seq4_Ki_2021), p=p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021, p.adj=p_adj_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021)

#print results
print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4)
print(p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4)
print(p_table_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki)


print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2020)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2020)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2021)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Go_2021)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2020)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2020)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2021)
print(p_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021)
print(p_table_pairwise_anosim_bray_microhabitat_psO_jki_seq2_seq4_Ki_2021)


## The end!!  : )


#### ANOSIM    -------> ROTATIONS
## testing of significance for the bray-curtis dissimilarity using ANOSIM
#calculate distance values between samples
dist_psO_jki_seq2_seq4 = phyloseq::distance(psO_jki_seq2_seq4, method = "bray")

dist_psO_jki_seq2_seq4_Go = phyloseq::distance(psO_jki_seq2_seq4_Go_filt, method = "bray")
dist_psO_jki_seq2_seq4_Ki = phyloseq::distance(psO_jki_seq2_seq4_Ki_filt, method = "bray")

dist_psO_jki_seq2_seq4_Go_2020 = phyloseq::distance(psO_jki_seq2_seq4_Go_2020_filt, method = "bray")
dist_psO_jki_seq2_seq4_Go_2021 = phyloseq::distance(psO_jki_seq2_seq4_Go_2021_filt, method = "bray")

dist_psO_jki_seq2_seq4_Ki_2020 = phyloseq::distance(psO_jki_seq2_seq4_Ki_2020_filt, method = "bray")
dist_psO_jki_seq2_seq4_Ki_2021 = phyloseq::distance(psO_jki_seq2_seq4_Ki_2021_filt, method = "bray")

#create a dataframe of the metadata
metadata_psO_jki_seq2_seq4 <- data.frame(sample_data(psO_jki_seq2_seq4))
metadata_psO_jki_seq2_seq4_Go <- data.frame(sample_data(psO_jki_seq2_seq4_Go_filt))
metadata_psO_jki_seq2_seq4_Ki <- data.frame(sample_data(psO_jki_seq2_seq4_Ki_filt))

metadata_psO_jki_seq2_seq4_Go_2020 <- data.frame(sample_data(psO_jki_seq2_seq4_Go_2020_filt))
metadata_psO_jki_seq2_seq4_Go_2021 <- data.frame(sample_data(psO_jki_seq2_seq4_Go_2021_filt))

metadata_psO_jki_seq2_seq4_Ki_2020 <- data.frame(sample_data(psO_jki_seq2_seq4_Ki_2020_filt))
metadata_psO_jki_seq2_seq4_Ki_2021 <- data.frame(sample_data(psO_jki_seq2_seq4_Ki_2021_filt))


#run Anosim using 10000 permutations
anosim_dist_psO_jki_seq2_seq4_site_year <- anosim(dist_psO_jki_seq2_seq4, metadata_psO_jki_seq2_seq4$Site_year, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Go_Rotation <- anosim(dist_psO_jki_seq2_seq4_Go, metadata_psO_jki_seq2_seq4_Go$Rotation, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Ki_Rotation <- anosim(dist_psO_jki_seq2_seq4_Ki, metadata_psO_jki_seq2_seq4_Ki$Rotation, permutations = 10000)

anosim_dist_psO_jki_seq2_seq4_Go_2020_Rotation <- anosim(dist_psO_jki_seq2_seq4_Go_2020, metadata_psO_jki_seq2_seq4_Go_2020$Rotation, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Go_2021_Rotation <- anosim(dist_psO_jki_seq2_seq4_Go_2021, metadata_psO_jki_seq2_seq4_Go_2021$Rotation, permutations = 10000)

anosim_dist_psO_jki_seq2_seq4_Ki_2020_Rotation <- anosim(dist_psO_jki_seq2_seq4_Ki_2020, metadata_psO_jki_seq2_seq4_Ki_2020$Rotation, permutations = 10000)
anosim_dist_psO_jki_seq2_seq4_Ki_2021_Rotation <- anosim(dist_psO_jki_seq2_seq4_Ki_2021, metadata_psO_jki_seq2_seq4_Ki_2021$Rotation, permutations = 10000)

#print results
print(anosim_dist_psO_jki_seq2_seq4_site_year)
print(anosim_dist_psO_jki_seq2_seq4_Go_Rotation)
print(anosim_dist_psO_jki_seq2_seq4_Ki_Rotation)

print(anosim_dist_psO_jki_seq2_seq4_Go_2020_Rotation)
print(anosim_dist_psO_jki_seq2_seq4_Go_2021_Rotation)

print(anosim_dist_psO_jki_seq2_seq4_Ki_2020_Rotation)
print(anosim_dist_psO_jki_seq2_seq4_Ki_2021_Rotation)


#### PAIRWISE ANOSIM
## testing of pairwise significance for the bray-curtis dissimilarity using ANOSIM
#Create metadata with an unique variable and run Anosim using 10000 permutations
metadata_site_year_psO_jki_seq2_seq4 <- combn(x=unique(metadata_psO_jki_seq2_seq4$Site_year), m=2)
p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- c()

metadata_rotation_psO_jki_seq2_seq4_Go <- combn(x=unique(metadata_psO_jki_seq2_seq4_Go$Rotation), m=2)
p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go <- c()

metadata_rotation_psO_jki_seq2_seq4_Ki <- combn(x=unique(metadata_psO_jki_seq2_seq4_Ki$Rotation), m=2)
p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki <- c()

metadata_rotation_psO_jki_seq2_seq4_Go_2020 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Go_2020$Rotation), m=2)
p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020 <- c()
metadata_rotation_psO_jki_seq2_seq4_Go_2021 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Go_2021$Rotation), m=2)
p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021 <- c()

metadata_rotation_psO_jki_seq2_seq4_Ki_2020 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Ki_2020$Rotation), m=2)
p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020 <- c()
metadata_rotation_psO_jki_seq2_seq4_Ki_2021 <- combn(x=unique(metadata_psO_jki_seq2_seq4_Ki_2021$Rotation), m=2)
p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021 <- c()



for (i in 1:ncol(metadata_site_year_psO_jki_seq2_seq4)){
  ps_subs_psO_jki_seq2_seq4 <- subset_samples(psO_jki_seq2_seq4, Site_year %in% metadata_site_year_psO_jki_seq2_seq4 [,i])
  metadata_subs_psO_jki_seq2_seq4 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4, method= "bray"), metadata_subs_psO_jki_seq2_seq4$Site_year)
  p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- c(p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4, pairwise_anosim_bray_dist_psO_jki_seq2_seq4$signif[1])
}

for (i in 1:ncol(metadata_rotation_psO_jki_seq2_seq4_Go)){
  ps_subs_psO_jki_seq2_seq4_Go <- subset_samples(psO_jki_seq2_seq4_Go, Rotation %in% metadata_rotation_psO_jki_seq2_seq4_Go [,i])
  metadata_subs_psO_jki_seq2_seq4_Go <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Go))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Go, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Go$Rotation)
  p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go <- c(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go$signif[1])
}

for (i in 1:ncol(metadata_rotation_psO_jki_seq2_seq4_Ki)){
  ps_subs_psO_jki_seq2_seq4_Ki <- subset_samples(psO_jki_seq2_seq4_Ki, Rotation %in% metadata_rotation_psO_jki_seq2_seq4_Ki [,i])
  metadata_subs_psO_jki_seq2_seq4_Ki <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Ki))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Ki, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Ki$Rotation)
  p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki <- c(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki$signif[1])
}


for (i in 1:ncol(metadata_rotation_psO_jki_seq2_seq4_Go_2020)){
  ps_subs_psO_jki_seq2_seq4_Go_2020 <- subset_samples(psO_jki_seq2_seq4_Go_2020, Rotation %in% metadata_rotation_psO_jki_seq2_seq4_Go_2020 [,i])
  metadata_subs_psO_jki_seq2_seq4_Go_2020 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Go_2020))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2020 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Go_2020, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Go_2020$Rotation)
  p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020 <- c(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2020$signif[1])
}

for (i in 1:ncol(metadata_rotation_psO_jki_seq2_seq4_Go_2021)){
  ps_subs_psO_jki_seq2_seq4_Go_2021 <- subset_samples(psO_jki_seq2_seq4_Go_2021, Rotation %in% metadata_rotation_psO_jki_seq2_seq4_Go_2021 [,i])
  metadata_subs_psO_jki_seq2_seq4_Go_2021 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Go_2021))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2021 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Go_2021, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Go_2021$Rotation)
  p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021 <- c(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2021$signif[1])
}

for (i in 1:ncol(metadata_rotation_psO_jki_seq2_seq4_Ki_2020)){
  ps_subs_psO_jki_seq2_seq4_Ki_2020 <- subset_samples(psO_jki_seq2_seq4_Ki_2020, Rotation %in% metadata_rotation_psO_jki_seq2_seq4_Ki_2020 [,i])
  metadata_subs_psO_jki_seq2_seq4_Ki_2020 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Ki_2020))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2020 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Ki_2020, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Ki_2020$Rotation)
  p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020 <- c(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2020$signif[1])
}

for (i in 1:ncol(metadata_rotation_psO_jki_seq2_seq4_Ki_2021)){
  ps_subs_psO_jki_seq2_seq4_Ki_2021 <- subset_samples(psO_jki_seq2_seq4_Ki_2021, Rotation %in% metadata_rotation_psO_jki_seq2_seq4_Ki_2021 [,i])
  metadata_subs_psO_jki_seq2_seq4_Ki_2021 <- data.frame(sample_data(ps_subs_psO_jki_seq2_seq4_Ki_2021))
  pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2021 <- anosim(phyloseq::distance(ps_subs_psO_jki_seq2_seq4_Ki_2021, method= "bray"), metadata_subs_psO_jki_seq2_seq4_Ki_2021$Rotation)
  p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021 <- c(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021, pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2021$signif[1])
}


#Adjust statistics values using BH
p_adj_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- p.adjust(p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4, method = "BH")
p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go <- p.adjust(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go, method = "BH")
p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki <- p.adjust(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki, method = "BH")

p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020 <- p.adjust(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020, method = "BH")
p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021 <- p.adjust(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021, method = "BH")

p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020 <- p.adjust(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020, method = "BH")
p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021 <- p.adjust(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021, method = "BH")


#Create a table with statistics results
p_table_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4 <- cbind.data.frame(t(metadata_site_year_psO_jki_seq2_seq4), p=p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4, p.adj=p_adj_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4)
p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go <- cbind.data.frame(t(metadata_rotation_psO_jki_seq2_seq4_Go), p=p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go, p.adj=p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go)
p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki <- cbind.data.frame(t(metadata_rotation_psO_jki_seq2_seq4_Ki), p=p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki, p.adj=p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki)

p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020 <- cbind.data.frame(t(metadata_rotation_psO_jki_seq2_seq4_Go_2020), p=p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020, p.adj=p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020)
p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021 <- cbind.data.frame(t(metadata_rotation_psO_jki_seq2_seq4_Go_2021), p=p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021, p.adj=p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021)

p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020 <- cbind.data.frame(t(metadata_rotation_psO_jki_seq2_seq4_Ki_2020), p=p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020, p.adj=p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020)
p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021 <- cbind.data.frame(t(metadata_rotation_psO_jki_seq2_seq4_Ki_2021), p=p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021, p.adj=p_adj_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021)

#print results
print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4)
print(p_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4)
print(p_table_pairwise_anosim_bray_site_year_psO_jki_seq2_seq4)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go)
print(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go)
print(p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki)
print(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki)
print(p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki)


print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2020)
print(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020)
print(p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2020)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Go_2021)
print(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021)
print(p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Go_2021)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2020)
print(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020)
print(p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2020)

print(pairwise_anosim_bray_dist_psO_jki_seq2_seq4_Ki_2021)
print(p_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021)
print(p_table_pairwise_anosim_bray_rotation_psO_jki_seq2_seq4_Ki_2021)



#### The end! Have fun!!
