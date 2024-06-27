##Creating MDS plot for selected dataset
#load package
library("phyloseq")
library("vegan")
library("ggplot2")
library("dplyr")
library("microbiome")
library("ggpubr")
library("ggordiplots")

# #source("https://raw.githubusercontent.com/devanmcg/IntroRangeR/master/11_IntroMultivariate/ordinationsGGplot.R")


################################   MDS  ################################
#set seed
set.seed(2022)
psO_jki_seq10
psO_jki_seq10_rarefied


#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
MDS_bray_psO_jki_seq10_rarefied_sqr<-ordinate(psO_jki_seq10_rarefied, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq10_rarefied_Go_filt_sqr<-ordinate(psO_jki_seq10_rarefied_Go_filt, "MDS","bray", autotransform=TRUE)
MDS_bray_psO_jki_seq10_rarefied_Ki_filt_sqr<-ordinate(psO_jki_seq10_rarefied_Ki_filt, "MDS","bray", autotransform=TRUE)

#Print stress data, dimensions and number of tries
head(MDS_bray_psO_jki_seq10_rarefied_sqr)
head(MDS_bray_psO_jki_seq10_rarefied_Go_filt_sqr)
head(MDS_bray_psO_jki_seq10_rarefied_Ki_filt_sqr)

#Create a MDS plot
####
plot_MDS_bray_psO_jki_seq10_rarefied_sqr <- plot_ordination(psO_jki_seq10_rarefied, MDS_bray_psO_jki_seq10_rarefied_sqr, type = "sample", color = "Site_rot") +
  geom_point(size = 4) +
  stat_ellipse(type = "norm", linetype = 2, geom = "polygon", alpha=0) +
  theme_bw() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 20)) +
  scale_color_manual(values = c("Site" <- c("#6baea5", "#216974", "#ebaf7a", "#D1711F")))
plot_MDS_bray_psO_jki_seq10_rarefied_sqr

ggsave("plot_MDS_bray_psO_jki_seq10_rarefied_sqr.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Ordination_jki_seq10/", width = 15, height = 11, units = "cm",dpi = 300)



## Individual soils
plot_MDS_bray_psO_jki_seq10_rarefied_Go_filt_sqr <- plot_ordination(psO_jki_seq10_rarefied_Go_filt, MDS_bray_psO_jki_seq10_rarefied_Go_filt_sqr, type = "sample", color = "Site_rot", shape = "Site_rot_tre") +
  geom_point(size = 4) +
  stat_ellipse(type = "norm", linetype = 2, geom = "polygon", alpha=0) +
  theme_bw() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 20)) +
  scale_color_manual(
    values = c("#6baea5", "#216974"),
    limits = c( "Go_W1", "Go_WM")) +
  scale_shape_manual(values=c(16, 17, 15, 18))
plot_MDS_bray_psO_jki_seq10_rarefied_Go_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq10_rarefied_Go_filt_sqr.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Ordination_jki_seq10/", width = 13, height = 7, units = "cm",dpi = 300)

plot_MDS_bray_psO_jki_seq10_rarefied_Ki_filt_sqr <- plot_ordination(psO_jki_seq10_rarefied_Ki_filt, MDS_bray_psO_jki_seq10_rarefied_Ki_filt_sqr, type = "sample", color = "Site_rot", shape = "Site_rot_tre") +
  geom_point(size = 4) +
  stat_ellipse(type = "norm", linetype = 2, geom = "polygon", alpha=0) +
  theme_bw() +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 20)) +
  scale_color_manual(
    values = c("#ebaf7a", "#D1711F"),
    limits = c("Ki_W1", "Ki_W3")) +
  scale_shape_manual(values=c(16, 17, 15, 18))
plot_MDS_bray_psO_jki_seq10_rarefied_Ki_filt_sqr

ggsave("plot_MDS_bray_psO_jki_seq10_rarefied_Ki_filt_sqr.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Ordination_jki_seq10/", width = 13, height = 7, units = "cm",dpi = 300)




################################## STATISTICS
set.seed(2022)
##Multivariate Anova for ordinate files

#Transform abundances to square root
#psO_jki_seq10 <- microbiome::transform(psO_jki_seq10, "hellinger")

#Convert phyloseq object to dataframe using abundances function and meta function from microbiome package
psO_jki_seq10_rarefied_abundances <- abundances(psO_jki_seq10_rarefied)
psO_jki_seq10_rarefied_Go_filt_abundances <- abundances(psO_jki_seq10_rarefied_Go_filt)
psO_jki_seq10_rarefied_Ki_filt_abundances <- abundances(psO_jki_seq10_rarefied_Ki_filt)

psO_jki_seq10_rarefied_meta <- meta(psO_jki_seq10_rarefied)
psO_jki_seq10_rarefied_Go_filt_meta <- meta(psO_jki_seq10_rarefied_Go_filt)
psO_jki_seq10_rarefied_Ki_filt_meta <- meta(psO_jki_seq10_rarefied_Ki_filt)

#######Permanova using adonis function from vegan and print p value
Permanova_psO_jki_seq10_rarefied_abundances <- adonis2(t(psO_jki_seq10_rarefied_abundances) ~ Site*Rotation*Site_rot*Treatment, data = psO_jki_seq10_rarefied_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq10_rarefied_abundances

write.csv(Permanova_psO_jki_seq10_rarefied_abundances, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/Permanova_psO_jki_seq10_rarefied_abundances_MDS.csv")


#######Permanova using adonis function from vegan and print p value
Permanova_psO_jki_seq10_rarefied_Go_filt_abundances <- adonis2(t(psO_jki_seq10_rarefied_Go_filt_abundances) ~ Rotation*Site_rot*Treatment, data = psO_jki_seq10_rarefied_Go_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq10_rarefied_Go_filt_abundances

write.csv(Permanova_psO_jki_seq10_rarefied_Go_filt_abundances, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/Permanova_psO_jki_seq10_rarefied_Go_filt_abundances.csv")


#######Permanova using adonis function from vegan and print p value
Permanova_psO_jki_seq10_rarefied_Ki_filt_abundances <- adonis2(t(psO_jki_seq10_rarefied_Ki_filt_abundances) ~ Rotation*Site_rot*Treatment, data = psO_jki_seq10_rarefied_Ki_filt_meta, permutations = 10000, method = "bray")
Permanova_psO_jki_seq10_rarefied_Ki_filt_abundances

write.csv(Permanova_psO_jki_seq10_rarefied_Ki_filt_abundances, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/Permanova_psO_jki_seq10_rarefied_Ki_filt_abundances.csv")



#### ANOSIM    -------> Site_rots
## testing of significance for the bray-curtis dissimilarity using ANOSIM
#calculate distance values between samples
dist_psO_jki_seq10_rarefied = phyloseq::distance(psO_jki_seq10_rarefied, method = "bray")

dist_psO_jki_seq10_rarefied_Go = phyloseq::distance(psO_jki_seq10_rarefied_Go_filt, method = "bray")
dist_psO_jki_seq10_rarefied_Ki = phyloseq::distance(psO_jki_seq10_rarefied_Ki_filt, method = "bray")

#create a dataframe of the metadata
metadata_psO_jki_seq10_rarefied <- data.frame(sample_data(psO_jki_seq10_rarefied))
metadata_psO_jki_seq10_rarefied_Go <- data.frame(sample_data(psO_jki_seq10_rarefied_Go_filt))
metadata_psO_jki_seq10_rarefied_Ki <- data.frame(sample_data(psO_jki_seq10_rarefied_Ki_filt))

#run Anosim using 10000 permutations
anosim_dist_psO_jki_seq10_rarefied_Site_rot <- anosim(dist_psO_jki_seq10_rarefied, metadata_psO_jki_seq10_rarefied$Site_rot, permutations = 10000)
anosim_dist_psO_jki_seq10_rarefied_Go_Site_rot <- anosim(dist_psO_jki_seq10_rarefied_Go, metadata_psO_jki_seq10_rarefied_Go$Site_rot, permutations = 10000)
anosim_dist_psO_jki_seq10_rarefied_Ki_Site_rot <- anosim(dist_psO_jki_seq10_rarefied_Ki, metadata_psO_jki_seq10_rarefied_Ki$Site_rot, permutations = 10000)

#print results
print(anosim_dist_psO_jki_seq10_rarefied_Site_rot)
print(anosim_dist_psO_jki_seq10_rarefied_Go_Site_rot)
print(anosim_dist_psO_jki_seq10_rarefied_Ki_Site_rot)

#### PAIRWISE ANOSIM
## testing of pairwise significance for the bray-curtis dissimilarity using ANOSIM
#Create metadata with an unique variable and run Anosim using 10000 permutations
metadata_Site_rot_tre_psO_jki_seq10_rarefied <- combn(x=unique(metadata_psO_jki_seq10_rarefied$Site_rot_tre), m=2)
p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied <- c()

metadata_Site_rot_tre_psO_jki_seq10_rarefied_Go <- combn(x=unique(metadata_psO_jki_seq10_rarefied_Go$Site_rot_tre), m=2)
p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go <- c()

metadata_Site_rot_tre_psO_jki_seq10_rarefied_Ki <- combn(x=unique(metadata_psO_jki_seq10_rarefied_Ki$Site_rot_tre), m=2)
p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki <- c()


for (i in 1:ncol(metadata_Site_rot_tre_psO_jki_seq10_rarefied)){
  ps_subs_psO_jki_seq10_rarefied <- subset_samples(psO_jki_seq10_rarefied, Site_rot_tre %in% metadata_Site_rot_tre_psO_jki_seq10_rarefied [,i])
  metadata_subs_psO_jki_seq10_rarefied <- data.frame(sample_data(ps_subs_psO_jki_seq10_rarefied))
  pairwise_anosim_bray_dist_psO_jki_seq10_rarefied <- anosim(phyloseq::distance(ps_subs_psO_jki_seq10_rarefied, method= "bray"), metadata_subs_psO_jki_seq10_rarefied$Site_rot_tre)
  p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied <- c(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied, pairwise_anosim_bray_dist_psO_jki_seq10_rarefied$signif[1])
}


for (i in 1:ncol(metadata_Site_rot_tre_psO_jki_seq10_rarefied_Go)){
  ps_subs_psO_jki_seq10_rarefied_Go <- subset_samples(psO_jki_seq10_rarefied_Go, Site_rot_tre %in% metadata_Site_rot_tre_psO_jki_seq10_rarefied_Go [,i])
  metadata_subs_psO_jki_seq10_rarefied_Go <- data.frame(sample_data(ps_subs_psO_jki_seq10_rarefied_Go))
  pairwise_anosim_bray_dist_psO_jki_seq10_rarefied_Go <- anosim(phyloseq::distance(ps_subs_psO_jki_seq10_rarefied_Go, method= "bray"), metadata_subs_psO_jki_seq10_rarefied_Go$Site_rot_tre)
  p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go <- c(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go, pairwise_anosim_bray_dist_psO_jki_seq10_rarefied_Go$signif[1])
}

for (i in 1:ncol(metadata_Site_rot_tre_psO_jki_seq10_rarefied_Ki)){
  ps_subs_psO_jki_seq10_rarefied_Ki <- subset_samples(psO_jki_seq10_rarefied_Ki, Site_rot_tre %in% metadata_Site_rot_tre_psO_jki_seq10_rarefied_Ki [,i])
  metadata_subs_psO_jki_seq10_rarefied_Ki <- data.frame(sample_data(ps_subs_psO_jki_seq10_rarefied_Ki))
  pairwise_anosim_bray_dist_psO_jki_seq10_rarefied_Ki <- anosim(phyloseq::distance(ps_subs_psO_jki_seq10_rarefied_Ki, method= "bray"), metadata_subs_psO_jki_seq10_rarefied_Ki$Site_rot_tre)
  p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki <- c(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki, pairwise_anosim_bray_dist_psO_jki_seq10_rarefied_Ki$signif[1])
}

#Adjust statistics values using BH
p_adj_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied <- p.adjust(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied, method = "BH")
p_adj_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go <- p.adjust(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go, method = "BH")
p_adj_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki <- p.adjust(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki, method = "BH")

#Create a table with statistics results
p_table_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied <- cbind.data.frame(t(metadata_Site_rot_tre_psO_jki_seq10_rarefied), p=p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied, p.adj=p_adj_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied)
p_table_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go <- cbind.data.frame(t(metadata_Site_rot_tre_psO_jki_seq10_rarefied_Go), p=p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go, p.adj=p_adj_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go)
p_table_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki <- cbind.data.frame(t(metadata_Site_rot_tre_psO_jki_seq10_rarefied_Ki), p=p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki, p.adj=p_adj_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki)

#print results
print(pairwise_anosim_bray_dist_psO_jki_seq10_rarefied)
print(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied)
print(p_table_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied)

print(pairwise_anosim_bray_dist_psO_jki_seq10_rarefied_Go)
print(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go)
print(p_table_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Go)

print(pairwise_anosim_bray_dist_psO_jki_seq10_rarefied_Ki)
print(p_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki)
print(p_table_pairwise_anosim_bray_Site_rot_tre_psO_jki_seq10_rarefied_Ki)

# The end! Enjoy!


