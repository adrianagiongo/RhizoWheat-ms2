###Select groups and subset data from the cleaned dataset
##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

#Calculate number of genera present after filtration
length(get_taxa_unique(psO_jki_seq10, taxonomic.rank = "Phylum"))
length(get_taxa_unique(psO_jki_seq10, taxonomic.rank = "Annotation"))

#Agglomerate taxa using tax glom function and removing NAs
psO_jki_seq10_phylum <- tax_glom(psO_jki_seq10, "Phylum", NArm= TRUE)
psO_jki_seq10_phylum
psO_jki_seq10_annotation <- tax_glom(psO_jki_seq10, "Annotation", NArm= TRUE)
psO_jki_seq10_annotation
ntaxa(psO_jki_seq10); ntaxa(psO_jki_seq10_annotation)

#### Create tables 
df_psO_jki_seq10_phylum <- data.frame(tax_table(psO_jki_seq10_phylum),otu_table(psO_jki_seq10_phylum))
write.csv(df_psO_jki_seq10_phylum, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_phylum.csv")
psO_jki_seq10_phylum_rel<-transform_sample_counts(psO_jki_seq10_phylum, function(x) (x*100)/sum(x))
psO_jki_seq10_phylum_rel
head(otu_table(psO_jki_seq10_phylum_rel))
df_psO_jki_seq10_phylum_rel <- data.frame(tax_table(psO_jki_seq10_phylum_rel),otu_table(psO_jki_seq10_phylum_rel))
write.csv(df_psO_jki_seq10_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_phylum_rel.csv")

#For psO_jki_seq10_annotation
df_psO_jki_seq10_annotation <- data.frame(tax_table(psO_jki_seq10_annotation),otu_table(psO_jki_seq10_annotation))
write.csv(df_psO_jki_seq10_annotation, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_annotation.csv")
psO_jki_seq10_annotation_rel<-transform_sample_counts(psO_jki_seq10_annotation, function(x) (x*100)/sum(x))
psO_jki_seq10_annotation_rel
head(otu_table(psO_jki_seq10_annotation_rel))
df_psO_jki_seq10_annotation_rel <- data.frame(tax_table(psO_jki_seq10_annotation_rel),otu_table(psO_jki_seq10_annotation_rel))
write.csv(df_psO_jki_seq10_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_annotation_rel.csv")

#######################################
#### Separated by Site
#######################################
#Go
psO_jki_seq10_Go <- subset_samples(psO_jki_seq10, Site == "Go")
psO_jki_seq10_Go
psO_jki_seq10_Go_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Go) > 0, psO_jki_seq10_Go) ### Take the samples with 0 reads off
psO_jki_seq10_Go_filt
colnames(tax_table(psO_jki_seq10_Go_filt))

##Phylum
psO_jki_seq10_Go_filt_phylum <- tax_glom(psO_jki_seq10_Go_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Go_filt); ntaxa(psO_jki_seq10_Go_filt_phylum)
psO_jki_seq10_Go_filt_phylum
psO_jki_seq10_Go_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Go_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Go_filt_phylum_rel
head(otu_table(psO_jki_seq10_Go_filt_phylum_rel))
df_psO_jki_seq10_Go_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Go_filt_phylum_rel),otu_table(psO_jki_seq10_Go_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Go_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Go_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Go_filt_annotation <- tax_glom(psO_jki_seq10_Go_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Go_filt); ntaxa(psO_jki_seq10_Go_filt_annotation)
psO_jki_seq10_Go_filt_annotation
psO_jki_seq10_Go_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Go_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Go_filt_annotation_rel
head(otu_table(psO_jki_seq10_Go_filt_annotation_rel))
df_psO_jki_seq10_Go_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Go_filt_annotation_rel),otu_table(psO_jki_seq10_Go_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Go_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Go_filt_annotation_rel.csv")


#Ki
psO_jki_seq10_Ki <- subset_samples(psO_jki_seq10, Site == "Ki")
psO_jki_seq10_Ki
psO_jki_seq10_Ki_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Ki) > 0, psO_jki_seq10_Ki) ### Take the samples with 0 reads off
psO_jki_seq10_Ki_filt
colnames(tax_table(psO_jki_seq10_Ki_filt))

##Phylum
psO_jki_seq10_Ki_filt_phylum <- tax_glom(psO_jki_seq10_Ki_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_filt); ntaxa(psO_jki_seq10_Ki_filt_phylum)
psO_jki_seq10_Ki_filt_phylum
psO_jki_seq10_Ki_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Ki_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_filt_phylum_rel
head(otu_table(psO_jki_seq10_Ki_filt_phylum_rel))
df_psO_jki_seq10_Ki_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Ki_filt_phylum_rel),otu_table(psO_jki_seq10_Ki_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Ki_filt_annotation <- tax_glom(psO_jki_seq10_Ki_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_filt); ntaxa(psO_jki_seq10_Ki_filt_annotation)
psO_jki_seq10_Ki_filt_annotation
psO_jki_seq10_Ki_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Ki_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_filt_annotation_rel
head(otu_table(psO_jki_seq10_Ki_filt_annotation_rel))
df_psO_jki_seq10_Ki_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Ki_filt_annotation_rel),otu_table(psO_jki_seq10_Ki_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_filt_annotation_rel.csv")

### I stopped here (12.02)
#######################################
#### Separated by Site and Soil
#######################################
#Go
#psO_jki_seq10_Go <- subset_samples(psO_jki_seq10, Site == "Go")
psO_jki_seq10_Go_W1 <- subset_samples(psO_jki_seq10_Go, Rotation == "W1")
psO_jki_seq10_Go_WM <- subset_samples(psO_jki_seq10_Go, Rotation == "WM")

##Go_W1
psO_jki_seq10_Go_W1
psO_jki_seq10_Go_W1_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Go_W1) > 0, psO_jki_seq10_Go_W1) ### Take the samples with 0 reads off
psO_jki_seq10_Go_W1_filt
colnames(tax_table(psO_jki_seq10_Go_W1_filt))

##Phylum
psO_jki_seq10_Go_W1_filt_phylum <- tax_glom(psO_jki_seq10_Go_W1_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Go_W1_filt); ntaxa(psO_jki_seq10_Go_W1_filt_phylum)
psO_jki_seq10_Go_W1_filt_phylum
psO_jki_seq10_Go_W1_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Go_W1_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Go_W1_filt_phylum_rel
head(otu_table(psO_jki_seq10_Go_W1_filt_phylum_rel))
df_psO_jki_seq10_Go_W1_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Go_W1_filt_phylum_rel),otu_table(psO_jki_seq10_Go_W1_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Go_W1_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Go_W1_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Go_W1_filt_annotation <- tax_glom(psO_jki_seq10_Go_W1_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Go_W1_filt); ntaxa(psO_jki_seq10_Go_W1_filt_annotation)
psO_jki_seq10_Go_W1_filt_annotation
psO_jki_seq10_Go_W1_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Go_W1_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Go_W1_filt_annotation_rel
head(otu_table(psO_jki_seq10_Go_W1_filt_annotation_rel))
df_psO_jki_seq10_Go_W1_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Go_W1_filt_annotation_rel),otu_table(psO_jki_seq10_Go_W1_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Go_W1_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Go_W1_filt_annotation_rel.csv")

##Go_WM
psO_jki_seq10_Go_WM
psO_jki_seq10_Go_WM_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Go_WM) > 0, psO_jki_seq10_Go_WM) ### Take the samples with 0 reads off
psO_jki_seq10_Go_WM_filt
colnames(tax_table(psO_jki_seq10_Go_WM_filt))

##Phylum
psO_jki_seq10_Go_WM_filt_phylum <- tax_glom(psO_jki_seq10_Go_WM_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Go_WM_filt); ntaxa(psO_jki_seq10_Go_WM_filt_phylum)
psO_jki_seq10_Go_WM_filt_phylum
psO_jki_seq10_Go_WM_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Go_WM_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Go_WM_filt_phylum_rel
head(otu_table(psO_jki_seq10_Go_WM_filt_phylum_rel))
df_psO_jki_seq10_Go_WM_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Go_WM_filt_phylum_rel),otu_table(psO_jki_seq10_Go_WM_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Go_WM_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Go_WM_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Go_WM_filt_annotation <- tax_glom(psO_jki_seq10_Go_WM_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Go_WM_filt); ntaxa(psO_jki_seq10_Go_WM_filt_annotation)
psO_jki_seq10_Go_WM_filt_annotation
psO_jki_seq10_Go_WM_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Go_WM_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Go_WM_filt_annotation_rel
head(otu_table(psO_jki_seq10_Go_WM_filt_annotation_rel))
df_psO_jki_seq10_Go_WM_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Go_WM_filt_annotation_rel),otu_table(psO_jki_seq10_Go_WM_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Go_WM_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Go_WM_filt_annotation_rel.csv")


#Ki
psO_jki_seq10_Ki <- subset_samples(psO_jki_seq10, Site == "Ki")
psO_jki_seq10_Ki
psO_jki_seq10_Ki_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Ki) > 0, psO_jki_seq10_Ki) ### Take the samples with 0 reads off
psO_jki_seq10_Ki_filt
colnames(tax_table(psO_jki_seq10_Ki_filt))

##Phylum
psO_jki_seq10_Ki_filt_phylum <- tax_glom(psO_jki_seq10_Ki_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_filt); ntaxa(psO_jki_seq10_Ki_filt_phylum)
psO_jki_seq10_Ki_filt_phylum
psO_jki_seq10_Ki_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Ki_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_filt_phylum_rel
head(otu_table(psO_jki_seq10_Ki_filt_phylum_rel))
df_psO_jki_seq10_Ki_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Ki_filt_phylum_rel),otu_table(psO_jki_seq10_Ki_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Ki_filt_annotation <- tax_glom(psO_jki_seq10_Ki_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_filt); ntaxa(psO_jki_seq10_Ki_filt_annotation)
psO_jki_seq10_Ki_filt_annotation
psO_jki_seq10_Ki_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Ki_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_filt_annotation_rel
head(otu_table(psO_jki_seq10_Ki_filt_annotation_rel))
df_psO_jki_seq10_Ki_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Ki_filt_annotation_rel),otu_table(psO_jki_seq10_Ki_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_filt_annotation_rel.csv")


#Ki
#psO_jki_seq10_Ki <- subset_samples(psO_jki_seq10, Site == "Ki")
psO_jki_seq10_Ki_W1 <- subset_samples(psO_jki_seq10_Ki, Rotation == "W1")
psO_jki_seq10_Ki_W3 <- subset_samples(psO_jki_seq10_Ki, Rotation == "W3")

##Ki_W1
psO_jki_seq10_Ki_W1
psO_jki_seq10_Ki_W1_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Ki_W1) > 0, psO_jki_seq10_Ki_W1) ### Take the samples with 0 reads off
psO_jki_seq10_Ki_W1_filt
colnames(tax_table(psO_jki_seq10_Ki_W1_filt))

##Phylum
psO_jki_seq10_Ki_W1_filt_phylum <- tax_glom(psO_jki_seq10_Ki_W1_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_W1_filt); ntaxa(psO_jki_seq10_Ki_W1_filt_phylum)
psO_jki_seq10_Ki_W1_filt_phylum
psO_jki_seq10_Ki_W1_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Ki_W1_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_W1_filt_phylum_rel
head(otu_table(psO_jki_seq10_Ki_W1_filt_phylum_rel))
df_psO_jki_seq10_Ki_W1_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Ki_W1_filt_phylum_rel),otu_table(psO_jki_seq10_Ki_W1_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_W1_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_W1_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Ki_W1_filt_annotation <- tax_glom(psO_jki_seq10_Ki_W1_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_W1_filt); ntaxa(psO_jki_seq10_Ki_W1_filt_annotation)
psO_jki_seq10_Ki_W1_filt_annotation
psO_jki_seq10_Ki_W1_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Ki_W1_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_W1_filt_annotation_rel
head(otu_table(psO_jki_seq10_Ki_W1_filt_annotation_rel))
df_psO_jki_seq10_Ki_W1_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Ki_W1_filt_annotation_rel),otu_table(psO_jki_seq10_Ki_W1_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_W1_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_W1_filt_annotation_rel.csv")

##Ki_W3
psO_jki_seq10_Ki_W3
psO_jki_seq10_Ki_W3_filt<-prune_taxa(taxa_sums(psO_jki_seq10_Ki_W3) > 0, psO_jki_seq10_Ki_W3) ### Take the samples with 0 reads off
psO_jki_seq10_Ki_W3_filt
colnames(tax_table(psO_jki_seq10_Ki_W3_filt))

##Phylum
psO_jki_seq10_Ki_W3_filt_phylum <- tax_glom(psO_jki_seq10_Ki_W3_filt, taxrank = "Phylum")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_W3_filt); ntaxa(psO_jki_seq10_Ki_W3_filt_phylum)
psO_jki_seq10_Ki_W3_filt_phylum
psO_jki_seq10_Ki_W3_filt_phylum_rel<-transform_sample_counts(psO_jki_seq10_Ki_W3_filt_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_W3_filt_phylum_rel
head(otu_table(psO_jki_seq10_Ki_W3_filt_phylum_rel))
df_psO_jki_seq10_Ki_W3_filt_phylum_rel <- data.frame(tax_table(psO_jki_seq10_Ki_W3_filt_phylum_rel),otu_table(psO_jki_seq10_Ki_W3_filt_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_W3_filt_phylum_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_W3_filt_phylum_rel.csv")

##Annotation
psO_jki_seq10_Ki_W3_filt_annotation <- tax_glom(psO_jki_seq10_Ki_W3_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_jki_seq10_Ki_W3_filt); ntaxa(psO_jki_seq10_Ki_W3_filt_annotation)
psO_jki_seq10_Ki_W3_filt_annotation
psO_jki_seq10_Ki_W3_filt_annotation_rel<-transform_sample_counts(psO_jki_seq10_Ki_W3_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq10_Ki_W3_filt_annotation_rel
head(otu_table(psO_jki_seq10_Ki_W3_filt_annotation_rel))
df_psO_jki_seq10_Ki_W3_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq10_Ki_W3_filt_annotation_rel),otu_table(psO_jki_seq10_Ki_W3_filt_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq10_Ki_W3_filt_annotation_rel, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_Ki_W3_filt_annotation_rel.csv")



## The end!  : )