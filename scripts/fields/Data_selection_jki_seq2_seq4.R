###Select groups and subset data from the cleaned dataset

##use on a phyloseq object (large phyloseq data)
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

###Subsetting samples from cleaned dataset to keep only samples that represent every case on variables  
##Variables / colors
#(1)Site = Go or Ki
#(2)Year = 2020 or 2021 
#(3)Plant time = T1 or T2 or T3 
#(4)Rotation = W1, W2, WM (Go) or W1, W3 (Ki) - / (MetBrewer:"VanGogh3") W1= "#e7e5cc", W2= "#9cc184", W3= "#447243", WM= "#1e3d14"
#(5)Microhabitats = BS, RH, RP

#### All samples
psO_jki_seq2_seq4
sample_variables(psO_jki_seq2_seq4)

## Agglomerate in Phylum
colnames(tax_table(psO_jki_seq2_seq4))
psO_jki_seq2_seq4_phylum <- tax_glom(psO_jki_seq2_seq4, taxrank = "Phylum")
ntaxa(psO_jki_seq2_seq4); ntaxa(psO_jki_seq2_seq4_phylum)

df_psO_jki_seq2_seq4_phylum <- data.frame(tax_table(psO_jki_seq2_seq4_phylum), otu_table(psO_jki_seq2_seq4_phylum))
head(df_psO_jki_seq2_seq4_phylum)
write.csv(df_psO_jki_seq2_seq4_phylum, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_phylum.csv")

psO_jki_seq2_seq4_phylum_rel<-transform_sample_counts(psO_jki_seq2_seq4_phylum, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq2_seq4_phylum_rel
head(otu_table(psO_jki_seq2_seq4_phylum_rel))
df_psO_jki_seq2_seq4_phylum_rel <- data.frame(tax_table(psO_jki_seq2_seq4_phylum_rel),otu_table(psO_jki_seq2_seq4_phylum_rel)) ### Create tables
write.csv(df_psO_jki_seq2_seq4_phylum_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_phylum_rel.csv")

## Agglomerate in Annotation
colnames(tax_table(psO_jki_seq2_seq4))
psO_jki_seq2_seq4_annotation <- tax_glom(psO_jki_seq2_seq4, taxrank = "Annotation")
ntaxa(psO_jki_seq2_seq4); ntaxa(psO_jki_seq2_seq4_annotation)

df_psO_jki_seq2_seq4_annotation <- data.frame(tax_table(psO_jki_seq2_seq4_annotation), otu_table(psO_jki_seq2_seq4_annotation))
head(df_psO_jki_seq2_seq4_annotation)
write.csv(df_psO_jki_seq2_seq4_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_annotation.csv")

psO_jki_seq2_seq4_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_jki_seq2_seq4_annotation_rel
head(otu_table(psO_jki_seq2_seq4_annotation_rel))
df_psO_jki_seq2_seq4_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_annotation_rel),otu_table(psO_jki_seq2_seq4_annotation_rel)) ### Create tables
write.csv(df_psO_jki_seq2_seq4_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_annotation_rel.csv")


#######################################
#### Separated by Site
#######################################
psO_jki_seq2_seq4_Go <- subset_samples(psO_jki_seq2_seq4, Site == "Go")
  psO_jki_seq2_seq4_Go
psO_jki_seq2_seq4_Go_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go) > 0, psO_jki_seq2_seq4_Go) ### Take the samples with 0 reads off
  psO_jki_seq2_seq4_Go_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_filt))
psO_jki_seq2_seq4_Go_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
  ntaxa(psO_jki_seq2_seq4_Go_filt); ntaxa(psO_jki_seq2_seq4_Go_filt_annotation)
  psO_jki_seq2_seq4_Go_filt_annotation
psO_jki_seq2_seq4_Go_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
  psO_jki_seq2_seq4_Go_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_filt_annotation_rel)) ### Create tables
  write.csv(df_psO_jki_seq2_seq4_Go_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki <- subset_samples(psO_jki_seq2_seq4, Site == "Ki")
  psO_jki_seq2_seq4_Ki
psO_jki_seq2_seq4_Ki_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki) > 0, psO_jki_seq2_seq4_Ki) ### Take the samples with 0 reads off
  psO_jki_seq2_seq4_Ki_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_filt))
psO_jki_seq2_seq4_Ki_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
  ntaxa(psO_jki_seq2_seq4_Ki_filt); ntaxa(psO_jki_seq2_seq4_Ki_filt_annotation)
  psO_jki_seq2_seq4_Ki_filt_annotation
psO_jki_seq2_seq4_Ki_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
  psO_jki_seq2_seq4_Ki_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_filt_annotation_rel)) ### Create tables
  write.csv(df_psO_jki_seq2_seq4_Ki_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_filt_annotation_rel.csv")
  
  
#######################################
#### Separated by Site and Year
#######################################
psO_jki_seq2_seq4_Go_2020 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020")
  psO_jki_seq2_seq4_Go_2020
psO_jki_seq2_seq4_Go_2020_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020) > 0, psO_jki_seq2_seq4_Go_2020) ### Take the samples with 0 reads off
  psO_jki_seq2_seq4_Go_2020_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_filt))
psO_jki_seq2_seq4_Go_2020_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_filt, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
  ntaxa(psO_jki_seq2_seq4_Go_2020_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_filt_annotation
psO_jki_seq2_seq4_Go_2020_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_filt_annotation, function(x) (x*100)/sum(x)) ### Transform to relative abundance
  psO_jki_seq2_seq4_Go_2020_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2020_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_filt_annotation_rel)) ### Create tables
  write.csv(df_psO_jki_seq2_seq4_Go_2020_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Go_2021 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021")
  psO_jki_seq2_seq4_Go_2021
psO_jki_seq2_seq4_Go_2021_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021) > 0, psO_jki_seq2_seq4_Go_2021)
  psO_jki_seq2_seq4_Go_2021_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_filt))
psO_jki_seq2_seq4_Go_2021_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_filt_annotation
psO_jki_seq2_seq4_Go_2021_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2021_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2020 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020")
  psO_jki_seq2_seq4_Ki_2020
psO_jki_seq2_seq4_Ki_2020_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020) > 0, psO_jki_seq2_seq4_Ki_2020)
  psO_jki_seq2_seq4_Ki_2020_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_filt))
psO_jki_seq2_seq4_Ki_2020_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_filt_annotation
psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2021 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021")
  psO_jki_seq2_seq4_Ki_2021
psO_jki_seq2_seq4_Ki_2021_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021) > 0, psO_jki_seq2_seq4_Ki_2021)
  psO_jki_seq2_seq4_Ki_2021_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_filt))
psO_jki_seq2_seq4_Ki_2021_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_filt_annotation
psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_filt_annotation_rel.csv")

  
#######################################
#### Separated by Site and Year and T2
#######################################
  
psO_jki_seq2_seq4_Go_2020_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year == "Y2020" & Stage == "T2")
  psO_jki_seq2_seq4_Go_2020_T2
psO_jki_seq2_seq4_Go_2020_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_T2) > 0, psO_jki_seq2_seq4_Go_2020_T2)
  psO_jki_seq2_seq4_Go_2020_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_T2_filt))
psO_jki_seq2_seq4_Go_2020_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_T2_filt_annotation
psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Go_2021_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year == "Y2021" & Stage == "T2")
  psO_jki_seq2_seq4_Go_2021_T2
psO_jki_seq2_seq4_Go_2021_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_T2) > 0, psO_jki_seq2_seq4_Go_2021_T2)
  psO_jki_seq2_seq4_Go_2021_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_T2_filt))
psO_jki_seq2_seq4_Go_2021_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_T2_filt_annotation
psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2020_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year == "Y2020" & Stage == "T2")
  psO_jki_seq2_seq4_Ki_2020_T2
psO_jki_seq2_seq4_Ki_2020_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_T2) > 0, psO_jki_seq2_seq4_Ki_2020_T2)
  psO_jki_seq2_seq4_Ki_2020_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_T2_filt))
psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2021_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year == "Y2021" & Stage == "T2")
  psO_jki_seq2_seq4_Ki_2021_T2
psO_jki_seq2_seq4_Ki_2021_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_T2) > 0, psO_jki_seq2_seq4_Ki_2021_T2)
  psO_jki_seq2_seq4_Ki_2021_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_T2_filt))
psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_T2_filt_annotation_rel.csv")
  
  
#######################################
#### Separated by Site and Year and Microhabitat and Stage (T2)
#### Usage = DESeq2
#######################################

## Go 
psO_jki_seq2_seq4_Go_2020_BS_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020" & Stage=="T2" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Go_2020_BS_T2
psO_jki_seq2_seq4_Go_2020_BS_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_BS_T2) > 0, psO_jki_seq2_seq4_Go_2020_BS_T2)
  psO_jki_seq2_seq4_Go_2020_BS_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt))
psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_BS_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_BS_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation
psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel.csv")

psO_jki_seq2_seq4_Go_2020_RH_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020" & Stage=="T2" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Go_2020_RH_T2
psO_jki_seq2_seq4_Go_2020_RH_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_RH_T2) > 0, psO_jki_seq2_seq4_Go_2020_RH_T2)
  psO_jki_seq2_seq4_Go_2020_RH_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt))
psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_RH_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_RH_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation
psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Go_2020_RP_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020" & Stage=="T2" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Go_2020_RP_T2
psO_jki_seq2_seq4_Go_2020_RP_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_RP_T2) > 0, psO_jki_seq2_seq4_Go_2020_RP_T2)
  psO_jki_seq2_seq4_Go_2020_RP_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt))
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_RP_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation
psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Go_2021_BS_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021" & Stage=="T2" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Go_2021_BS_T2
psO_jki_seq2_seq4_Go_2021_BS_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_BS_T2) > 0, psO_jki_seq2_seq4_Go_2021_BS_T2)
  psO_jki_seq2_seq4_Go_2021_BS_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_BS_T2_filt))
psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_BS_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_BS_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation
psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_BS_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Go_2021_RH_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021" & Stage=="T2" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Go_2021_RH_T2
psO_jki_seq2_seq4_Go_2021_RH_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_RH_T2) > 0, psO_jki_seq2_seq4_Go_2021_RH_T2)
  psO_jki_seq2_seq4_Go_2021_RH_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_RH_T2_filt))
psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_RH_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_RH_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation
psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_RH_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Go_2021_RP_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021" & Stage=="T2" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Go_2021_RP_T2
psO_jki_seq2_seq4_Go_2021_RP_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_RP_T2) > 0, psO_jki_seq2_seq4_Go_2021_RP_T2)
  psO_jki_seq2_seq4_Go_2021_RP_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_RP_T2_filt))
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_RP_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation
psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_RP_T2_filt_annotation_rel.csv")
  
### Ki
psO_jki_seq2_seq4_Ki_2020_BS_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020" & Stage=="T2" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Ki_2020_BS_T2
psO_jki_seq2_seq4_Ki_2020_BS_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_BS_T2) > 0, psO_jki_seq2_seq4_Ki_2020_BS_T2)
  psO_jki_seq2_seq4_Ki_2020_BS_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt))
psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2020_RH_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020" & Stage=="T2" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Ki_2020_RH_T2
psO_jki_seq2_seq4_Ki_2020_RH_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_RH_T2) > 0, psO_jki_seq2_seq4_Ki_2020_RH_T2)
  psO_jki_seq2_seq4_Ki_2020_RH_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt))
psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2020_RP_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020" & Stage=="T2" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Ki_2020_RP_T2
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_RP_T2) > 0, psO_jki_seq2_seq4_Ki_2020_RP_T2)
  psO_jki_seq2_seq4_Ki_2020_RP_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt))
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2021_BS_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021" & Stage=="T2" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Ki_2021_BS_T2
psO_jki_seq2_seq4_Ki_2021_BS_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_BS_T2) > 0, psO_jki_seq2_seq4_Ki_2021_BS_T2)
  psO_jki_seq2_seq4_Ki_2021_BS_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt))
psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_BS_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2021_RH_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021" & Stage=="T2" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Ki_2021_RH_T2
psO_jki_seq2_seq4_Ki_2021_RH_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_RH_T2) > 0, psO_jki_seq2_seq4_Ki_2021_RH_T2)
  psO_jki_seq2_seq4_Ki_2021_RH_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt))
psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation
psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_RH_T2_filt_annotation_rel.csv")
  
psO_jki_seq2_seq4_Ki_2021_RP_T2 <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021" & Stage=="T2" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Ki_2021_RP_T2
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_RP_T2) > 0, psO_jki_seq2_seq4_Ki_2021_RP_T2)
  psO_jki_seq2_seq4_Ki_2021_RP_T2_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt))
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation)
psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel))
df_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_RP_T2_filt_annotation_rel.csv")

  
  
  
  
  #######################################
  #### Separated by Site and Year and Microhabitat and Stage (T2)
  #### Usage = DESeq2
  #######################################
  
  ## Go 
  psO_jki_seq2_seq4_Go_2020_BS <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Go_2020_BS
  psO_jki_seq2_seq4_Go_2020_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_BS) > 0, psO_jki_seq2_seq4_Go_2020_BS)
  psO_jki_seq2_seq4_Go_2020_BS_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_BS_filt))
  psO_jki_seq2_seq4_Go_2020_BS_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_BS_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_BS_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_BS_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_BS_filt_annotation
  psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_BS_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_BS_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Go_2020_RH <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Go_2020_RH
  psO_jki_seq2_seq4_Go_2020_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_RH) > 0, psO_jki_seq2_seq4_Go_2020_RH)
  psO_jki_seq2_seq4_Go_2020_RH_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_RH_filt))
  psO_jki_seq2_seq4_Go_2020_RH_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_RH_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_RH_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_RH_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_RH_filt_annotation
  psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_RH_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_RH_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Go_2020_RP <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2020" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Go_2020_RP
  psO_jki_seq2_seq4_Go_2020_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2020_RP) > 0, psO_jki_seq2_seq4_Go_2020_RP)
  psO_jki_seq2_seq4_Go_2020_RP_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2020_RP_filt))
  psO_jki_seq2_seq4_Go_2020_RP_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2020_RP_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2020_RP_filt); ntaxa(psO_jki_seq2_seq4_Go_2020_RP_filt_annotation)
  psO_jki_seq2_seq4_Go_2020_RP_filt_annotation
  psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2020_RP_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2020_RP_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Go_2021_BS <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Go_2021_BS
  psO_jki_seq2_seq4_Go_2021_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_BS) > 0, psO_jki_seq2_seq4_Go_2021_BS)
  psO_jki_seq2_seq4_Go_2021_BS_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_BS_filt))
  psO_jki_seq2_seq4_Go_2021_BS_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_BS_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_BS_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_BS_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_BS_filt_annotation
  psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_BS_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_BS_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Go_2021_RH <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Go_2021_RH
  psO_jki_seq2_seq4_Go_2021_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_RH) > 0, psO_jki_seq2_seq4_Go_2021_RH)
  psO_jki_seq2_seq4_Go_2021_RH_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_RH_filt))
  psO_jki_seq2_seq4_Go_2021_RH_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_RH_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_RH_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_RH_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_RH_filt_annotation
  psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_RH_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_RH_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Go_2021_RP <- subset_samples(psO_jki_seq2_seq4, Site == "Go" & Year=="Y2021" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Go_2021_RP
  psO_jki_seq2_seq4_Go_2021_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Go_2021_RP) > 0, psO_jki_seq2_seq4_Go_2021_RP)
  psO_jki_seq2_seq4_Go_2021_RP_filt
  colnames(tax_table(psO_jki_seq2_seq4_Go_2021_RP_filt))
  psO_jki_seq2_seq4_Go_2021_RP_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Go_2021_RP_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Go_2021_RP_filt); ntaxa(psO_jki_seq2_seq4_Go_2021_RP_filt_annotation)
  psO_jki_seq2_seq4_Go_2021_RP_filt_annotation
  psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Go_2021_RP_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Go_2021_RP_filt_annotation_rel.csv")
  
  ### Ki
  psO_jki_seq2_seq4_Ki_2020_BS <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Ki_2020_BS
  psO_jki_seq2_seq4_Ki_2020_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_BS) > 0, psO_jki_seq2_seq4_Ki_2020_BS)
  psO_jki_seq2_seq4_Ki_2020_BS_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_BS_filt))
  psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_BS_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_BS_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation
  psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_BS_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Ki_2020_RH <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Ki_2020_RH
  psO_jki_seq2_seq4_Ki_2020_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_RH) > 0, psO_jki_seq2_seq4_Ki_2020_RH)
  psO_jki_seq2_seq4_Ki_2020_RH_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_RH_filt))
  psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_RH_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_RH_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation
  psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_RH_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Ki_2020_RP <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2020" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Ki_2020_RP
  psO_jki_seq2_seq4_Ki_2020_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2020_RP) > 0, psO_jki_seq2_seq4_Ki_2020_RP)
  psO_jki_seq2_seq4_Ki_2020_RP_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2020_RP_filt))
  psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2020_RP_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2020_RP_filt); ntaxa(psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation)
  psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation
  psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2020_RP_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Ki_2021_BS <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021" & Microhabitat=="BS")
  psO_jki_seq2_seq4_Ki_2021_BS
  psO_jki_seq2_seq4_Ki_2021_BS_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_BS) > 0, psO_jki_seq2_seq4_Ki_2021_BS)
  psO_jki_seq2_seq4_Ki_2021_BS_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_BS_filt))
  psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_BS_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_BS_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation
  psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_BS_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Ki_2021_RH <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021" & Microhabitat=="RH")
  psO_jki_seq2_seq4_Ki_2021_RH
  psO_jki_seq2_seq4_Ki_2021_RH_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_RH) > 0, psO_jki_seq2_seq4_Ki_2021_RH)
  psO_jki_seq2_seq4_Ki_2021_RH_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_RH_filt))
  psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_RH_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_RH_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation
  psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_RH_filt_annotation_rel.csv")
  
  psO_jki_seq2_seq4_Ki_2021_RP <- subset_samples(psO_jki_seq2_seq4, Site == "Ki" & Year=="Y2021" & Microhabitat=="RP")
  psO_jki_seq2_seq4_Ki_2021_RP
  psO_jki_seq2_seq4_Ki_2021_RP_filt<-prune_taxa(taxa_sums(psO_jki_seq2_seq4_Ki_2021_RP) > 0, psO_jki_seq2_seq4_Ki_2021_RP)
  psO_jki_seq2_seq4_Ki_2021_RP_filt
  colnames(tax_table(psO_jki_seq2_seq4_Ki_2021_RP_filt))
  psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation <- tax_glom(psO_jki_seq2_seq4_Ki_2021_RP_filt, taxrank = "Annotation")
  ntaxa(psO_jki_seq2_seq4_Ki_2021_RP_filt); ntaxa(psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation)
  psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel<-transform_sample_counts(psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation, function(x) (x*100)/sum(x))
  psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel
  head(otu_table(psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel))
  df_psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel <- data.frame(tax_table(psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel),otu_table(psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel))
  write.csv(df_psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/df_psO_jki_seq2_seq4_Ki_2021_RP_filt_annotation_rel.csv")
  
  
  
  