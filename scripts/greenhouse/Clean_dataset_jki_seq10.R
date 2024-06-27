##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

####################################### 
#Show available ranks in the dataset 
rank_names(jki_seq10_data)
sample_names(jki_seq10_data)
jki_seq10_data

#Removing 2 bad quality samples (total 32 - 2 = 30 samples left) 
jki_seq10_data_30 <- subset_samples(jki_seq10_data, !(Sample_id)%in%c("GH06_89", "GH06_90"))
jki_seq10_data_30

#Filter
jki_seq10_data_30_filt <- prune_taxa(taxa_sums(jki_seq10_data_30) > 0, jki_seq10_data_30)
jki_seq10_data_30_filt

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level (Keeping Archaea)
#Kingdom
table(tax_table(jki_seq10_data_30_filt)[,"Kingdom"], exclude = NULL)
psO_jki_seq10 <- subset_taxa(jki_seq10_data_30_filt, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota"))
table(tax_table(psO_jki_seq10)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_jki_seq10)[,"Phylum"], exclude = NULL)
psO_jki_seq10 <- subset_taxa(psO_jki_seq10, !Phylum %in% c("Archaea")) #No NAs anymore (Archaea instead)
table(tax_table(psO_jki_seq10)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_jki_seq10)[,"Order"], exclude = NULL)
psO_jki_seq10 <- subset_taxa(psO_jki_seq10, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(psO_jki_seq10)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_jki_seq10)[,"Family"], exclude = NULL)
psO_jki_seq10 <- subset_taxa(psO_jki_seq10, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(psO_jki_seq10)[,"Family"], exclude = NULL)

#Genus
#table(tax_table(psO_jki_seq10)[,"Genus"], exclude = NULL)

#Annotation
#table(tax_table(psO_jki_seq10)[,"Annotation"], exclude = NULL)

##Observe psO of selected after clean spurious taxa
psO_jki_seq10

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(psO_jki_seq10),
               MARGIN = ifelse(taxa_are_rows(psO_jki_seq10), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_jki_seq10 = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(psO_jki_seq10), tax_table(psO_jki_seq10), otu_table(psO_jki_seq10))

#Write data frame in csv format
write.csv(prevdf_jki_seq10, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/prevdf_psO_jki_seq10.csv")

## The end : )

