##Supervised cleaning process to clean the raw dataset

#load packages
library("phyloseq")
library("ggplot2")
library("dplyr")

####################################### 
#Show available ranks in the dataset  --> original dataset = 153 Go + 136 Ki = 289 samples
rank_names(jki_seq2_seq4_data)
sample_names(jki_seq2_seq4_data)
jki_seq2_seq4_data

#Removing 3 bad quality samples + 34 RA from 2020 (2 are the same of bad quality) (total 254 samples left) + 51 samples T1 = 
jki_seq2_seq4_data_286 <- subset_samples(jki_seq2_seq4_data, !(Sample_id)%in%c("Go018", "Ki028", "Ki054"))
jki_seq2_seq4_data_286
jki_seq2_seq4_data_254 <- subset_samples(jki_seq2_seq4_data_286, !(Microhabitat)%in%c("RA"))
jki_seq2_seq4_data_254
jki_seq2_seq4_data_203 <- subset_samples(jki_seq2_seq4_data_254, !(Stage)%in%c("T1"))
jki_seq2_seq4_data_203

#Filter
jki_seq2_seq4_data_203_filt <- prune_taxa(taxa_sums(jki_seq2_seq4_data_203) > 0, jki_seq2_seq4_data_203)
jki_seq2_seq4_data_203_filt

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level (Keeping Archaea)
#Kingdom
table(tax_table(jki_seq2_seq4_data_203_filt)[,"Kingdom"], exclude = NULL)
psO_jki_seq2_seq4 <- subset_taxa(jki_seq2_seq4_data_203_filt, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota"))
table(tax_table(psO_jki_seq2_seq4)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_jki_seq2_seq4)[,"Phylum"], exclude = NULL)
psO_jki_seq2_seq4 <- subset_taxa(psO_jki_seq2_seq4, !Phylum %in% c("Archaea")) #No NAs anymore (Archaea instead)
table(tax_table(psO_jki_seq2_seq4)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_jki_seq2_seq4)[,"Order"], exclude = NULL)
psO_jki_seq2_seq4 <- subset_taxa(psO_jki_seq2_seq4, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(psO_jki_seq2_seq4)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_jki_seq2_seq4)[,"Family"], exclude = NULL)
psO_jki_seq2_seq4 <- subset_taxa(psO_jki_seq2_seq4, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(psO_jki_seq2_seq4)[,"Family"], exclude = NULL)

#Genus
#table(tax_table(psO_jki_seq2_seq4)[,"Genus"], exclude = NULL)

#Annotation
#table(tax_table(psO_jki_seq2_seq4)[,"Annotation"], exclude = NULL)

##Observe psO of selected after clean spurious taxa
jki_seq2_seq4_data_286
psO_jki_seq2_seq4

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(psO_jki_seq2_seq4),
               MARGIN = ifelse(taxa_are_rows(psO_jki_seq2_seq4), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf_jki_seq2_seq4 = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(psO_jki_seq2_seq4), tax_table(psO_jki_seq2_seq4), otu_table(psO_jki_seq2_seq4))

#Write data frame in csv format
write.csv(prevdf_jki_seq2_seq4, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prevdf_psO_jki_seq2_seq4.csv")

## The end : )

