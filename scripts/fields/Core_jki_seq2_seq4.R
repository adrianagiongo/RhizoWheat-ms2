## Core - jki_seq2_seq4
#### Adriana Giongo
#### (21.05.2023) 

#Load packages
library("phyloseq")
library("microbiome")
library("dplyr")

# Calculate compositional version of the data (relative abundances)
#collapse the taxa into "Annotation" level
rank_names(psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt)
aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt, "Annotation")


rank_names(psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt)
aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt, "Annotation")




##Compute prevalence of each feature
prev_Go_2020_BS_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2020_BS_T2_W2 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2020_BS_T2_WM = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2020_BS_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2020_BS_T2_W3 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})


prev_Go_2021_BS_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2021_BS_T2_W2 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2021_BS_T2_WM = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2021_BS_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2021_BS_T2_W3 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})


#Add taxonomy and total reads counts to a data frame
prev_Go_2020_BS_T2_W1 = data.frame(Prevalence=prev_Go_2020_BS_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation))
prev_Go_2020_BS_T2_W2 = data.frame(Prevalence=prev_Go_2020_BS_T2_W2, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation))
prev_Go_2020_BS_T2_WM = data.frame(Prevalence=prev_Go_2020_BS_T2_WM, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation))

prev_Ki_2020_BS_T2_W1 = data.frame(Prevalence=prev_Ki_2020_BS_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation))
prev_Ki_2020_BS_T2_W3 = data.frame(Prevalence=prev_Ki_2020_BS_T2_W3, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation))


prev_Go_2021_BS_T2_W1 = data.frame(Prevalence=prev_Go_2021_BS_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation))
prev_Go_2021_BS_T2_W2 = data.frame(Prevalence=prev_Go_2021_BS_T2_W2, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation))
prev_Go_2021_BS_T2_WM = data.frame(Prevalence=prev_Go_2021_BS_T2_WM, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation))

prev_Ki_2021_BS_T2_W1 = data.frame(Prevalence=prev_Ki_2021_BS_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation))
prev_Ki_2021_BS_T2_W3 = data.frame(Prevalence=prev_Ki_2021_BS_T2_W3, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation))

#Write data frame in csv format
write.csv(prev_Go_2020_BS_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation.csv")
write.csv(prev_Go_2020_BS_T2_W2, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation.csv")
write.csv(prev_Go_2020_BS_T2_WM, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation.csv")

write.csv(prev_Ki_2020_BS_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation.csv")
write.csv(prev_Ki_2020_BS_T2_W3, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation.csv")

write.csv(prev_Go_2021_BS_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation.csv")
write.csv(prev_Go_2021_BS_T2_W2, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation.csv")
write.csv(prev_Go_2021_BS_T2_WM, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation.csv")

write.csv(prev_Ki_2021_BS_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation.csv")
write.csv(prev_Ki_2021_BS_T2_W3, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation.csv")


#taxa that exceed the given prevalence and detection thresholds
core_taxa_Go_2020_BS_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2020_BS_T2_W2_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_W2_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2020_BS_T2_WM_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_BS_T2_WM_filt_annotation, detection = 0, prevalence = 99/100)

core_taxa_Ki_2020_BS_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Ki_2020_BS_T2_W3_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2020_BS_T2_W3_filt_annotation, detection = 0, prevalence = 99/100)


core_taxa_Go_2021_BS_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2021_BS_T2_W2_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_W2_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2021_BS_T2_WM_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_BS_T2_WM_filt_annotation, detection = 0, prevalence = 99/100)

core_taxa_Ki_2021_BS_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Ki_2021_BS_T2_W3_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2021_BS_T2_W3_filt_annotation, detection = 0, prevalence = 99/100)


#write core taxa as csv table
write.csv(core_taxa_Go_2020_BS_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_BS_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Go_2020_BS_T2_W2_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_BS_T2_W2_filt_annotation.csv")
write.csv(core_taxa_Go_2020_BS_T2_WM_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_BS_T2_WM_filt_annotation.csv")

write.csv(core_taxa_Ki_2020_BS_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2020_BS_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Ki_2020_BS_T2_W3_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2020_BS_T2_W3_filt_annotation.csv")

write.csv(core_taxa_Go_2021_BS_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_BS_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Go_2021_BS_T2_W2_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_BS_T2_W2_filt_annotation.csv")
write.csv(core_taxa_Go_2021_BS_T2_WM_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_BS_T2_WM_filt_annotation.csv")

write.csv(core_taxa_Ki_2021_BS_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2021_BS_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Ki_2021_BS_T2_W3_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2021_BS_T2_W3_filt_annotation.csv")

### RH
#collapse the taxa into "Annotation" level
rank_names(psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt)
aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt, "Annotation")


rank_names(psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt)
aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt, "Annotation")


##Compute prevalence of each feature
prev_Go_2020_RH_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2020_RH_T2_W2 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2020_RH_T2_WM = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2020_RH_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2020_RH_T2_W3 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})


prev_Go_2021_RH_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2021_RH_T2_W2 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2021_RH_T2_WM = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2021_RH_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2021_RH_T2_W3 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})


#Add taxonomy and total reads counts to a data frame
prev_Go_2020_RH_T2_W1 = data.frame(Prevalence=prev_Go_2020_RH_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation))
prev_Go_2020_RH_T2_W2 = data.frame(Prevalence=prev_Go_2020_RH_T2_W2, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation))
prev_Go_2020_RH_T2_WM = data.frame(Prevalence=prev_Go_2020_RH_T2_WM, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation))

prev_Ki_2020_RH_T2_W1 = data.frame(Prevalence=prev_Ki_2020_RH_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation))
prev_Ki_2020_RH_T2_W3 = data.frame(Prevalence=prev_Ki_2020_RH_T2_W3, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation))


prev_Go_2021_RH_T2_W1 = data.frame(Prevalence=prev_Go_2021_RH_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation))
prev_Go_2021_RH_T2_W2 = data.frame(Prevalence=prev_Go_2021_RH_T2_W2, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation))
prev_Go_2021_RH_T2_WM = data.frame(Prevalence=prev_Go_2021_RH_T2_WM, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation))

prev_Ki_2021_RH_T2_W1 = data.frame(Prevalence=prev_Ki_2021_RH_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation))
prev_Ki_2021_RH_T2_W3 = data.frame(Prevalence=prev_Ki_2021_RH_T2_W3, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation))

#Write data frame in csv format
write.csv(prev_Go_2020_RH_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation.csv")
write.csv(prev_Go_2020_RH_T2_W2, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation.csv")
write.csv(prev_Go_2020_RH_T2_WM, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation.csv")

write.csv(prev_Ki_2020_RH_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation.csv")
write.csv(prev_Ki_2020_RH_T2_W3, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation.csv")

write.csv(prev_Go_2021_RH_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation.csv")
write.csv(prev_Go_2021_RH_T2_W2, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation.csv")
write.csv(prev_Go_2021_RH_T2_WM, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation.csv")

write.csv(prev_Ki_2021_RH_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation.csv")
write.csv(prev_Ki_2021_RH_T2_W3, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation.csv")


#taxa that exceed the given prevalence and detection thresholds
core_taxa_Go_2020_RH_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2020_RH_T2_W2_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_W2_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2020_RH_T2_WM_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_RH_T2_WM_filt_annotation, detection = 0, prevalence = 99/100)

core_taxa_Ki_2020_RH_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Ki_2020_RH_T2_W3_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2020_RH_T2_W3_filt_annotation, detection = 0, prevalence = 99/100)


core_taxa_Go_2021_RH_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2021_RH_T2_W2_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_W2_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2021_RH_T2_WM_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_RH_T2_WM_filt_annotation, detection = 0, prevalence = 99/100)

core_taxa_Ki_2021_RH_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Ki_2021_RH_T2_W3_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2021_RH_T2_W3_filt_annotation, detection = 0, prevalence = 99/100)

#write core taxa as csv table
write.csv(core_taxa_Go_2020_RH_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_RH_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Go_2020_RH_T2_W2_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_RH_T2_W2_filt_annotation.csv")
write.csv(core_taxa_Go_2020_RH_T2_WM_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_RH_T2_WM_filt_annotation.csv")

write.csv(core_taxa_Ki_2020_RH_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2020_RH_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Ki_2020_RH_T2_W3_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2020_RH_T2_W3_filt_annotation.csv")

write.csv(core_taxa_Go_2021_RH_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_RH_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Go_2021_RH_T2_W2_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_RH_T2_W2_filt_annotation.csv")
write.csv(core_taxa_Go_2021_RH_T2_WM_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_RH_T2_WM_filt_annotation.csv")

write.csv(core_taxa_Ki_2021_RH_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2021_RH_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Ki_2021_RH_T2_W3_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2021_RH_T2_W3_filt_annotation.csv")


##RP
#collapse the taxa into "Annotation" level
rank_names(psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt)
aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt)
aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt, "Annotation")


rank_names(psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt)
aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt)
aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt, "Annotation")

rank_names(psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt)
aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation <- aggregate_taxa(psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt, "Annotation")


##Compute prevalence of each feature
prev_Go_2020_RP_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2020_RP_T2_W2 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2020_RP_T2_WM = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2020_RP_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2020_RP_T2_W3 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})


prev_Go_2021_RP_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2021_RP_T2_W2 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Go_2021_RP_T2_WM = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2021_RP_T2_W1 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})

prev_Ki_2021_RP_T2_W3 = apply(X = otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation),
                           MARGIN = ifelse(taxa_are_rows(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation), yes = 1, no = 2),
                           FUN = function (x) {sum(x>0)})


#Add taxonomy and total reads counts to a data frame
prev_Go_2020_RP_T2_W1 = data.frame(Prevalence=prev_Go_2020_RP_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation))
prev_Go_2020_RP_T2_W2 = data.frame(Prevalence=prev_Go_2020_RP_T2_W2, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation))
prev_Go_2020_RP_T2_WM = data.frame(Prevalence=prev_Go_2020_RP_T2_WM, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation))

prev_Ki_2020_RP_T2_W1 = data.frame(Prevalence=prev_Ki_2020_RP_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation))
prev_Ki_2020_RP_T2_W3 = data.frame(Prevalence=prev_Ki_2020_RP_T2_W3, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation))


prev_Go_2021_RP_T2_W1 = data.frame(Prevalence=prev_Go_2021_RP_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation))
prev_Go_2021_RP_T2_W2 = data.frame(Prevalence=prev_Go_2021_RP_T2_W2, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation))
prev_Go_2021_RP_T2_WM = data.frame(Prevalence=prev_Go_2021_RP_T2_WM, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation))

prev_Ki_2021_RP_T2_W1 = data.frame(Prevalence=prev_Ki_2021_RP_T2_W1, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation))
prev_Ki_2021_RP_T2_W3 = data.frame(Prevalence=prev_Ki_2021_RP_T2_W3, TotalAbundance = taxa_sums(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation), tax_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation), otu_table(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation))

#Write data frame in csv format
write.csv(prev_Go_2020_RP_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation.csv")
write.csv(prev_Go_2020_RP_T2_W2, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation.csv")
write.csv(prev_Go_2020_RP_T2_WM, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation.csv")

write.csv(prev_Ki_2020_RP_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation.csv")
write.csv(prev_Ki_2020_RP_T2_W3, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation.csv")

write.csv(prev_Go_2021_RP_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation.csv")
write.csv(prev_Go_2021_RP_T2_W2, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation.csv")
write.csv(prev_Go_2021_RP_T2_WM, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation.csv")

write.csv(prev_Ki_2021_RP_T2_W1, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation.csv")
write.csv(prev_Ki_2021_RP_T2_W3, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/prev_aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation.csv")


#taxa that exceed the given prevalence and detection thresholds
core_taxa_Go_2020_RP_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2020_RP_T2_W2_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_W2_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2020_RP_T2_WM_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2020_RP_T2_WM_filt_annotation, detection = 0, prevalence = 99/100)

core_taxa_Ki_2020_RP_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Ki_2020_RP_T2_W3_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2020_RP_T2_W3_filt_annotation, detection = 0, prevalence = 99/100)


core_taxa_Go_2021_RP_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2021_RP_T2_W2_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_W2_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Go_2021_RP_T2_WM_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Go_2021_RP_T2_WM_filt_annotation, detection = 0, prevalence = 99/100)

core_taxa_Ki_2021_RP_T2_W1_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W1_filt_annotation, detection = 0, prevalence = 99/100)
core_taxa_Ki_2021_RP_T2_W3_filt_annotation <- core_members(aggreg_psO_jki_seq2_seq4_Ki_2021_RP_T2_W3_filt_annotation, detection = 0, prevalence = 99/100)

#write core taxa as csv table
write.csv(core_taxa_Go_2020_RP_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_RP_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Go_2020_RP_T2_W2_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_RP_T2_W2_filt_annotation.csv")
write.csv(core_taxa_Go_2020_RP_T2_WM_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2020_RP_T2_WM_filt_annotation.csv")

write.csv(core_taxa_Ki_2020_RP_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2020_RP_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Ki_2020_RP_T2_W3_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2020_RP_T2_W3_filt_annotation.csv")

write.csv(core_taxa_Go_2021_RP_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_RP_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Go_2021_RP_T2_W2_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_RP_T2_W2_filt_annotation.csv")
write.csv(core_taxa_Go_2021_RP_T2_WM_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Go_2021_RP_T2_WM_filt_annotation.csv")

write.csv(core_taxa_Ki_2021_RP_T2_W1_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2021_RP_T2_W1_filt_annotation.csv")
write.csv(core_taxa_Ki_2021_RP_T2_W3_filt_annotation, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/core_taxa_Ki_2021_RP_T2_W3_filt_annotation.csv")


## The end! : )


