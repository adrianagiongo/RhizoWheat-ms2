##Define taxa differential abundant using DESeq2

#Loading package (to install microbiomeSeq use library("devtools")
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")


##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation)$Rotation,10)
head(sample_data(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation)$Rotation,10)
head(sample_data(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation)$Rotation,10)

head(sample_data(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation)$Rotation,10)
head(sample_data(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation)$Rotation,10)
head(sample_data(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation)$Rotation,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
#by year
# diagdds_BS_W1W22020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation, ~ Rotation)
# diagdds_BS_W1WM2020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation, ~ Rotation)
# diagdds_BS_W2WM2020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation, ~ Rotation)
 
# diagdds_RH_W1W22020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation, ~ Rotation)
# diagdds_RH_W1WM2020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation, ~ Rotation)
# diagdds_RH_W2WM2020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation, ~ Rotation)

diagdds_RP_W1W22020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation, ~ Rotation)
diagdds_RP_W1WM2020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation, ~ Rotation)
diagdds_RP_W2WM2020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation, ~ Rotation)

# diagdds_BS_W1W32020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation, ~ Rotation)
 
# diagdds_RH_W1W32020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation, ~ Rotation)

diagdds_RP_W1W32020 = phyloseq_to_deseq2(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation, ~ Rotation)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
#by year
# diagdds_BS_rotW1W22020 = DESeq(diagdds_BS_W1W22020, test="Wald", fitType="parametric")
# diagdds_BS_rotW1W22020$Rotation <- relevel(diagdds_BS_rotW1W22020$Rotation, ref = "W1")
# 
# diagdds_BS_rotW1WM2020 = DESeq(diagdds_BS_W1WM2020, test="Wald", fitType="parametric")
# diagdds_BS_rotW1WM2020$Rotation <- relevel(diagdds_BS_rotW1WM2020$Rotation, ref = "W1")
# 
# diagdds_BS_rotW2WM2020 = DESeq(diagdds_BS_W2WM2020, test="Wald", fitType="parametric")
# diagdds_BS_rotW2WM2020$Rotation <- relevel(diagdds_BS_rotW2WM2020$Rotation, ref = "W2")
# 
# diagdds_RH_rotW1W22020 = DESeq(diagdds_RH_W1W22020, test="Wald", fitType="parametric")
# diagdds_RH_rotW1W22020$Rotation <- relevel(diagdds_RH_rotW1W22020$Rotation, ref = "W1")
# 
# diagdds_RH_rotW1WM2020 = DESeq(diagdds_RH_W1WM2020, test="Wald", fitType="parametric")
# diagdds_RH_rotW1WM2020$Rotation <- relevel(diagdds_RH_rotW1WM2020$Rotation, ref = "W1")
# 
# diagdds_RH_rotW2WM2020 = DESeq(diagdds_RH_W2WM2020, test="Wald", fitType="parametric")
# diagdds_RH_rotW2WM2020$Rotation <- relevel(diagdds_RH_rotW2WM2020$Rotation, ref = "W2")

diagdds_RP_rotW1W22020 = DESeq(diagdds_RP_W1W22020, test="Wald", fitType="parametric")
#diagdds_RP_rotW1W22020$Rotation <- relevel(diagdds_RP_rotW1W22020$Rotation, ref = "W1")

diagdds_RP_rotW1WM2020 = DESeq(diagdds_RP_W1WM2020, test="Wald", fitType="parametric")
#diagdds_RP_rotW1WM2020$Rotation <- relevel(diagdds_RP_rotW1WM2020$Rotation, ref = "W1")

diagdds_RP_rotW2WM2020 = DESeq(diagdds_RP_W2WM2020, test="Wald", fitType="parametric")
#diagdds_RP_rotW2WM2020$Rotation <- relevel(diagdds_RP_rotW2WM2020$Rotation, ref = "W2")

# diagdds_BS_rotW1W32020 = DESeq(diagdds_BS_W1W32020, test="Wald", fitType="parametric")
# diagdds_BS_rotW1W32020$Rotation <- relevel(diagdds_BS_rotW1W32020$Rotation, ref = "W1")

#samples are too dispersed for parametric fitType
# diagdds_RH_rotW1W32020 = DESeq(diagdds_RH_W1W32020, test="Wald", fitType="parametric")
# diagdds_RH_rotW1W32020$Rotation <- relevel(diagdds_RH_rotW1W32020$Rotation, ref = "W1")
# 
diagdds_RP_rotW1W32020 = DESeq(diagdds_RP_W1W32020, test="Wald", fitType="parametric")
#diagdds_RP_rotW1W32020$Rotation <- relevel(diagdds_RP_rotW1W32020$Rotation, ref = "W1")


##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
#rotW1W22020
# res_BS_rotW1W22020 = results(diagdds_BS_rotW1W22020, contrast=c("Rotation","W2","W1"))
# summary(res_BS_rotW1W22020)
# 
# res_RH_rotW1W22020 = results(diagdds_RH_rotW1W22020, contrast=c("Rotation","W2","W1"))
# summary(res_RH_rotW1W22020)

res_RP_rotW1W22020 = results(diagdds_RP_rotW1W22020, contrast=c("Rotation","W2","W1"))
summary(res_RP_rotW1W22020)

#rotW1WM2020
# res_BS_rotW1WM2020 = results(diagdds_BS_rotW1WM2020, contrast=c("Rotation","WM","W1"))
# summary(res_BS_rotW1WM2020)
# 
# res_RH_rotW1WM2020 = results(diagdds_RH_rotW1WM2020, contrast=c("Rotation","WM","W1"))
# summary(res_RH_rotW1WM2020)

res_RP_rotW1WM2020 = results(diagdds_RP_rotW1WM2020, contrast=c("Rotation","WM","W1"))
summary(res_RP_rotW1WM2020)

# rotW2WM2020
# res_BS_rotW2WM2020 = results(diagdds_BS_rotW2WM2020, contrast=c("Rotation","WM","W2"))
# summary(res_BS_rotW2WM2020)

# res_RH_rotW2WM2020 = results(diagdds_RH_rotW2WM2020, contrast=c("Rotation","WM","W2"))
# summary(res_RH_rotW2WM2020)
# 
# res_RP_rotW2WM2020 = results(diagdds_RP_rotW2WM2020, contrast=c("Rotation","WM","W2"))
# summary(res_RP_rotW2WM2020)


#rotW1W32020
# res_BS_rotW1W32020 = results(diagdds_BS_rotW1W32020, contrast=c("Rotation","W3","W1"))
# summary(res_BS_rotW1W32020)
# 
# res_RH_rotW1W32020 = results(diagdds_RH_rotW1W32020, contrast=c("Rotation","W3","W1"))
# summary(res_RH_rotW1W32020)

res_RP_rotW1W32020 = results(diagdds_RP_rotW1W32020, contrast=c("Rotation","W3","W1"))
summary(res_RP_rotW1W32020)

##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

#rotW1W22020
alpha_DESeq2 = 0.05
# sigtab_BS_rotW1W22020 = res_BS_rotW1W22020[which(res_BS_rotW1W22020$padj < alpha_DESeq2), ]
# 
# sigtab_RH_rotW1W22020 = res_RH_rotW1W22020[which(res_RH_rotW1W22020$padj < alpha_DESeq2), ]

sigtab_RP_rotW1W22020 = res_RP_rotW1W22020[which(res_RP_rotW1W22020$padj < alpha_DESeq2), ]

#Formating results on table and save
# sigtab_BS_rotW1W22020 = cbind(as(sigtab_BS_rotW1W22020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation)[rownames(sigtab_BS_rotW1W22020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel)[rownames(sigtab_BS_rotW1W22020), ], "matrix"))
# head(sigtab_BS_rotW1W22020)
# dim(sigtab_BS_rotW1W22020)
# 
# sigtab_RH_rotW1W22020 = cbind(as(sigtab_RH_rotW1W22020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation)[rownames(sigtab_RH_rotW1W22020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel)[rownames(sigtab_RH_rotW1W22020), ], "matrix"))
# head(sigtab_RH_rotW1W22020)
# dim(sigtab_RH_rotW1W22020)

sigtab_RP_rotW1W22020 = cbind(as(sigtab_RP_rotW1W22020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation)[rownames(sigtab_RP_rotW1W22020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel)[rownames(sigtab_RP_rotW1W22020), ], "matrix"))
head(sigtab_RP_rotW1W22020)
dim(sigtab_RP_rotW1W22020)

# write.csv(sigtab_BS_rotW1W22020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_BS_rotW1W22020.csv")
# write.csv(sigtab_RH_rotW1W22020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_RH_rotW1W22020.csv")
write.csv(sigtab_RP_rotW1W22020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_RP_rotW1W22020.csv")

#rotW1WM2020
alpha_DESeq2 = 0.05
# sigtab_BS_rotW1WM2020 = res_BS_rotW1WM2020[which(res_BS_rotW1WM2020$padj < alpha_DESeq2), ]
# 
# sigtab_RH_rotW1WM2020 = res_RH_rotW1WM2020[which(res_RH_rotW1WM2020$padj < alpha_DESeq2), ]

sigtab_RP_rotW1WM2020 = res_RP_rotW1WM2020[which(res_RP_rotW1WM2020$padj < alpha_DESeq2), ]

#Formating results on table and save
# sigtab_BS_rotW1WM2020 = cbind(as(sigtab_BS_rotW1WM2020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation)[rownames(sigtab_BS_rotW1WM2020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel)[rownames(sigtab_BS_rotW1WM2020), ], "matrix"))
# head(sigtab_BS_rotW1WM2020)
# dim(sigtab_BS_rotW1WM2020)
# 
# sigtab_RH_rotW1WM2020 = cbind(as(sigtab_RH_rotW1WM2020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation)[rownames(sigtab_RH_rotW1WM2020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel)[rownames(sigtab_RH_rotW1WM2020), ], "matrix"))
# head(sigtab_RH_rotW1WM2020)
# dim(sigtab_RH_rotW1WM2020)

sigtab_RP_rotW1WM2020 = cbind(as(sigtab_RP_rotW1WM2020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation)[rownames(sigtab_RP_rotW1WM2020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel)[rownames(sigtab_RP_rotW1WM2020), ], "matrix"))
head(sigtab_RP_rotW1WM2020)
dim(sigtab_RP_rotW1WM2020)

# write.csv(sigtab_BS_rotW1WM2020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_BS_rotW1WM2020.csv")
# write.csv(sigtab_RH_rotW1WM2020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_RH_rotW1WM2020.csv")
write.csv(sigtab_RP_rotW1WM2020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_RP_rotW1WM2020.csv")

#rotW2WM2020
alpha_DESeq2 = 0.05
# sigtab_BS_rotW2WM2020 = res_BS_rotW2WM2020[which(res_BS_rotW2WM2020$padj < alpha_DESeq2), ]
# 
# sigtab_RH_rotW2WM2020 = res_RH_rotW2WM2020[which(res_RH_rotW2WM2020$padj < alpha_DESeq2), ]
# 
# sigtab_RP_rotW2WM2020 = res_RP_rotW2WM2020[which(res_RP_rotW2WM2020$padj < alpha_DESeq2), ]
# 
# #Formating results on table and save
# sigtab_BS_rotW2WM2020 = cbind(as(sigtab_BS_rotW2WM2020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation)[rownames(sigtab_BS_rotW2WM2020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation_rel)[rownames(sigtab_BS_rotW2WM2020), ], "matrix"))
# head(sigtab_BS_rotW2WM2020)
# dim(sigtab_BS_rotW2WM2020)
# 
# sigtab_RH_rotW2WM2020 = cbind(as(sigtab_RH_rotW2WM2020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation)[rownames(sigtab_RH_rotW2WM2020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation_rel)[rownames(sigtab_RH_rotW2WM2020), ], "matrix"))
# head(sigtab_RH_rotW2WM2020)
# dim(sigtab_RH_rotW2WM2020)
# 
# sigtab_RP_rotW2WM2020 = cbind(as(sigtab_RP_rotW2WM2020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation)[rownames(sigtab_RP_rotW2WM2020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation_rel)[rownames(sigtab_RP_rotW2WM2020), ], "matrix"))
# head(sigtab_RP_rotW2WM2020)
# dim(sigtab_RP_rotW2WM2020)
# 
# write.csv(sigtab_BS_rotW2WM2020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_rotW2WM2020_deseq2_psO_jki_seq2_seq4_Go_2020_BS_T2_filt_annotation.csv")
# write.csv(sigtab_RH_rotW2WM2020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_rotW2WM2020_deseq2_psO_jki_seq2_seq4_Go_2020_RH_T2_filt_annotation.csv")
# write.csv(sigtab_RP_rotW2WM2020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_rotW2WM2020_deseq2_psO_jki_seq2_seq4_Go_2020_RP_T2_filt_annotation.csv")


#rotW1W3_2020
alpha_DESeq2 = 0.05
# sigtab_BS_rotW1W32020 = res_BS_rotW1W32020[which(res_BS_rotW1W32020$padj < alpha_DESeq2), ]
# 
# sigtab_RH_rotW1W32020 = res_RH_rotW1W32020[which(res_RH_rotW1W32020$padj < alpha_DESeq2), ]
# 
sigtab_RP_rotW1W32020 = res_RP_rotW1W32020[which(res_RP_rotW1W32020$padj < alpha_DESeq2), ]

#Formating results on table and save
# sigtab_BS_rotW1W32020 = cbind(as(sigtab_BS_rotW1W32020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation)[rownames(sigtab_BS_rotW1W32020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Ki_2020_BS_T2_filt_annotation_rel)[rownames(sigtab_BS_rotW1W32020), ], "matrix"))
# head(sigtab_BS_rotW1W32020)
# dim(sigtab_BS_rotW1W32020)
# 
# sigtab_RH_rotW1W32020 = cbind(as(sigtab_RH_rotW1W32020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation)[rownames(sigtab_RH_rotW1W32020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Ki_2020_RH_T2_filt_annotation_rel)[rownames(sigtab_RH_rotW1W32020), ], "matrix"))
# head(sigtab_RH_rotW1W32020)
# dim(sigtab_RH_rotW1W32020)

sigtab_RP_rotW1W32020 = cbind(as(sigtab_RP_rotW1W32020, "data.frame"), as(tax_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation)[rownames(sigtab_RP_rotW1W32020), ], "matrix"), as(otu_table(psO_jki_seq2_seq4_Ki_2020_RP_T2_filt_annotation_rel)[rownames(sigtab_RP_rotW1W32020), ], "matrix"))
head(sigtab_RP_rotW1W32020)
dim(sigtab_RP_rotW1W32020)

# write.csv(sigtab_BS_rotW1W32020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_BS_rotW1W32020.csv")
# write.csv(sigtab_RH_rotW1W32020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_RH_rotW1W32020.csv")
write.csv(sigtab_RP_rotW1W32020, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_RP_rotW1W32020.csv")



##The end!

