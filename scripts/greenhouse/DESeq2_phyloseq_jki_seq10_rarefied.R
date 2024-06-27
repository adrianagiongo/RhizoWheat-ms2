##Define taxa differential abundant using DESeq2

#Loading package (to install microbiomeSeq use library("devtools")
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")


##Differential analysis
#Summarize the first few entries (10) of the factor
head(sample_data(psO_jki_seq10_rarefied_Go_W1_filt_annotation)$Treatment,10)
head(sample_data(psO_jki_seq10_rarefied_Go_WM_filt_annotation)$Treatment,10)

head(sample_data(psO_jki_seq10_rarefied_Ki_W1_filt_annotation)$Treatment,10)
head(sample_data(psO_jki_seq10_rarefied_Ki_W3_filt_annotation)$Treatment,10)

diagdds_rarefied_Go_W1 = phyloseq_to_deseq2(psO_jki_seq10_rarefied_Go_W1_filt_annotation, ~ Treatment)
diagdds_rarefied_Go_WM = phyloseq_to_deseq2(psO_jki_seq10_rarefied_Go_WM_filt_annotation, ~ Treatment)

diagdds_rarefied_Ki_W1 = phyloseq_to_deseq2(psO_jki_seq10_rarefied_Ki_W1_filt_annotation, ~ Treatment)
diagdds_rarefied_Ki_W3 = phyloseq_to_deseq2(psO_jki_seq10_rarefied_Ki_W3_filt_annotation, ~ Treatment)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
diagdds_rarefied_Go_W1 = DESeq(diagdds_rarefied_Go_W1, test="Wald", fitType="parametric")
diagdds_rarefied_Go_WM = DESeq(diagdds_rarefied_Go_WM, test="Wald", fitType="parametric")
diagdds_rarefied_Ki_W1 = DESeq(diagdds_rarefied_Ki_W1, test="Wald", fitType="parametric")
diagdds_rarefied_Ki_W3 = DESeq(diagdds_rarefied_Ki_W3, test="Wald", fitType="parametric")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
res_rarefied_Go_W1 = results(diagdds_rarefied_Go_W1, contrast=c("Treatment","Ggt","Control"))
summary(res_rarefied_Go_W1)

res_rarefied_Go_WM = results(diagdds_rarefied_Go_WM, contrast=c("Treatment","Ggt","Control"))
summary(res_rarefied_Go_WM)

res_rarefied_Ki_W1 = results(diagdds_rarefied_Ki_W1, contrast=c("Treatment","Ggt","Control"))
summary(res_rarefied_Ki_W1)

res_rarefied_Ki_W3 = results(diagdds_rarefied_Ki_W3, contrast=c("Treatment","Ggt","Control"))
summary(res_rarefied_Ki_W3)


##############################
#for Phyloseq plot and tables
##############################
#Indicate alpha and adjust values

#rotW1W2
alpha_DESeq2 = 0.05

### Go_W1
sigtab_rarefied_Go_W1 = res_rarefied_Go_W1[which(res_rarefied_Go_W1$padj < alpha_DESeq2), ]
sigtab_rarefied_Go_WM = res_rarefied_Go_WM[which(res_rarefied_Go_WM$padj < alpha_DESeq2), ]

sigtab_rarefied_Ki_W1 = res_rarefied_Ki_W1[which(res_rarefied_Ki_W1$padj < alpha_DESeq2), ]
sigtab_rarefied_Ki_W3 = res_rarefied_Ki_W3[which(res_rarefied_Ki_W3$padj < alpha_DESeq2), ]

#Formating results on table and save
sigtab_rarefied_Go_W1 = cbind(as(sigtab_rarefied_Go_W1, "data.frame"), as(tax_table(psO_jki_seq10_rarefied_Go_W1_filt_annotation)[rownames(sigtab_rarefied_Go_W1), ], "matrix"), as(otu_table(psO_jki_seq10_rarefied_Go_W1_filt_annotation_rel)[rownames(sigtab_rarefied_Go_W1), ], "matrix"))
head(sigtab_rarefied_Go_W1)
dim(sigtab_rarefied_Go_W1)

write.csv(sigtab_rarefied_Go_W1, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/sigtab_rarefied_Go_W1.csv")


#Formating results on table and save
sigtab_rarefied_Go_WM = cbind(as(sigtab_rarefied_Go_WM, "data.frame"), as(tax_table(psO_jki_seq10_rarefied_Go_WM_filt_annotation)[rownames(sigtab_rarefied_Go_WM), ], "matrix"), as(otu_table(psO_jki_seq10_rarefied_Go_WM_filt_annotation_rel)[rownames(sigtab_rarefied_Go_WM), ], "matrix"))
head(sigtab_rarefied_Go_WM)
dim(sigtab_rarefied_Go_WM)

write.csv(sigtab_rarefied_Go_WM, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/sigtab_rarefied_Go_WM.csv")


#Formating results on table and save
sigtab_rarefied_Ki_W1 = cbind(as(sigtab_rarefied_Ki_W1, "data.frame"), as(tax_table(psO_jki_seq10_rarefied_Ki_W1_filt_annotation)[rownames(sigtab_rarefied_Ki_W1), ], "matrix"), as(otu_table(psO_jki_seq10_rarefied_Ki_W1_filt_annotation_rel)[rownames(sigtab_rarefied_Ki_W1), ], "matrix"))
head(sigtab_rarefied_Ki_W1)
dim(sigtab_rarefied_Ki_W1)

write.csv(sigtab_rarefied_Ki_W1, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/sigtab_rarefied_Ki_W1.csv")


#Formating results on table and save
sigtab_rarefied_Ki_W3 = cbind(as(sigtab_rarefied_Ki_W3, "data.frame"), as(tax_table(psO_jki_seq10_rarefied_Ki_W3_filt_annotation)[rownames(sigtab_rarefied_Ki_W3), ], "matrix"), as(otu_table(psO_jki_seq10_rarefied_Ki_W3_filt_annotation_rel)[rownames(sigtab_rarefied_Ki_W3), ], "matrix"))
head(sigtab_rarefied_Ki_W3)
dim(sigtab_rarefied_Ki_W3)

write.csv(sigtab_rarefied_Ki_W3, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/sigtab_rarefied_Ki_W3.csv")


##The end!

