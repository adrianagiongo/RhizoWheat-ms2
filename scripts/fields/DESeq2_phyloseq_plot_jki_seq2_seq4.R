##Define taxa deferentially abundant using DESeq2 on phyloseq pipeline

#Loading package
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")


# # colors_RP <- list(
#     Phylum =  c("Acidobacteriota" = "#ffc331", 
#                 "Actinobacteriota" = "#ff8b60", 
#                 "Bacteroidota"  = "#Bff4be", 
#                 "Campilobacterota" = "#800080",
#                 "Chloroflexi" = "#2b8f22",
#                 "Cyanobacteria" = "#185113",                 
#                 "Deferribacterota" = "#d4ff31", 
#                 "Desulfobacterota" = "#964B00", 
#                 "Firmicutes" = "#4B69F5",
#                 "Gemmatimonadota" = "#363781",
#                 "Latescibacterota" = "#4a2500", 
#                 "Myxococcota" = "#FFC0CB", 
#                 "NB1-j" = "#b5814c", 
#                 "Nitrospirota" = "#1361cf",
#                 "Patescibacteria" = "#4b0096", 
#                 "Planctomycetota" = "#9fc5e8", 
#                 "Proteobacteria" = "#ca250f", 
#                 "Verrucomicrobiota" = "#ced2ce",
#                 "Zixibacteria" = "#e89fc5")




########## W1W2 2020

#BS_rotW1W22020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW1W22020$log2FoldChange, sigtab_BS_rotW1W22020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W22020$Phylum = factor(as.character(sigtab_BS_rotW1W22020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW1W22020$log2FoldChange, sigtab_BS_rotW1W22020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W22020$Annotation = factor(as.character(sigtab_BS_rotW1W22020$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW1W22020 <- c("Acidobacteriota" = "#ffc331",
                                 "Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Campilobacterota" = "#800080",
                                 "Chloroflexi" = "#2b8f22",
                                 "Crenarchaeota" = "#e89fc5",
                                 "Desulfobacterota" = "#964B00",
                                 "Firmicutes" = "#4B69F5",
                                 "Gemmatimonadota" = "#363781",
                                 "Myxococcota" = "#FFC0CB",
                                 "Proteobacteria" = "#ca250f",
                                 "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_BS_rotW1W22020 <- ggplot(sigtab_BS_rotW1W22020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW1W22020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW1W22020

ggsave("plot_sigtab_BS_rotW1W22020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 32, units = "cm", dpi = 300, device = "png")

#RH_rotW1W22020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RH_rotW1W22020$log2FoldChange, sigtab_RH_rotW1W22020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW1W22020$Phylum = factor(as.character(sigtab_RH_rotW1W22020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RH_rotW1W22020$log2FoldChange, sigtab_RH_rotW1W22020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW1W22020$Annotation = factor(as.character(sigtab_RH_rotW1W22020$Annotation), levels=names(x))

# Define color palette
color.phylum_RH_rotW1W22020 <- c("Acidobacteriota" = "#ffc331",
                                 "Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Deferribacterota" = "#d4ff31",
                                 "Desulfobacterota" = "#964B00",
                                 "Firmicutes" = "#4B69F5",
                                 "Myxococcota" = "#FFC0CB",
                                 "Proteobacteria" = "#ca250f",
                                 "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_RH_rotW1W22020 <- ggplot(sigtab_RH_rotW1W22020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RH_rotW1W22020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RH_rotW1W22020

ggsave("plot_sigtab_RH_rotW1W22020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 23, units = "cm", dpi = 300, device = "png")


#RP_rotW1W22020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RP_rotW1W22020$log2FoldChange, sigtab_RP_rotW1W22020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1W22020$Phylum = factor(as.character(sigtab_RP_rotW1W22020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RP_rotW1W22020$log2FoldChange, sigtab_RP_rotW1W22020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1W22020$Annotation = factor(as.character(sigtab_RP_rotW1W22020$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1W22020 <- c("Actinobacteriota" = "#ff8b60",
                                 "Firmicutes" = "#4B69F5",
                                 "Proteobacteria" = "#ca250f")

plot_sigtab_RP_rotW1W22020 <- ggplot(sigtab_RP_rotW1W22020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1W22020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RP_rotW1W22020

ggsave("plot_sigtab_RP_rotW1W22020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 9, units = "cm", dpi = 300, device = "png")


########## W1WM 2020

#BS_rotW1WM2020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW1WM2020$log2FoldChange, sigtab_BS_rotW1WM2020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1WM2020$Phylum = factor(as.character(sigtab_BS_rotW1WM2020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW1WM2020$log2FoldChange, sigtab_BS_rotW1WM2020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1WM2020$Annotation = factor(as.character(sigtab_BS_rotW1WM2020$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW1WM2020 <- color.phylum_RH_rotW1W22020 <- c("Acidobacteriota" = "#ffc331",
                                                                "Actinobacteriota" = "#ff8b60",
                                                                "Bacteroidota"  = "#Bff4be",
                                                                "Campilobacterota" = "#800080",
                                                                "Chloroflexi" = "#2b8f22",
                                                                "Desulfobacterota" = "#964B00",
                                                                "Firmicutes" = "#4B69F5",
                                                                "Proteobacteria" = "#ca250f",
                                                                "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_BS_rotW1WM2020 <- ggplot(sigtab_BS_rotW1WM2020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW1WM2020) +
  scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW1WM2020

ggsave("plot_sigtab_BS_rotW1WM2020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 25, units = "cm", dpi = 300, device = "png")

#RH_rotW1WM2020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RH_rotW1WM2020$log2FoldChange, sigtab_RH_rotW1WM2020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW1WM2020$Phylum = factor(as.character(sigtab_RH_rotW1WM2020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RH_rotW1WM2020$log2FoldChange, sigtab_RH_rotW1WM2020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW1WM2020$Annotation = factor(as.character(sigtab_RH_rotW1WM2020$Annotation), levels=names(x))

# Define color palette
color.phylum_RH_rotW1WM2020 <- c("Acidobacteriota" = "#ffc331",
                                 "Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Chloroflexi" = "#2b8f22",
                                 "Cyanobacteria" = "#185113",                 
                                 "Deferribacterota" = "#d4ff31", 
                                 "Desulfobacterota" = "#964B00",
                                 "Firmicutes" = "#4B69F5",
                                 "Gemmatimonadota" = "#363781",
                                 "Myxococcota" = "#FFC0CB",
                                 "Patescibacteria" = "#4b0096", 
                                 "Proteobacteria" = "#ca250f",
                                 "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_RH_rotW1WM2020 <- ggplot(sigtab_RH_rotW1WM2020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RH_rotW1WM2020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RH_rotW1WM2020

ggsave("plot_sigtab_RH_rotW1WM2020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 20, units = "cm", dpi = 300, device = "png")


#RP_rotW1WM2020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RP_rotW1WM2020$log2FoldChange, sigtab_RP_rotW1WM2020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1WM2020$Phylum = factor(as.character(sigtab_RP_rotW1WM2020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RP_rotW1WM2020$log2FoldChange, sigtab_RP_rotW1WM2020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1WM2020$Annotation = factor(as.character(sigtab_RP_rotW1WM2020$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1WM2020 <- c("Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Firmicutes" = "#4B69F5",
                                 "Proteobacteria" = "#ca250f")

plot_sigtab_RP_rotW1WM2020 <- ggplot(sigtab_RP_rotW1WM2020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1WM2020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RP_rotW1WM2020

ggsave("plot_sigtab_RP_rotW1WM2020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 12, units = "cm", dpi = 300, device = "png")





########## W1W3 2020

#BS_rotW1W32020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW1W32020$log2FoldChange, sigtab_BS_rotW1W32020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W32020$Phylum = factor(as.character(sigtab_BS_rotW1W32020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW1W32020$log2FoldChange, sigtab_BS_rotW1W32020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W32020$Annotation = factor(as.character(sigtab_BS_rotW1W32020$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW1W32020 <- c("Acidobacteriota" = "#ffc331",
                                 "Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Chloroflexi" = "#2b8f22",
                                 "Desulfobacterota" = "#964B00",
                                 "Firmicutes" = "#4B69F5",
                                 "Proteobacteria" = "#ca250f",
                                 "Verrucomicrobiota" = "#ced2ce")

plot_sigtab_BS_rotW1W32020 <- ggplot(sigtab_BS_rotW1W32020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW1W32020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW1W32020

ggsave("plot_sigtab_BS_rotW1W32020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 9, units = "cm", dpi = 300, device = "png")

#RH_rotW1W32020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RH_rotW1W32020$log2FoldChange, sigtab_RH_rotW1W32020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW1W32020$Phylum = factor(as.character(sigtab_RH_rotW1W32020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RH_rotW1W32020$log2FoldChange, sigtab_RH_rotW1W32020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW1W32020$Annotation = factor(as.character(sigtab_RH_rotW1W32020$Annotation), levels=names(x))

# Define color palette
color.phylum_RH_rotW1W32020 <- c("Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Desulfobacterota" = "#964B00",
                                 "Firmicutes" = "#4B69F5",
                                 "Proteobacteria" = "#ca250f")

plot_sigtab_RH_rotW1W32020 <- ggplot(sigtab_RH_rotW1W32020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RH_rotW1W32020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RH_rotW1W32020

ggsave("plot_sigtab_RH_rotW1W32020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 8, units = "cm", dpi = 300, device = "png")


#RP_rotW1W32020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RP_rotW1W32020$log2FoldChange, sigtab_RP_rotW1W32020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1W32020$Phylum = factor(as.character(sigtab_RP_rotW1W32020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RP_rotW1W32020$log2FoldChange, sigtab_RP_rotW1W32020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1W32020$Annotation = factor(as.character(sigtab_RP_rotW1W32020$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1W32020 <- c("Actinobacteriota" = "#ff8b60", 
                                "Bacteroidota"  = "#Bff4be", 
                                "Patescibacteria" = "#4b0096", 
                                "Proteobacteria" = "#ca250f") 
                                

plot_sigtab_RP_rotW1W32020 <- ggplot(sigtab_RP_rotW1W32020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1W32020) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RP_rotW1W32020

ggsave("plot_sigtab_RP_rotW1W32020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 7, units = "cm", dpi = 300, device = "png")


########## W1W3 2021

#BS_rotW1W32021
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW1W32021$log2FoldChange, sigtab_BS_rotW1W32021$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W32021$Phylum = factor(as.character(sigtab_BS_rotW1W32021$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW1W32021$log2FoldChange, sigtab_BS_rotW1W32021$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W32021$Annotation = factor(as.character(sigtab_BS_rotW1W32021$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW1W32021 <- c("Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Firmicutes" = "#4B69F5",
                                 "Proteobacteria" = "#ca250f")

plot_sigtab_BS_rotW1W32021 <- ggplot(sigtab_BS_rotW1W32021, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW1W32021) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW1W32021

ggsave("plot_sigtab_BS_rotW1W32021.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 6, units = "cm", dpi = 300, device = "png")

#RH_rotW1W32021
## no taxa

#RP_rotW1W32021
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RP_rotW1W32021$log2FoldChange, sigtab_RP_rotW1W32021$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1W32021$Phylum = factor(as.character(sigtab_RP_rotW1W32021$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RP_rotW1W32021$log2FoldChange, sigtab_RP_rotW1W32021$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RP_rotW1W32021$Annotation = factor(as.character(sigtab_RP_rotW1W32021$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1W32021 <- c("Acidobacteriota" = "#ffc331",
                                 "Chloroflexi" = "#2b8f22",
                                 "Desulfobacterota" = "#964B00")

plot_sigtab_RP_rotW1W32021 <- ggplot(sigtab_RP_rotW1W32021, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1W32021) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RP_rotW1W32021

ggsave("plot_sigtab_RP_rotW1W32021.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 20, height = 4, units = "cm", dpi = 300, device = "png")







## The end! : )







########## W1W2 2021

#BS_rotW1W22021
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW1W22021$log2FoldChange, sigtab_BS_rotW1W22021$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W22021$Phylum = factor(as.character(sigtab_BS_rotW1W22021$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW1W22021$log2FoldChange, sigtab_BS_rotW1W22021$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1W22021$Annotation = factor(as.character(sigtab_BS_rotW1W22021$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW1W22021 <- c( "Firmicutes" = "#4B69F5",
                                  "Proteobacteria" = "#ca250f")

plot_sigtab_BS_rotW1W22021 <- ggplot(sigtab_BS_rotW1W22021, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW1W22021) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW1W22021

ggsave("plot_sigtab_BS_rotW1W22021.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 4, units = "cm", dpi = 300, device = "png")

#RH_rotW1W22021
## no taxa

#RP_rotW1W22021
## no taxa

########## W1WM 2021

#BS_rotW1WM2021
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW1WM2021$log2FoldChange, sigtab_BS_rotW1WM2021$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1WM2021$Phylum = factor(as.character(sigtab_BS_rotW1WM2021$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW1WM2021$log2FoldChange, sigtab_BS_rotW1WM2021$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW1WM2021$Annotation = factor(as.character(sigtab_BS_rotW1WM2021$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW1WM2021 <- c("Proteobacteria" = "#ca250f")

plot_sigtab_BS_rotW1WM2021 <- ggplot(sigtab_BS_rotW1WM2021, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW1WM2021) +
  scale_y_discrete("Annotation", labels = c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' = "ANPR group")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW1WM2021

ggsave("plot_sigtab_BS_rotW1WM2021.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 22, height = 4, units = "cm", dpi = 300, device = "png")

#RH_rotW1WM2021
## no taxa

#RP_rotW1WM2021
## no taxa

########## W2WM 2020

#BS_rotW2WM2021
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_BS_rotW2WM2021$log2FoldChange, sigtab_BS_rotW2WM2021$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW2WM2021$Phylum = factor(as.character(sigtab_BS_rotW2WM2021$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_BS_rotW2WM2021$log2FoldChange, sigtab_BS_rotW2WM2021$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_BS_rotW2WM2021$Annotation = factor(as.character(sigtab_BS_rotW2WM2021$Annotation), levels=names(x))

# Define color palette
color.phylum_BS_rotW2WM2021 <- c("Firmicutes" = "#A6CEE3")

plot_sigtab_BS_rotW2WM2021 <- ggplot(sigtab_BS_rotW2WM2021, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_BS_rotW2WM2021) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_BS_rotW2WM2021

ggsave("plot_sigtab_BS_rotW2WM2021.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 20, height = 4, units = "cm", dpi = 300, device = "png")

#RH_rotW2WM2021
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_RH_rotW2WM2021$log2FoldChange, sigtab_RH_rotW2WM2021$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW2WM2021$Phylum = factor(as.character(sigtab_RH_rotW2WM2021$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_RH_rotW2WM2021$log2FoldChange, sigtab_RH_rotW2WM2021$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_RH_rotW2WM2021$Annotation = factor(as.character(sigtab_RH_rotW2WM2021$Annotation), levels=names(x))

# Define color palette
color.phylum_RH_rotW2WM2021 <- c("Actinobacteriota" = "#1F78B4")

plot_sigtab_RH_rotW2WM2021 <- ggplot(sigtab_RH_rotW2WM2021, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum_RH_rotW2WM2021) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_RH_rotW2WM2021

ggsave("plot_sigtab_RH_rotW2WM2021.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 20, height = 4, units = "cm", dpi = 300, device = "png")


#RP_rotW2WM2021
## no taxa
