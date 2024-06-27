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

##### This script uses a file from the output that was modified to contain only what presented more than 0.1% (relative abundance) 
sigtab_RP_rotW1W22020 <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_01_RP_rotW1W22020.csv", row.names = 1)
sigtab_RP_rotW1WM2020 <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_01_RP_rotW1WM2020.csv", row.names = 1)
sigtab_RP_rotW1W32020 <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/sigtab_01_RP_rotW1W32020.csv", row.names = 1)


########## W1W2 2020

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
                                 "Proteobacteria" = "#363781",
                                 "Bacteroidota"  = "#Bff4be")

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




## The end! : )



