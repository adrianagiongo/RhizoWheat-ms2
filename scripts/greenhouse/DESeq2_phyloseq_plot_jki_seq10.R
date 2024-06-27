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

"Actinobacteriota" = "#ff8b60"
"Proteobacteria" = "#363781"
"Acidobacteriota" = "#ffc331"
"Firmicutes" = "#979191"
"Chloroflexi" = "#2b8f22"
"Nitrospirota" = "#1361cf"
"Crenarchaeota" = "#90c541"
"Gemmatimonadota" = "#ca250f"
"Verrucomicrobiota" = "#ced2ce"
"Bacteroidota"  = "#Bff4be"
"Myxococcota" = "#FFC0CB"
"Others" = "#3a3845"
##### RP 
########## W1W2 2020

##### This script uses a file from the output that was modified to contain only what presented more than 0.1% (relative abundance) 
sigtab_Go_W1_03 <- read.csv ("~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/sigtab_Go_W1_03.csv", row.names = 1)

##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_Go_W1_03$log2FoldChange, sigtab_Go_W1_03$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_Go_W1_03$Phylum = factor(as.character(sigtab_Go_W1_03$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_Go_W1_03$log2FoldChange, sigtab_Go_W1_03$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_Go_W1_03$Annotation = factor(as.character(sigtab_Go_W1_03$Annotation), levels=names(x))

# Define color palette
color.phylum <- c("Actinobacteriota" = "#ff8b60",
                  "Proteobacteria" = "#363781",
                  "Acidobacteriota" = "#ffc331",
                  "Firmicutes" = "#979191",
                  "Chloroflexi" = "#2b8f22",
                  "Nitrospirota" = "#1361cf",
                  "Crenarchaeota" = "#90c541",
                  "Gemmatimonadota" = "#ca250f",
                  "Verrucomicrobiota" = "#ced2ce",
                  "Bacteroidota"  = "#Bff4be",
                  "Desulfobacterota" = "#FFC0CB",
                  "Others" = "#3a3845")

plot_sigtab_Go_W1_03 <- ggplot(sigtab_Go_W1_03, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=5, alpha = 0.75) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum, breaks = c('Actinobacteriota',
                                                       'Proteobacteria')) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 3), label = c("-3", "0", "3")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_Go_W1_03

ggsave("plot_sigtab_Go_W1_03.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/DESeq2_plot_jki_seq10/", width = 16, height = 12, units = "cm", dpi = 300, device = "png")


