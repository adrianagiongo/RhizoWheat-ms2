##Define taxa deferentially abundant using DESeq2 on phyloseq pipeline

#Loading package
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")


# colors_RP <- list(
    Phylum =  c("Acidobacteriota" = "#ffc331",
                "Actinobacteriota" = "#ff8b60",
                "Bacteroidota"  = "#Bff4be",
                "Campilobacterota" = "#800080",
                "Chloroflexi" = "#2b8f22",
                "Cyanobacteria" = "#185113",
                "Deferribacterota" = "#d4ff31",
                "Desulfobacterota" = "#964B00",
                "Firmicutes" = "#4B69F5",
                "Gemmatimonadota" = "#363781",
                "Latescibacterota" = "#4a2500",
                "Myxococcota" = "#FFC0CB",
                "NB1-j" = "#b5814c",
                "Nitrospirota" = "#1361cf",
                "Patescibacteria" = "#4b0096",
                "Planctomycetota" = "#9fc5e8",
                "Proteobacteria" = "#ca250f",
                "Verrucomicrobiota" = "#ced2ce",
                "Zixibacteria" = "#e89fc5")




########## W1W2 2020

#RP_rotW1W22020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_RP_rotW1W22020$log2FoldChange, sigtab_rarefied_RP_rotW1W22020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_RP_rotW1W22020$Phylum = factor(as.character(sigtab_rarefied_RP_rotW1W22020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_RP_rotW1W22020$log2FoldChange, sigtab_rarefied_RP_rotW1W22020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_RP_rotW1W22020$Annotation = factor(as.character(sigtab_rarefied_RP_rotW1W22020$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1W22020 <- c("Actinobacteriota" = "#ff8b60",
                                 "Proteobacteria" = "#363781")

plot_sigtab_rarefied_RP_rotW1W22020 <- ggplot(sigtab_rarefied_RP_rotW1W22020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 20, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(legend.title = element_text(size = 24, color = "black")) +
  theme(legend.text = element_text(size = 24, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1W22020) +
  scale_x_continuous(limits = c(-7, 7), breaks = seq(-7, 7, by = 3.5), label = c("-7", "-3.5", "-0", "3.5", "7")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_RP_rotW1W22020

ggsave("plot_sigtab_rarefied_RP_rotW1W22020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 20, height = 10, units = "cm", dpi = 300, device = "png")


########## W1WM 2020

#RP_rotW1WM2020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_RP_rotW1WM2020$log2FoldChange, sigtab_rarefied_RP_rotW1WM2020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_RP_rotW1WM2020$Phylum = factor(as.character(sigtab_rarefied_RP_rotW1WM2020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_RP_rotW1WM2020$log2FoldChange, sigtab_rarefied_RP_rotW1WM2020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_RP_rotW1WM2020$Annotation = factor(as.character(sigtab_rarefied_RP_rotW1WM2020$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1WM2020 <- c("Actinobacteriota" = "#ff8b60",
                                 "Bacteroidota"  = "#Bff4be",
                                 "Proteobacteria" = "#363781")

plot_sigtab_rarefied_RP_rotW1WM2020 <- ggplot(sigtab_rarefied_RP_rotW1WM2020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 20, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 24, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(legend.title = element_text(size = 24, color = "black")) +
  theme(legend.text = element_text(size = 24, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1WM2020) +
  scale_x_continuous(limits = c(-7, 7), breaks = seq(-7, 7, by = 3.5), label = c("-7", "-3.5", "-0", "3.5", "7")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_RP_rotW1WM2020

ggsave("plot_sigtab_rarefied_RP_rotW1WM2020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 18, height = 9.5, units = "cm", dpi = 300, device = "png")


########## W1W3 2020

#RP_rotW1W32020
##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_RP_rotW1W32020$log2FoldChange, sigtab_rarefied_RP_rotW1W32020$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_RP_rotW1W32020$Phylum = factor(as.character(sigtab_rarefied_RP_rotW1W32020$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_RP_rotW1W32020$log2FoldChange, sigtab_rarefied_RP_rotW1W32020$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_RP_rotW1W32020$Annotation = factor(as.character(sigtab_rarefied_RP_rotW1W32020$Annotation), levels=names(x))

# Define color palette
color.phylum_RP_rotW1W32020 <- c("Actinobacteriota" = "#ff8b60",
                                 "Patescibacteria" = "#e89fc5", 
                                 "Proteobacteria" = "#363781") 
                                

plot_sigtab_rarefied_RP_rotW1W32020 <- ggplot(sigtab_rarefied_RP_rotW1W32020, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=4) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 20, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 24, color = "black")) +
  theme(axis.text.y = element_text(size = 24, color = "black")) +
  theme(legend.title = element_text(size = 24, color = "black")) +
  theme(legend.text = element_text(size = 24, color = "black")) +
  scale_color_manual(values = color.phylum_RP_rotW1W32020) +
  scale_x_continuous(limits = c(-7, 7), breaks = seq(-7, 7, by = 3.5), label = c("-7", "-3.5", "-0", "3.5", "7")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_RP_rotW1W32020

ggsave("plot_sigtab_rarefied_RP_rotW1W32020.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 19.8, height = 14, units = "cm", dpi = 300, device = "png")



###########################################
#arrange multiple ggplots in one single figure ALL
deseq2_plots_all <- ggarrange(plot_sigtab_rarefied_RP_rotW1W22020,
                                      plot_sigtab_rarefied_RP_rotW1W32020,
                                      plot_sigtab_rarefied_RP_rotW1WM2020,
                                      ncol =2, nrow =2, legend="right", widths = c(1, 1), heights = c(1, 1), common.legend = TRUE, label.y = 1)
deseq2_plots_all

ggsave("deseq2_plots_all.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/DESeq2_plot_jki_seq2_seq4/", width = 40, height = 20, units = "cm",dpi = 300, device = "svg")
###########################################


## The end! : )

