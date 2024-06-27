##Define taxa deferentially abundant using DESeq2 on phyloseq pipeline

#Loading package
library("phyloseq")
library("ggplot2")
library("DESeq2")
library("RColorBrewer")
library("svglite")

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


##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_Go_W1$log2FoldChange, sigtab_rarefied_Go_W1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Go_W1$Phylum = factor(as.character(sigtab_rarefied_Go_W1$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_Go_W1$log2FoldChange, sigtab_rarefied_Go_W1$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Go_W1$Annotation = factor(as.character(sigtab_rarefied_Go_W1$Annotation), levels=names(x))


### Go W1
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
                  "Myxococcota" = "#FFC0CB",
                  "Others" = "#3a3845")

plot_sigtab_rarefied_Go_W1 <- ggplot(sigtab_rarefied_Go_W1, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=5, alpha = 0.80) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum, breaks = c('Actinobacteriota', 'Bacteroidota', 'Chloroflexi',
                                                        'Firmicutes', 'Gemmatimonadota', 'Myxococcota',
                                                       'Proteobacteria', 'Verrucomicrobiota')) +
  scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, by = 2), label = c("-2", "0", "2", "4")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_Go_W1

ggsave("plot_sigtab_rarefied_Go_W1.svg", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/DESeq2_plot_jki_seq10/", width = 20, height = 24, units = "cm", dpi = 300, device = "svg")



##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_Go_WM$log2FoldChange, sigtab_rarefied_Go_WM$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Go_WM$Phylum = factor(as.character(sigtab_rarefied_Go_WM$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_Go_WM$log2FoldChange, sigtab_rarefied_Go_WM$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Go_WM$Annotation = factor(as.character(sigtab_rarefied_Go_WM$Annotation), levels=names(x))

# Define color palette
color.phylum <- c("Actinobacteriota" = "#ff8b60",
                  "Proteobacteria" = "#363781",
                  "Acidobacteriota" = "#ffc331",
                  "Firmicutes" = "#979191",
                  "Chloroflexi" = "#2b8f22",
                  "Nitrospirota" = "#036ab1",
                  "Crenarchaeota" = "#90c541",
                  "Gemmatimonadota" = "#ca250f",
                  "Verrucomicrobiota" = "#ced2ce",
                  "Bacteroidota"  = "#Bff4be",
                  "Myxococcota" = "#FFC0CB",
                  "NB1_j" = "#1361cf")

plot_sigtab_rarefied_Go_WM <- ggplot(sigtab_rarefied_Go_WM, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=5, alpha = 0.80) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum, breaks = c('Actinobacteriota', 'Bacteroidota', 'Chloroflexi',
                                                        'Firmicutes', 'Gemmatimonadota', 'Myxococcota', 'NB1_j',
                                                        'Proteobacteria', 'Verrucomicrobiota')) +
  #scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 3), label = c("-3", "0", "3")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_Go_WM

ggsave("plot_sigtab_rarefied_Go_WM.svg", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/DESeq2_plot_jki_seq10/", width = 20, height = 24, units = "cm", dpi = 300, device = "svg")





##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_Ki_W1$log2FoldChange, sigtab_rarefied_Ki_W1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Ki_W1$Phylum = factor(as.character(sigtab_rarefied_Ki_W1$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_Ki_W1$log2FoldChange, sigtab_rarefied_Ki_W1$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Ki_W1$Annotation = factor(as.character(sigtab_rarefied_Ki_W1$Annotation), levels=names(x))


### Ki W1
# Define color palette
color.phylum <- c("Actinobacteriota" = "#ff8b60",
                  "Proteobacteria" = "#363781",
                  "Acidobacteriota" = "#ffc331",
                  "Firmicutes" = "#979191",
                  "Chloroflexi" = "#2b8f22",
                  "Desulfobacterota" = "#1361cf",
                  "Crenarchaeota" = "#90c541",
                  "Gemmatimonadota" = "#ca250f",
                  "Verrucomicrobiota" = "#ced2ce",
                  "Bacteroidota"  = "#Bff4be",
                  "Myxococcota" = "#FFC0CB",
                  "Others" = "#3a3845")

plot_sigtab_rarefied_Ki_W1 <- ggplot(sigtab_rarefied_Ki_W1, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=5, alpha = 0.80) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum, breaks = c('Acidobacteriota', 'Actinobacteriota', 'Bacteroidota', 
                                                        'Chloroflexi', 'Desulfobacterota', 'Firmicutes', 
                                                        'Myxococcota', 'Proteobacteria', 'Verrucomicrobiota')) +
  #scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 3), label = c("-3", "0", "3")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_Ki_W1

ggsave("plot_sigtab_rarefied_Ki_W1.svg", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/DESeq2_plot_jki_seq10/", width = 20, height = 24, units = "cm", dpi = 300, device = "svg")



##Convert variables characters in factors
# Phylum order
x = tapply(sigtab_rarefied_Ki_W3$log2FoldChange, sigtab_rarefied_Ki_W3$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Ki_W3$Phylum = factor(as.character(sigtab_rarefied_Ki_W3$Phylum), levels=names(x))

# Annotation order
x = tapply(sigtab_rarefied_Ki_W3$log2FoldChange, sigtab_rarefied_Ki_W3$Annotation, function(x) max(x))
x = sort(x, TRUE)
sigtab_rarefied_Ki_W3$Annotation = factor(as.character(sigtab_rarefied_Ki_W3$Annotation), levels=names(x))

# Define color palette
color.phylum <- c("Actinobacteriota" = "#ff8b60",
                  "Proteobacteria" = "#363781",
                  "Acidobacteriota" = "#ffc331",
                  "Firmicutes" = "#979191",
                  "Chloroflexi" = "#2b8f22",
                  "Nitrospirota" = "#036ab1",
                  "Crenarchaeota" = "#90c541",
                  "Gemmatimonadota" = "#ca250f",
                  "Verrucomicrobiota" = "#ced2ce",
                  "Bacteroidota"  = "#Bff4be",
                  "Myxococcota" = "#FFC0CB",
                  "NB1_j" = "#1361cf")

plot_sigtab_rarefied_Ki_W3 <- ggplot(sigtab_rarefied_Ki_W3, aes(y=Annotation, x=log2FoldChange, color=Phylum)) +
  geom_point(size=5, alpha = 0.80) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5, size = 16, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 16, color = "black")) +
  theme(axis.text.y = element_text(size = 16, color = "black")) +
  theme(legend.title = element_text(size = 16, color = "black")) +
  theme(legend.text = element_text(size = 16, color = "black")) +
  scale_color_manual(values = color.phylum, breaks = c('Actinobacteriota', 'Bacteroidota', 'Chloroflexi',
                                                       'Firmicutes', 'Proteobacteria', 'Verrucomicrobiota')) +
  #scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 3), label = c("-3", "0", "3")) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray", size=1)
plot_sigtab_rarefied_Ki_W3

ggsave("plot_sigtab_rarefied_Ki_W3.svg", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/DESeq2_plot_jki_seq10/", width = 20, height = 24, units = "cm", dpi = 300, device = "svg")
