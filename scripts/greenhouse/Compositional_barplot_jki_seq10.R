##Creating barplot for selected dataset

#Loading package
library("microbiome")
library("phyloseq")
library("dplyr")
library("ggpubr")
library("vegan")

#count number of samples on dataset using "nsamples" function, and check the number of OTUs on the dataset
nsamples(psO_jki_seq10)
psO_jki_seq10

#plot histogram using "hist" function
hist(log10(taxa_sums(psO_jki_seq10)))

#create function aggregate_rare_taxa
aggregate_rare_taxa <- function (x, level, detection, prevalence, include.lowest=FALSE, ...) {
  x <- aggregate_taxa(x, level)
  rare <- rare_members(x, detection, prevalence, include.lowest)
  tax <- tax_table(x)
  inds <- which(rownames(tax) %in% rare)
  tax[inds, level] <- "Others"
  tax_table(x) <- tax
  tt <- tax_table(x)[, level]
  tax_table(x) <- tax_table(tt)
  aggregate_taxa(x, level)
}

##Aggregate rare taxa
#define total sums
total_samples_psO_jki_seq10 <- phyloseq::nsamples(psO_jki_seq10)
total_sum_psO_psO_jki_seq10 = sample_sums(psO_jki_seq10)

#aggregate based on phylum (detection higher than 40 copies and present in more than 85% of samples)
psO_jki_seq10_aggreg_rare_phy <- aggregate_rare_taxa(psO_jki_seq10, level="Phylum", detection = 600, prevalence =15/total_samples_psO_jki_seq10)
psO_jki_seq10_aggreg_rare_phy

##transform absolute to relative abundance using "transform" function
#for phylum
psO_jki_seq10_aggreg_rare_phy_rel <- microbiome::transform(psO_jki_seq10_aggreg_rare_phy, "compositional")
psO_jki_seq10_aggreg_rare_phy_rel

phy_colors <- c("Actinobacteriota" = "#ff8b60", 
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

#######
plot_bar_average_psO_jki_seq10_aggreg_rare_phy_rel <- plot_composition(psO_jki_seq10_aggreg_rare_phy_rel, otu.sort = "abundance", sample.sort = "Firmicutes", x.label = "Site_rot_tre", plot.type = "barplot", average_by = "Site_rot_tre", verbose = FALSE) +
  theme_bw() +
  scale_fill_manual("Phylum", values =  phy_colors) +
  #scale_fill_brewer("Phylum", palette = "Paired") +
  labs(x=" ", y= "Relative abundance") +
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust =0.5, hjust = 0.5, colour = "black"), axis.text.y = element_text(size = 20, colour = "black")) +
  theme(axis.title.x = element_text(size = 20, colour = "black"), axis.title.y = element_text(size = 20, colour = "black")) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size = 20, colour = "black", face = "italic")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), label = c("0", "25", "50", "75", "100"))

plot_bar_average_psO_jki_seq10_aggreg_rare_phy_rel

ggsave("plot_bar_average_psO_jki_seq10_aggreg_rare_phy_rel.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Compositional_barplot_jki_seq10/", width = 16, height = 15, units = "cm",dpi = 300, device = "png")



# The end! : )
