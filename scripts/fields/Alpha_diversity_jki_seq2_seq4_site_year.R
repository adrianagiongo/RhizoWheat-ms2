###########
##### Correspondent to Figure 1A and Suppl Figure S1A
###########
#Go = "#1d5b65"   Go 2020 = "#216974"  Go 2021 = "#6baea5"
#Ki = "#A34828"   Ki 2020 = "#D1711F"  Ki 2021 = "#ebaf7a"

#Create rarefied data from phyloseq object using vegan package or microbiome package
#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("ggplot2")

##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq2_seq4_rarefied_1 <- richness(psO_jki_seq2_seq4_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq2_seq4_rarefied_2 <- evenness(psO_jki_seq2_seq4_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq2_seq4_rarefied_3 <- microbiome::diversity(psO_jki_seq2_seq4_rarefied, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity
#Prepare file using meta function
psO_jki_seq2_seq4_rarefied.meta <- meta(psO_jki_seq2_seq4_rarefied)

#select column from output file to be plotted and tested (only one by time)
#for Go_2020
psO_jki_seq2_seq4_rarefied.meta$observed <- alfa_div_psO_jki_seq2_seq4_rarefied_1$observed
psO_jki_seq2_seq4_rarefied.meta$chao1    <- alfa_div_psO_jki_seq2_seq4_rarefied_1$chao1
psO_jki_seq2_seq4_rarefied.meta$pielou   <- alfa_div_psO_jki_seq2_seq4_rarefied_2$pielou
psO_jki_seq2_seq4_rarefied.meta$shannon  <- alfa_div_psO_jki_seq2_seq4_rarefied_3$shannon

head(psO_jki_seq2_seq4_rarefied.meta)
write.csv(psO_jki_seq2_seq4_rarefied.meta, "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Tables_jki_seq2_seq4/psO_jki_seq2_seq4_rarefied_meta.csv")


#Select variable of comparison
#alfa_div_comparison_Site <- list(c("Go","Ki"))


#Plots SHANNON
#Select variable of comparison
#alfa_div_comparison_Site <- list(c("Go","Ki"))

########## Figure 1A
#to compare 
plot_psO_jki_seq2_seq4_rarefied_shannon_Site <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site", "shannon", fill = "Site", color = "Site", palette = c("#1d5b65", "#A34828"), ylab = "Shannon index") +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  scale_y_continuous(limits = c(5.5, 7.5), breaks = seq(5.5, 7.5, by = 1.0), label = c("5.5", "6.5", "7.5")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5 ), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 5.5, label.x = 1, size=7)
 plot_psO_jki_seq2_seq4_rarefied_shannon_Site

ggsave("plot_psO_jki_seq2_seq4_rarefied_shannon_Site.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 10, height = 12, units = "cm", dpi = 300, device = "svg")


#Plots CHAO1
plot_psO_jki_seq2_seq4_rarefied_chao1_Site <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site", "chao1", fill = "Site", color = "Site", palette = c("#1d5b65", "#A34828"), ylab = "Chao1 index") +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  scale_y_continuous(limits = c(1500, 4000), breaks = seq(1500, 4000, by = 1000), label = c("1500", "2500", "3500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 1500, label.x = 1, size=7)
plot_psO_jki_seq2_seq4_rarefied_chao1_Site

ggsave("plot_psO_jki_seq2_seq4_rarefied_chao1_Site.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 10, height = 12, units = "cm", dpi = 300, device = "svg")

#Plots PIELOU
plot_psO_jki_seq2_seq4_rarefied_pielou_Site <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site", "pielou", fill = "Site", color = "Site", palette = c("#1d5b65", "#A34828"), ylab = "Pielou index") +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  scale_y_continuous(limits = c(0.85, 0.95), breaks = seq(0.85, 0.95, by = 0.05), label = c("0.85", "0.9", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 0,vjust =0.5, hjust = 0.5), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  stat_compare_means(method = "wilcox.test", label= "p", label.y = 0.85, label.x = 1, size=7)
plot_psO_jki_seq2_seq4_rarefied_pielou_Site

ggsave("plot_psO_jki_seq2_seq4_rarefied_pielou_Site.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 10, height = 12, units = "cm", dpi = 300, device = "svg")


###########################################
#arrange multiple ggplots in one single figure ALL
alpha_plots_all <- ggarrange(plot_psO_jki_seq2_seq4_rarefied_shannon_Site, plot_psO_jki_seq2_seq4_rarefied_chao1_Site, 
                             plot_psO_jki_seq2_seq4_rarefied_pielou_Site, 
                             ncol =2, nrow =2, legend="right", widths = c(1, 1), common.legend = TRUE, label.y = 1)
alpha_plots_all

ggsave("alpha_plots_all.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 25, height = 20, units = "cm",dpi = 300, device = "svg")
###########################################





### I stopped here (04.03.2024)!




#### Separated by year
alfa_div_comparison_all <- list(c("Go_2020", "Go_2021"),c("Ki_2020", "Ki_2021"))

#Plots SHANNON
#to compare 
plot_psO_jki_seq2_seq4_rarefied_shannon_Site_year <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site_year", "shannon", fill = "Site_year", color = "Site_year", palette = c("#216974", "#6baea5", "#D1711F", "#ebaf7a"), ylab = "Shannon index", order = c("Go_2020", "Go_2021", "Ki_2020", "Ki_2021")) +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  scale_y_continuous(limits = c(6, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size = 22), strip.text.y = element_text(size = 22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) 
  #stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_all, label = "p.format",  bracket.size = .3, size=4, label.y = c(7.7,7.7)) 
plot_psO_jki_seq2_seq4_rarefied_shannon_Site_year

ggsave("plot_psO_jki_seq2_seq4_rarefied_shannon_Site_year.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 12, height = 9, units = "cm", dpi = 300, device = "svg")


#Plots Chao1
#to compare 
plot_psO_jki_seq2_seq4_rarefied_chao1_Site_year <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site_year", "chao1", fill = "Site_year", color = "Site_year", palette = c("#216974", "#6baea5", "#D1711F", "#ebaf7a"), ylab = "Chao1 index", order = c("Go_2020", "Go_2021", "Ki_2020", "Ki_2021")) +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  scale_y_continuous(limits = c(1500, 4500), breaks = seq(1500, 4500, by = 1000), label = c("1500", "2500", "3500", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size = 22), strip.text.y = element_text(size = 22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) 
  #stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_all, label = "p.format",  bracket.size = .3, size=4, label.y = c(4200,4200)) 
plot_psO_jki_seq2_seq4_rarefied_chao1_Site_year

ggsave("plot_psO_jki_seq2_seq4_rarefied_chao1_Site_year.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 12, height = 9, units = "cm", dpi = 300, device = "svg")

#Plots Pielou
#to compare 
plot_psO_jki_seq2_seq4_rarefied_pielou_Site_year <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site_year", "pielou", fill = "Site_year", color = "Site_year", palette = c("#216974", "#6baea5", "#D1711F", "#ebaf7a"), ylab = "Pielou index", order = c("Go_2020", "Go_2021", "Ki_2020", "Ki_2021")) +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  scale_y_continuous(limits = c(0.85, 0.95), breaks = seq(0.85, 0.95, by = 0.05), label = c("0.85", "0.9", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size = 22), strip.text.y = element_text(size = 22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) 
  #stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_all, label = "p.format",  bracket.size = .3, size=4, label.y = c(0.93,0.93)) 
plot_psO_jki_seq2_seq4_rarefied_pielou_Site_year

ggsave("plot_psO_jki_seq2_seq4_rarefied_pielou_Site_year.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 12, height = 9, units = "cm", dpi = 300, device = "svg")






#### Separated by year BUT ALSO by STAGE (in facet)
alfa_div_comparison_all <- list(c("Go_2020", "Go_2021"),c("Ki_2020", "Ki_2021"))

#Plots SHANNON
#to compare 
plot_psO_jki_seq2_seq4_rarefied_shannon_Site_year_facet <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site_year", "shannon", fill = "Site_year", color = "Site_year", palette = c("#216974", "#6baea5", "#D1711F", "#ebaf7a"), ylab = "Shannon index", order = c("Go_2020", "Go_2021", "Ki_2020", "Ki_2021")) +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  facet_grid(~Stage) +
  scale_y_continuous(limits = c(6, 8), breaks = seq(6, 8, by = 1.0), label = c("6.0", "7.0", "8.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size = 22), strip.text.y = element_text(size = 22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_all, label = "p.format",  bracket.size = .3, size=4, label.y = c(7.7,7.7)) 
plot_psO_jki_seq2_seq4_rarefied_shannon_Site_year_facet

ggsave("plot_psO_jki_seq2_seq4_rarefied_shannon_Site_year_facet.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 12, height = 12, units = "cm", dpi = 300, device = "svg")


#Plots Chao1
#to compare 
plot_psO_jki_seq2_seq4_rarefied_chao1_Site_year_facet <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site_year", "chao1", fill = "Site_year", color = "Site_year", palette = c("#216974", "#6baea5", "#D1711F", "#ebaf7a"), ylab = "Chao1 index", order = c("Go_2020", "Go_2021", "Ki_2020", "Ki_2021")) +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  facet_grid(~Stage) +
  scale_y_continuous(limits = c(1500, 4500), breaks = seq(1500, 4500, by = 1000), label = c("1500", "2500", "3500", "4500")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size = 22), strip.text.y = element_text(size = 22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_all, label = "p.format",  bracket.size = .3, size=4, label.y = c(4200,4200)) 
plot_psO_jki_seq2_seq4_rarefied_chao1_Site_year_facet

ggsave("plot_psO_jki_seq2_seq4_rarefied_chao1_Site_year_facet.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 13, height = 12, units = "cm", dpi = 300, device = "svg")

#Plots Pielou
#to compare 
plot_psO_jki_seq2_seq4_rarefied_pielou_Site_year_facet <- ggboxplot(psO_jki_seq2_seq4_rarefied.meta, "Site_year", "pielou", fill = "Site_year", color = "Site_year", palette = c("#216974", "#6baea5", "#D1711F", "#ebaf7a"), ylab = "Pielou index", order = c("Go_2020", "Go_2021", "Ki_2020", "Ki_2021")) +
  theme_bw() +  geom_boxplot(alpha = 0.1, width = 0.7) +
  facet_grid(~Stage) +
  scale_y_continuous(limits = c(0.85, 0.95), breaks = seq(0.85, 0.95, by = 0.05), label = c("0.85", "0.9", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 22, colour = "black")) +
  theme(axis.text.x = element_text(size=22, angle = 90,vjust =0.5, hjust = 1 ), axis.text.y = element_text(size=22)) +
  theme(strip.text.x = element_text(size = 22), strip.text.y = element_text(size = 22)) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison_all, label = "p.format",  bracket.size = .3, size=4, label.y = c(0.93,0.93)) 
plot_psO_jki_seq2_seq4_rarefied_pielou_Site_year_facet

gsave("plot_psO_jki_seq2_seq4_rarefied_pielou_Site_year_facet.svg", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Alpha_div_jki_seq2_seq4/", width = 12, height = 12, units = "cm", dpi = 300, device = "svg")

## Have a good day!  : )