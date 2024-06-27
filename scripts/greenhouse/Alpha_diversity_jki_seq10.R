#Create rarefied data from phyloseq object using vegan package or microbiome package
#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("agricolae")


##Calculates the alpha diversity of microbial communities
#Using dominance, richness, evenness, and diversity, all functions from microbiome package
alfa_div_psO_jki_seq10_rarefied_1 <- richness(psO_jki_seq10_rarefied, c("observed","chao1"), detection = 0)
alfa_div_psO_jki_seq10_rarefied_2 <- evenness(psO_jki_seq10_rarefied, index = 'pielou', zeroes = TRUE, detection = 0)
alfa_div_psO_jki_seq10_rarefied_3 <- microbiome::diversity(psO_jki_seq10_rarefied, index = 'shannon', zeroes = TRUE)

##################################################
##Boxplot alpha diversity 
#Prepare file using meta function
psO_jki_seq10_rarefied.meta <- meta(psO_jki_seq10_rarefied)

#select column from output file to be plotted and tested (only one by time)

psO_jki_seq10_rarefied.meta$observed <- alfa_div_psO_jki_seq10_rarefied_1$observed
psO_jki_seq10_rarefied.meta$chao1    <- alfa_div_psO_jki_seq10_rarefied_1$chao1
psO_jki_seq10_rarefied.meta$pielou   <- alfa_div_psO_jki_seq10_rarefied_2$pielou
psO_jki_seq10_rarefied.meta$shannon  <- alfa_div_psO_jki_seq10_rarefied_3$shannon

head(psO_jki_seq10_rarefied.meta)
tail(psO_jki_seq10_rarefied.meta)
write.csv(psO_jki_seq10_rarefied.meta, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/psO_jki_seq10_rarefied_meta.csv")

##check for the distribution of the diversity using "hist" function to plot, 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05) 
#then the null hypotesis (the population is normally distributed) is rejected
#for data with only two variable

hist(psO_jki_seq10_rarefied.meta$shannon)
shapiro.test(psO_jki_seq10_rarefied.meta$shannon)
qqnorm(psO_jki_seq10_rarefied.meta$shannon)

hist(psO_jki_seq10_rarefied.meta$chao1)
shapiro.test(psO_jki_seq10_rarefied.meta$chao1)
qqnorm(psO_jki_seq10_rarefied.meta$chao1)

hist(psO_jki_seq10_rarefied.meta$pielou)
shapiro.test(psO_jki_seq10_rarefied.meta$pielou)
qqnorm(psO_jki_seq10_rarefied.meta$pielou)


# Kruskal with Bonferoni using Agricolae  <-- isso gera as letras que precisam ser colocadas no topo de cada boxplot ####
Kruskal_shannon <-with(psO_jki_seq10_rarefied.meta,kruskal(shannon,Site_rot,group=TRUE, main="Site_rot", console=TRUE))
Kruskal_shannon
# Ki_W3_C     32.25      a
# BL_W1_C     31.25     ab
# Go_W1_C     29.00     ab
# Ki_W1_Fg1   20.50    abc
# Ki_W1_C     20.00    abc
# BL_W1_Fg1   18.75    abc
# Go_W1_Fg1   18.00     bc
# Ki_W3_Fg1   17.50     bc
# Go_WM_Fg1   10.75      c
# Go_WM_C      7.00      c

Kruskal_chao1 <-with(psO_jki_seq10_rarefied.meta,kruskal(chao1,Site_rot,group=TRUE, main="Site_rot", console=TRUE))
Kruskal_chao1
# BL_W1_C   37.50      a
# BL_W1_Fg1 34.00     ab
# Ki_W3_C   28.75    abc
# Go_W1_C   27.25     bc
# Ki_W3_Fg1 20.00     cd
# Ki_W1_C   14.75     de
# Ki_W1_Fg1 14.25     de
# Go_W1_Fg1 14.00     de
# Go_WM_C    7.25      e
# Go_WM_Fg1  7.25      e

Kruskal_pielou <-with(psO_jki_seq10_rarefied.meta,kruskal(pielou,Site_rot,group=TRUE, main="Site_rot", console=TRUE))
Kruskal_pielou
# Ki_W3_C    33.25      a
# Go_W1_C    29.50     ab
# BL_W1_C    25.50    abc
# Ki_W1_Fg1  23.00    abc
# Ki_W1_C    22.50    abc
# Go_W1_Fg1  19.25   abcd
# Ki_W3_Fg1  17.50    bcd
# BL_W1_Fg1  14.75     cd
# Go_WM_Fg1  12.25     cd
# Go_WM_C     7.50      d

#Select variable of comparison 
alfa_div_comparison <- list(c("Go_W1_C", "Go_W1_Ggt"), c("Go_WM_C", "Go_WM_Ggt"),
                            c("Ki_W1_C", "Ki_W1_Ggt"), c("Ki_W3_C", "Ki_W3_Ggt"))

## Shannon  (Statistics)
plot_psO_jki_seq10_rarefied_shannon <- ggboxplot(psO_jki_seq10_rarefied.meta, "Site_rot_tre", "shannon", fill = "Site_rot_tre",
                                                palette = c( "#6baea5","#6baea5", "#216974", "#216974", 
                                                             "#ebaf7a", "#ebaf7a","#D1711F","#D1711F"),
                                                ylab = "Shannon index", order = c("Go_W1_C", "Go_W1_Ggt", "Go_WM_C", "Go_WM_Ggt", 
                                                                                  "Ki_W1_C", "Ki_W1_Ggt", "Ki_W3_C", "Ki_W3_Ggt"),
                                                bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +
  #scale_y_continuous(limits = c(1, 5), breaks = seq(1, 5, by = 2), label = c("1.0", "3.0", "5.0")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  theme(axis.text.x = element_text(size=20, angle = 90, vjust =0.5, hjust = 1), axis.text.y = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size = 20, colour = "black")) 
  #annotate("text", x = 1:12, y = 5, label = c("abc", "de", "a", "ab", "bcd", "de", "e", "e", "cde", "cde", "abc", "de"), size = 5)
  # annotate("text", x = 1, y = 6.62, label = c("bc"), size = 5) +
  # annotate("text", x = 2, y = 7.15, label = c("abc"), size = 5) +
  # annotate("text", x = 3, y = 6.7, label = c("bc"), size = 5) +
  # annotate("text", x = 4, y = 7.1, label = c("ab"), size = 5) +
  # annotate("text", x = 5, y = 6.79, label = c("c"), size = 5) +
  # annotate("text", x = 6, y = 6.92, label = c("ab"), size = 5) +
  # annotate("text", x = 7, y = 6.88, label = c("abc"), size = 5) +
  # annotate("text", x = 8, y = 7.04, label = c("a"), size = 5) +
  #stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison, label= "p", bracket.size = .3, size=4, label.y = c(7.3,7.3,7.3,7.3)) 
plot_psO_jki_seq10_rarefied_shannon

ggsave("plot_psO_jki_seq10_rarefied_shannon.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Alpha_div_jki_seq10/", width = 20, height = 14, units = "cm",dpi = 200)

#Plots CHAO1
plot_psO_jki_seq10_rarefied_chao1 <- ggboxplot(psO_jki_seq10_rarefied.meta, "Site_rot_tre", "chao1", fill = "Site_rot_tre",
                                               palette = c( "#6baea5","#6baea5", "#216974", "#216974", 
                                                            "#ebaf7a", "#ebaf7a","#D1711F","#D1711F"),
                                                       ylab = "Chao1 index", order = c("Go_W1_C", "Go_W1_Ggt", "Go_WM_C", "Go_WM_Ggt",
                                                                                         "Ki_W1_C", "Ki_W1_Ggt", "Ki_W3_C", "Ki_W3_Ggt"),
                                                       bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +
  #scale_y_continuous(limits = c(2000, 3250), breaks = seq(2000, 3250, by = 250), label = c("2000", "2250", "2500", "2750", "3000", "3250")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  theme(axis.text.x = element_text(size=20, angle = 90, vjust =0.5, hjust = 1), axis.text.y = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size = 20, colour = "black")) +
  # #annotate("text", x = 1:6, y = 3180, label = c("a", "ab", "ab", "b", "ab", "ab"), size = 5)
  # annotate("text", x = 1, y = 2640, label = c("b"), size = 5) +
  # annotate("text", x = 2, y = 3000, label = c("b"), size = 5) +
  # annotate("text", x = 3, y = 2850, label = c("b"), size = 5) +
  # annotate("text", x = 4, y = 3020, label = c("ab"), size = 5) +
  # annotate("text", x = 5, y = 2820, label = c("b"), size = 5) +
  # annotate("text", x = 6, y = 2890, label = c("ab"), size = 5) +
  # annotate("text", x = 7, y = 2740, label = c("b"), size = 5) +
  # annotate("text", x = 8, y = 3100, label = c("a"), size = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison, label= "p", bracket.size = .3, size=4, label.y = c(3250, 3250, 3250, 3250))
plot_psO_jki_seq10_rarefied_chao1

ggsave("plot_psO_jki_seq10_rarefied_chao1.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Alpha_div_jki_seq10/", width = 20, height = 14, units = "cm",dpi = 200)


#Plots PIELOU
plot_psO_jki_seq10_rarefied_pielou <- ggboxplot(psO_jki_seq10_rarefied.meta, "Site_rot_tre", "pielou", fill = "Site_rot_tre",
                                                palette = c( "#6baea5","#6baea5", "#216974", "#216974", 
                                                             "#ebaf7a", "#ebaf7a","#D1711F","#D1711F"),
                                                        ylab = "Pielou index", order = c("BL_W1_C", "BL_W1_Ggt", 
                                                                                        "Go_W1_C", "Go_W1_Ggt", "Go_WM_C", "Go_WM_Ggt",
                                                                                        "Ki_W1_C", "Ki_W1_Ggt", "Ki_W3_C", "Ki_W3_Ggt"),
                                                        bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +
  #scale_y_continuous(limits = c(0.85, 0.95), breaks = seq(0.85, 0.95, by = 0.05), label = c("0.85", "0.90", "0.95")) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20, colour = "black")) +
  theme(axis.text.x = element_text(size=20, angle = 90, vjust =0.5, hjust = 1), axis.text.y = element_text(size=20)) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size = 20, colour = "black")) +
  # #annotate("text", x = 1:6, y = 0.90, label = c("ab", "bc", "abc", "c", "ab", "a"), size = 5) 
  # annotate("text", x = 1, y = 0.845, label = c("c"), size = 5) +
  # annotate("text", x = 2, y = 0.893, label = c("abc"), size = 5) +
  # annotate("text", x = 3, y = 0.846, label = c("bc"), size = 5) +
  # annotate("text", x = 4, y = 0.888, label = c("a"), size = 5) +
  # annotate("text", x = 5, y = 0.859, label = c("bc"), size = 5) +
  # annotate("text", x = 6, y = 0.875, label = c("ab"), size = 5) +
  # annotate("text", x = 7, y = 0.872, label = c("abc"), size = 5) +
  # annotate("text", x = 8, y = 0.882, label = c("a"), size = 5) +
  stat_compare_means(method = "wilcox.test", comparisons = alfa_div_comparison, label= "p", bracket.size = .3, size=4, label.y = c(0.91, 0.91, 0.91, 0.91))
plot_psO_jki_seq10_rarefied_pielou

ggsave("plot_psO_jki_seq10_rarefied_pielou.png", path = "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Alpha_div_jki_seq10/", width = 20, height = 14, units = "cm",dpi = 200)


## The end! : )

