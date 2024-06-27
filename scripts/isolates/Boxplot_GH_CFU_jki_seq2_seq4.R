#Load packages
library("readxl")
library("ggplot2")
library("ggpubr")
library("dplyr")

# GH13
########### CFU with non-autoclaved soil
GH13_CFU <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/CFU_GH13.csv")


##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypothesis (the population is normally distributed) is rejected
hist(GH13_CFU$CFU)
shapiro.test(GH13_CFU$CFU)
qqnorm(GH13_CFU$CFU)

compare_treats <- list(c("B_pum_before", "B_pum_after"),
                       c("P_bra_before", "P_bra_after"),
                       c("P_pro_before", "P_pro_after"),
                       c("S_pan_before", "S_pan_after"))


################ Plot with statistics
boxplot_GH13_CFU <- ggboxplot(GH13_CFU, "Treat", "CFU", fill = "#ebaf7a",
                                   palette = c("#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a"),
                                   ylab = "CFU",
                                   bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +
  #scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5), label = c("0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
  #stat_compare_means(method = "wilcox.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c(4,5,6,7)) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))+
  scale_x_discrete("Treat", labels = c('B_pum' = "B. pum",
                                           'P_bra' = "P. bra",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan"))

boxplot_GH13_CFU

ggsave("boxplot_GH13_CFU.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 12, height = 22, units = "cm",dpi = 200)

################ Plot with facet by Test
# Reorder levels of the Test variable
GH13_CFU$Test <- factor(GH13_CFU$Test,
                             levels = c("CFU_before", "CFU_after"))

boxplot_GH13_CFU_facet <- ggboxplot(GH13_CFU, "Treat", "CFU", fill = "Treat",
                                         palette = c("#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a"),
                                         # order = c("Silty_W1_Control", "Silty_W1_Ggt",
                                         #           "Silty_WM_Control", "Silty_WM_Ggt",
                                         #           "Sandy_W1_Control", "Sandy_W1_Ggt",
                                         #           "Sandy_W3_Control", "Sandy_W3_Ggt"),
                                         ylab = "CFU (log)",
                                         bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  
  #geom_boxplot(alpha = 0.1) +
  facet_grid(~Test, scales = "free", space="free_x") +
  # scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 2.5), label = c("0", "2.5", "5.0")) +
  #stat_compare_means(method = "wilcox.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c(5,5,5,5,5,5,5,5,5,5)) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  scale_x_discrete("Treat", labels = c('B_pum_before' = "B. pum",
                                       'P_bra_before' = "P. bra",
                                       'P_pro_before' = "P. pro",
                                       'S_pan_before' = "S. pan",
                                       'B_pum_after' = "B. pum",
                                       'P_bra_after' = "P. bra",
                                       'P_pro_after' = "P. pro",
                                       'S_pan_after' = "S. pan"))
boxplot_GH13_CFU_facet

ggsave("boxplot_GH13_CFU_facet.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 18, height = 10, units = "cm",dpi = 200)















# GH11_GH13
########### CFU with non-autoclaved soil
GH11_GH13_CFU <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/CFU_GH11_GH13.csv")


##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypothesis (the population is normally distributed) is rejected
hist(GH11_GH13_CFU$CFU)
shapiro.test(GH11_GH13_CFU$CFU)
qqnorm(GH11_GH13_CFU$CFU)

compare_treats <- list(c("B_pumilus_45_39_Ggt_before", "B_pumilus_45_39_Ggt_after"),
                       c("P_brassicacearum_37_15_Ggt_before", "P_brassicacearum_37_15_Ggt_after"),
                       c("P_proteolytica_68_48_Ggt_before", "P_proteolytica_68_48_Ggt_after"),
                       c("S_panacis_40_16_Ggt_before", "S_panacis_40_16_Ggt_after"))


################ Plot with statistics
boxplot_GH11_GH13_CFU <- ggboxplot(GH11_GH13_CFU, "Treat", "CFU", fill = "Treat",
                               palette = c("#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a"),
                               ylab = "CFU (log)",
                               bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +
  #scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5), label = c("0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
  #stat_compare_means(method = "wilcox.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c(4,5,6,7)) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))+
  scale_x_discrete("Treat", labels = c('B_pumilus_45_39_Ggt_before' = "B. pumilus 45.39 before",
                                          'P_brassicacearum_37_15_Ggt_before' = "P. brassicacearum 37.15 before",
                                          'P_proteolytica_68_48_Ggt_before' = "P. proteolytica 68.48 before",
                                          'S_panacis_40_16_Ggt_before' = "S. panacis 40.16 before",
                                       'B_pumilus_45_39_Ggt_after' = "B. pumilus 45.39 after",
                                       'P_brassicacearum_37_15_Ggt_after' = "P. brassicacearum 37.15 after",
                                       'P_proteolytica_68_48_Ggt_after' = "P. proteolytica 68.48 after",
                                       'S_panacis_40_16_Ggt_after' = "S. panacis 40.16 after"))
boxplot_GH11_GH13_CFU

ggsave("boxplot_GH11_GH13_CFU.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 12, height = 22, units = "cm",dpi = 200)

################ Plot with facet by Test
# Reorder levels of the Test variable
GH11_GH13_CFU$Test <- factor(GH11_GH13_CFU$Test,
                                    levels = c("CFU_before", "CFU_after"))

boxplot_GH11_GH13_CFU_facet <- ggboxplot(GH11_GH13_CFU, "Treat", "CFU", fill = "Treat",
                                     palette = c("#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a"),
                                     # order = c("Silty_W1_Control", "Silty_W1_Ggt",
                                     #           "Silty_WM_Control", "Silty_WM_Ggt",
                                     #           "Sandy_W1_Control", "Sandy_W1_Ggt",
                                     #           "Sandy_W3_Control", "Sandy_W3_Ggt"),
                                     ylab = "CFU (log)",
                                     bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  
  #geom_boxplot(alpha = 0.1) +
  facet_grid(~Test, scales = "free", space="free_x") +
  # scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 2.5), label = c("0", "2.5", "5.0")) +
  #stat_compare_means(method = "wilcox.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c(5,5,5,5,5,5,5,5,5,5)) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  scale_x_discrete("Treat", labels = c('B_pumilus_45_39_Ggt_before' = "B. pumilus 45.39",
                                       'P_brassicacearum_37_15_Ggt_before' = "P. brassicacearum 37.15",
                                       'P_proteolytica_68_48_Ggt_before' = "P. proteolytica 68.48",
                                       'S_panacis_40_16_Ggt_before' = "S. panacis 40.16",
                                       'B_pumilus_45_39_Ggt_after' = "B. pumilus 45.39",
                                       'P_brassicacearum_37_15_Ggt_after' = "P. brassicacearum 37.15",
                                       'P_proteolytica_68_48_Ggt_after' = "P. proteolytica 68.48",
                                       'S_panacis_40_16_Ggt_after' = "S. panacis 40.16"))
boxplot_GH11_GH13_CFU_facet

ggsave("boxplot_GH11_GH13_CFU_facet.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 18, height = 22, units = "cm",dpi = 200)



# GH11_GH12_GH13
########### CFU with non-autoclaved soil
GH11_GH12_GH13_CFU <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/CFU_GH11_GH12_GH13.csv")


##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypothesis (the population is normally distributed) is rejected
hist(GH11_GH12_GH13_CFU$CFU)
shapiro.test(GH11_GH12_GH13_CFU$CFU)
qqnorm(GH11_GH12_GH13_CFU$CFU)

compare_treats <- list(c("B_pumilus_45_39_Ggt_before", "B_pumilus_45_39_Ggt_after"),
                       c("P_brassicacearum_37_15_Ggt_before", "P_brassicacearum_37_15_Ggt_after"),
                       c("P_proteolytica_68_48_Ggt_before", "P_proteolytica_68_48_Ggt_after"),
                       c("S_panacis_40_16_Ggt_before", "S_panacis_40_16_Ggt_after"))


################ Plot with statistics
boxplot_GH11_GH12_GH13_CFU <- ggboxplot(GH11_GH12_GH13_CFU, "Treat", "CFU", fill = "Treat",
                                   palette = c("#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a"),
                                   ylab = "CFU (log)",
                                   bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +
  #scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5), label = c("0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
  #stat_compare_means(method = "wilcox.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c(4,5,6,7)) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))+
  scale_x_discrete("Treat", labels = c('B_pumilus_45_39_Ggt_before' = "B. pumilus 45.39 before",
                                       'P_brassicacearum_37_15_Ggt_before' = "P. brassicacearum 37.15 before",
                                       'P_proteolytica_68_48_Ggt_before' = "P. proteolytica 68.48 before",
                                       'S_panacis_40_16_Ggt_before' = "S. panacis 40.16 before",
                                       'B_pumilus_45_39_Ggt_after' = "B. pumilus 45.39 after",
                                       'P_brassicacearum_37_15_Ggt_after' = "P. brassicacearum 37.15 after",
                                       'P_proteolytica_68_48_Ggt_after' = "P. proteolytica 68.48 after",
                                       'S_panacis_40_16_Ggt_after' = "S. panacis 40.16 after"))
boxplot_GH11_GH12_GH13_CFU

ggsave("boxplot_GH11_GH12_GH13_CFU.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 12, height = 22, units = "cm",dpi = 200)

################ Plot with facet by Test
# Reorder levels of the Test variable
GH11_GH12_GH13_CFU$Test <- factor(GH11_GH12_GH13_CFU$Test,
                             levels = c("CFU_before", "CFU_after"))

boxplot_GH11_GH12_GH13_CFU_facet <- ggboxplot(GH11_GH12_GH13_CFU, "Treat", "CFU", fill = "Treat",
                                         palette = c("#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a","#ebaf7a", "#ebaf7a"),
                                         # order = c("Silty_W1_Control", "Silty_W1_Ggt",
                                         #           "Silty_WM_Control", "Silty_WM_Ggt",
                                         #           "Sandy_W1_Control", "Sandy_W1_Ggt",
                                         #           "Sandy_W3_Control", "Sandy_W3_Ggt"),
                                         ylab = "CFU (log)",
                                         bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  
  #geom_boxplot(alpha = 0.1) +
  facet_grid(~Test, scales = "free", space="free_x") +
  # scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 2.5), label = c("0", "2.5", "5.0")) +
  #stat_compare_means(method = "wilcox.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c(5,5,5,5,5,5,5,5,5,5)) +
  theme(legend.position = 'none') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black")) +
  scale_x_discrete("Treat", labels = c('B_pumilus_45_39_Ggt_before' = "B. pumilus 45.39",
                                       'P_brassicacearum_37_15_Ggt_before' = "P. brassicacearum 37.15",
                                       'P_proteolytica_68_48_Ggt_before' = "P. proteolytica 68.48",
                                       'S_panacis_40_16_Ggt_before' = "S. panacis 40.16",
                                       'B_pumilus_45_39_Ggt_after' = "B. pumilus 45.39",
                                       'P_brassicacearum_37_15_Ggt_after' = "P. brassicacearum 37.15",
                                       'P_proteolytica_68_48_Ggt_after' = "P. proteolytica 68.48",
                                       'S_panacis_40_16_Ggt_after' = "S. panacis 40.16"))
boxplot_GH11_GH12_GH13_CFU_facet

ggsave("boxplot_GH11_GH12_GH13_CFU_facet.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 18, height = 22, units = "cm",dpi = 200)





