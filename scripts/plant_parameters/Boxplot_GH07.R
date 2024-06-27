#Load packages
library("readxl")
library("ggplot2")
library("ggpubr")

########### Shoot_FW with non-autoclaved soil
GH07_results <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/GH07_results.csv")


##### SHOOT FW
##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(GH07_results$Shoot_FW)
shapiro.test(GH07_results$Shoot_FW)
qqnorm(GH07_results$Shoot_FW)

################ 
boxplot_Shoot_FW_GH07 <- ggboxplot(GH07_results, "Treatment", "Shoot_FW", fill = "#80c698",
                                   color = "Treatment", palette = c("#80c698", "#80c698","#80c698","#80c698","#80c698","#80c698",
                                                                    "#80c698","#80c698","#80c698","#80c698","#80c698"),         
                                  ylab = "Shoot_FW",
                                  bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  #facet_grid(~Site, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0, 0.32), breaks = seq(0, 0.32, by = 0.1), label = c("0", "0.1", "0.2", "0.3")) +
  #stat_compare_means(test = "kruskal.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c("4.0","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9")) +
  #stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "Negative", label.y = 0.31) +
  #stat_compare_means(label.y = 1) + 
  scale_x_discrete("Treatment", labels = c('B_myc' = "B. myc",
                                           'B_pum' = "B. pum",
                                           'E_mor' = "E. mor",
                                           'Pae_kri' = "Pae. kri",
                                           'P_bra' = "P. bra",
                                           'P_ext' = "P. ext",
                                           'P_lin' = "P. lin",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan",
                                           'B_vel' = "B. vel")) + 
  ylab("Shoot fresh weight (g)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_Shoot_FW_GH07

ggsave("boxplot_Shoot_FW_GH07.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 20, height = 16, units = "cm",dpi = 200)


##### ROOT FW
#check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(GH07_results$Root_FW)
shapiro.test(GH07_results$Root_FW)
qqnorm(GH07_results$Root_FW)

################ 
boxplot_Root_FW_GH07 <- ggboxplot(GH07_results, "Treatment", "Root_FW", fill = "#c6ba9c",
                                   color = "Treatment", palette = c("#c6ba9c", "#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c",
                                                                    "#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c"),         
                                   ylab = "Root_FW",
                                   bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  #facet_grid(~Site, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1), label = c("0", "0.1", "0.2", "0.3", "0.4", "0.5")) +
  #stat_compare_means(test = "kruskal.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c("4.0","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9")) +
  #stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "Negative", label.y = 0.49) +
  #stat_compare_means(label.y = 1) + 
  scale_x_discrete("Treatment", labels = c('B_myc' = "B. myc",
                                           'B_pum' = "B. pum",
                                           'E_mor' = "E. mor",
                                           'Pae_kri' = "Pae. kri",
                                           'P_bra' = "P. bra",
                                           'P_ext' = "P. ext",
                                           'P_lin' = "P. lin",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan",
                                           'B_vel' = "B. vel")) + 
  ylab("Root fresh weight (g)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_Root_FW_GH07

ggsave("boxplot_Root_FW_GH07.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 20, height = 16, units = "cm",dpi = 200)
