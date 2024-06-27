#Load packages
library("readxl")
library("ggplot2")
library("ggpubr")

########### Shoot_FW with non-autoclaved soil
GH13_results <- read.csv ("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/GH13_results.csv")


##### SHOOT FW
##check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(GH13_results$Shoot_FW)
shapiro.test(GH13_results$Shoot_FW)
qqnorm(GH13_results$Shoot_FW)

################ 
boxplot_Shoot_FW_GH13 <- ggboxplot(GH13_results, "Treatment", "Shoot_FW", fill = "#80c698",
                                   color = "Treatment", palette = c("#80c698", "#80c698","#80c698","#80c698","#80c698","#80c698",
                                                                    "#80c698","#80c698","#80c698","#80c698","#80c698"),         
                                  ylab = "Shoot_FW",
                                  bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  #facet_grid(~Site, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, by = 0.1), label = c("0", "0.1", "0.2", "0.3")) +
  #stat_compare_means(test = "kruskal.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c("4.0","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9")) +
  #stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "Positive", label.y = 0.24) +
  #stat_compare_means(label.y = 1) + 
  scale_x_discrete("Treatment", labels = c('B_pum' = "B. pum",
                                           'P_bra' = "P. bra",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan",
                                           'Consortium_1' = "Consortium 1",
                                           'Consortium_2' = "Consortium 2")) + 
  ylab("Shoot fresh weight (g)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_Shoot_FW_GH13

ggsave("boxplot_Shoot_FW_GH13.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 12, height = 16, units = "cm",dpi = 200)


##### ROOT FW
#check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(GH13_results$Root_FW)
shapiro.test(GH13_results$Root_FW)
qqnorm(GH13_results$Root_FW)

################ 
boxplot_Root_FW_GH13 <- ggboxplot(GH13_results, "Treatment", "Root_FW", fill = "#c6ba9c",
                                   color = "Treatment", palette = c("#c6ba9c", "#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c",
                                                                    "#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c"),         
                                   ylab = "Root_FW",
                                   bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  #facet_grid(~Site, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0, 0.4), breaks = seq(0, 0.4, by = 0.1), label = c("0", "0.1", "0.2", "0.3", "0.4")) +
  #stat_compare_means(test = "kruskal.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c("4.0","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9")) +
  #theme(axis.ticks.x = element_blank()) +stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "Positive", label.y = 0.4) +
  #stat_compare_means(label.y = 1) + 
  scale_x_discrete("Treatment", labels = c('B_pum' = "B. pum",
                                           'P_bra' = "P. bra",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan",
                                           'Consortium_1' = "Consortium 1",
                                           'Consortium_2' = "Consortium 2")) + 
  ylab("Root fresh weight (g)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_Root_FW_GH13

ggsave("boxplot_Root_FW_GH13.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 12, height = 16, units = "cm",dpi = 200)


##### Disease
#check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(GH13_results$Disease)
shapiro.test(GH13_results$Disease)
qqnorm(GH13_results$Disease)

################ 
boxplot_Disease_GH13 <- ggboxplot(GH13_results, "Treatment", "Disease", fill = "#f67280",
                                  color = "Treatment", palette = c("#f67280", "#f67280","#f67280","#f67280","#f67280","#f67280",
                                                                   "#f67280","#f67280","#f67280","#f67280","#f67280"),         
                                  ylab = "Disease",
                                  bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  #facet_grid(~Site, scales = "free", space="free_x") +
  scale_y_continuous(limits = c(0, 4.3), breaks = seq(0, 4.3, by = 1), label = c("0", "1", "2", "3", "4")) +
  #stat_compare_means(test = "kruskal.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c("4.0","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9")) +
  #stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "Positive", label.y = 4.2) +
  #stat_compare_means(label.y = 1) + 
  scale_x_discrete("Treatment", labels = c('B_pum' = "B. pum",
                                           'P_bra' = "P. bra",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan",
                                           'Consortium_1' = "Consortium 1",
                                           'Consortium_2' = "Consortium 2")) + 
  ylab("Disease severity index") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_Disease_GH13

ggsave("boxplot_Disease_GH13.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 11.5, height = 16, units = "cm",dpi = 200)


##### ROOT WR
#check for the distribution 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05),
# then the null hypotesis (the population is normally distributed) is rejected
hist(GH13_results$Root_WR)
shapiro.test(GH13_results$Root_WR)
qqnorm(GH13_results$Root_WR)

################ 
boxplot_Root_WR_GH13 <- ggboxplot(GH13_results, "Treatment", "Root_WR", fill = "#c6ba9c",
                                  color = "Treatment", palette = c("#c6ba9c", "#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c",
                                                                   "#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c","#c6ba9c"),         
                                  ylab = "Root_WR",
                                  bxp.errorbar = TRUE, bxp.errorbar.width = 0.2) +
  theme_bw() +  geom_boxplot(alpha = 0.2) +
  #facet_grid(~Site, scales = "free", space="free_x") +
  #scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5), label = c("0", "0.5", "1.0", "1.5", "2.0", "2.5")) +
  #stat_compare_means(test = "kruskal.test", comparisons = compare_treats, label= "p", bracket.size = .3, size=4, label.y = c("4.0","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9")) +
 # stat_compare_means(aes(label = after_stat(p.format)), method = "wilcox.test", ref.group = "Positive", label.y = 200) +
  #stat_compare_means(label.y = 1) + 
  scale_x_discrete("Treatment", labels = c('B_pum' = "B. pum",
                                           'P_bra' = "P. bra",
                                           'P_pro' = "P. pro",
                                           'S_pan' = "S. pan",
                                           'Consortium_1' = "Consortium 1",
                                           'Consortium_2' = "Consortium 2")) + 
  ylab("Root length (cm)") +
  theme(legend.position = 'none') +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=24, colour = "black")) +
  theme(axis.text.x = element_text(size=24, angle = 90,vjust =0.5, hjust = 1), axis.text.y = element_text(size=24)) +
  theme(strip.text.x = element_text(size=24), strip.text.y = element_text(size=24)) +
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24, colour = "black"))
boxplot_Root_WR_GH13

ggsave("boxplot_Root_WR_GH13.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 13, height = 16, units = "cm",dpi = 200)

