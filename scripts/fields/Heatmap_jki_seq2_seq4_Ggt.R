#Load packages
library("readxl")
library("ggplot2")
library("ggpubr")
library("dplyr")
#Site_colors <- c("Go" = "#1d5b65", "Ki" = "#A34828")

################ Load table.  --- Two YEARS

jki_seq2_seq4_Ggt <- read_excel("~/Documents/R_analysis/jki_seq2_seq4/data_jki_seq2_seq4/jki_seq2_seq4_qPCR_Ggt.xlsx")
jki_seq2_seq4_Ggt
str(jki_seq2_seq4_Ggt)

# Select only Go samples
jki_seq2_seq4_Ggt_Go <- jki_seq2_seq4_Ggt[jki_seq2_seq4_Ggt$Site == "Go", ]

#Convert variables on the x-axis to factors with desired order
jki_seq2_seq4_Ggt_Go$Microhabitat <- factor(jki_seq2_seq4_Ggt_Go$Microhabitat, levels = c("RP", "RH", "BS"))

heatmap_jki_seq2_seq4_Ggt_Go<-ggplot(data=jki_seq2_seq4_Ggt_Go,mapping=aes(x=Rotation, y=Microhabitat, fill=Ggt_ng_DNA_g_soil))
heatmap_jki_seq2_seq4_Ggt_Go+ 
  geom_tile(colour="black",size=0.25) +
  facet_grid(Year~Time, scales="free", space="free")+
  xlab(label="") +
  ylab(label="")+
  scale_fill_gradient(limits=c(0, 250), breaks=seq(0,250,by=100), name="Ggt (ng/g)",low="white", high="#1d5b65", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size=22, colour="black")) +
  theme(axis.text.y = element_text(size=22, colour="black")) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  theme(strip.text.x = element_text(size = 22, colour = "black", angle = 0)) +
  theme(strip.text.y = element_text(size = 22, colour = "black", angle = 90)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7), 
        strip.background = element_rect(color = "black", size = 0.7)) +
  #theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="grey", color="grey")) +
  ggtitle(label="")

ggsave("heatmap_jki_seq2_seq4_Ggt_Go.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4", width = 16, height = 13, units = "cm",dpi = 300)


# Select only Ki samples
jki_seq2_seq4_Ggt_Ki <- jki_seq2_seq4_Ggt[jki_seq2_seq4_Ggt$Site == "Ki", ]

#Convert variables on the x-axis to factors with desired order
jki_seq2_seq4_Ggt_Ki$Microhabitat <- factor(jki_seq2_seq4_Ggt_Ki$Microhabitat, levels = c("RP", "RH", "BS"))

heatmap_jki_seq2_seq4_Ggt_Ki<-ggplot(data=jki_seq2_seq4_Ggt_Ki,mapping=aes(x=Rotation, y=Microhabitat, fill=Ggt_ng_DNA_g_soil))
heatmap_jki_seq2_seq4_Ggt_Ki+ 
  geom_tile(colour="black",size=0.25) +
  facet_grid(Year~Time, scales="free", space="free")+
  xlab(label="") +
  ylab(label="")+
  scale_fill_gradient(limits=c(0, 250), breaks=seq(0,250,by=100), name="Ggt (ng/g)",low="white", high="#A34828", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size=22, colour="black")) +
  theme(axis.text.y = element_text(size=22, colour="black")) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  theme(strip.text.x = element_text(size = 22, colour = "black", angle = 0)) +
  theme(strip.text.y = element_text(size = 22, colour = "black", angle = 90)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7), 
        strip.background = element_rect(color = "black", size = 0.7)) +
  #theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="grey", color="grey")) +
  ggtitle(label="")

ggsave("heatmap_jki_seq2_seq4_Ggt_Ki.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4", width = 14, height = 13, units = "cm",dpi = 300)



################ Load table.  2021

jki_seq2_seq4_Ggt_2021 <- read_excel("~/Documents/R_analysis/jki_seq2_seq4/data_jki_seq2_seq4/jki_seq2_seq4_qPCR_Ggt_2021.xlsx")
jki_seq2_seq4_Ggt_2021
str(jki_seq2_seq4_Ggt_2021)

# Select only Go samples
jki_seq2_seq4_Ggt_2021_Go <- jki_seq2_seq4_Ggt_2021[jki_seq2_seq4_Ggt_2021$Site == "Go", ]

#Convert variables on the x-axis to factors with desired order
jki_seq2_seq4_Ggt_2021_Go$Microhabitat <- factor(jki_seq2_seq4_Ggt_2021_Go$Microhabitat, levels = c("RP", "RH", "BS"))

heatmap_jki_seq2_seq4_Ggt_2021_Go<-ggplot(data=jki_seq2_seq4_Ggt_2021_Go,mapping=aes(x=Rotation, y=Microhabitat, fill=Ggt_ng_DNA_g_soil))
heatmap_jki_seq2_seq4_Ggt_2021_Go+ 
  geom_tile(colour="black",size=0.25) +
  facet_grid(~Time, scales="free", space="free")+
  xlab(label="") +
  ylab(label="")+
  scale_fill_gradient(limits=c(0, 50), breaks=seq(0,50,by=25), name="Ggt_2021 (ng/g)",low="white", high="#1d5b65", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size=22, colour="black")) +
  theme(axis.text.y = element_text(size=22, colour="black")) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  theme(strip.text.x = element_text(size = 22, colour = "black", angle = 0)) +
  theme(strip.text.y = element_text(size = 22, colour = "black", angle = 90)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7), 
        strip.background = element_rect(color = "black", size = 0.7)) +
  #theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="grey", color="grey")) +
  ggtitle(label="")

ggsave("heatmap_jki_seq2_seq4_Ggt_2021_Go.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4", width = 14, height = 8, units = "cm",dpi = 300)


# Select only Ki samples
jki_seq2_seq4_Ggt_2021_Ki <- jki_seq2_seq4_Ggt_2021[jki_seq2_seq4_Ggt_2021$Site == "Ki", ]

#Convert variables on the x-axis to factors with desired order
jki_seq2_seq4_Ggt_2021_Ki$Microhabitat <- factor(jki_seq2_seq4_Ggt_2021_Ki$Microhabitat, levels = c("RP", "RH", "BS"))

heatmap_jki_seq2_seq4_Ggt_2021_Ki<-ggplot(data=jki_seq2_seq4_Ggt_2021_Ki,mapping=aes(x=Rotation, y=Microhabitat, fill=Ggt_ng_DNA_g_soil))
heatmap_jki_seq2_seq4_Ggt_2021_Ki+ 
  geom_tile(colour="black",size=0.25) +
  facet_grid(~Time, scales="free", space="free")+
  xlab(label="") +
  ylab(label="")+
  scale_fill_gradient(limits=c(0, 50), breaks=seq(0,50,by=25), name="Ggt_2021 (ng/g)",low="white", high="#A34828", na.value="lightgrey") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size=22, colour="black")) +
  theme(axis.text.y = element_text(size=22, colour="black")) +
  theme(legend.title = element_text(size=22), legend.text = element_text(size = 22, colour = "black")) +
  theme(strip.text.x = element_text(size = 22, colour = "black", angle = 0)) +
  theme(strip.text.y = element_text(size = 22, colour = "black", angle = 90)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7), 
        strip.background = element_rect(color = "black", size = 0.7)) +
  #theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="grey", color="grey")) +
  ggtitle(label="")

ggsave("heatmap_jki_seq2_seq4_Ggt_2021_Ki.png", path = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Heatmap_jki_seq2_seq4", width = 14, height = 8, units = "cm",dpi = 300)


## The end! Enjoy!