#load libraries
library("ggplot2")
library("tidyr")
library("readxl")
library("dplyr")
library("RColorBrewer")

#load data 
Selected9<-read_excel("~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/data/jki_seq2_seq4_isolates_selected9.xlsx")
Selected9

###Preparing tidyr format table (variables in columns, observation in rows, values in cells), collecting a set of column names and place them into a single value column
Selected9_processed<-gather(data=Selected9,key = Tests, value = Result,Cellulase:PRND,-ID,-Site,-Original_medium,-Sample,-Isolate, -Rotation,-Replicate,-Compartment,-SCM, -SC, -CM, -Closest_match, -Isolate_id, -Genus_id)
head(Selected9_processed)

###Convert data to square root 
#vt_long$Sqrt.abundance<-sqrt(vt.long$Absolute_Abundance)

#####Constructing the ggplot heatmap

### Heatmap GO_96  
#Convert variables on the x-axis to factors with desired order
Selected9_processed$Tests <- factor(Selected9_processed$Tests, levels = c("Ggt", "Fusarium", "Rhizoctonia", "Cellulase", "Chitinase", "Glucanase", "Protease", "ACC", "IAA", "PO4", "Siderophore", "DAPG", "PCA", "PRND"))
Selected9_processed$Compartment <- factor(Selected9_processed$Compartment, levels = c("BS", "RH", "RP"))

#Convert variables on the x-axis to factors with desired order
Selected9_processed$Genus_id <- factor(Selected9_processed$Genus_id, levels = rev(sort(unique(Selected9_processed$Genus_id))))

heatmap_selected9<-ggplot(data=Selected9_processed,mapping=aes(x=Tests, y=Genus_id, fill=Result))
heatmap_selected9 + 
  geom_tile(colour="black",size=0.25) +
  facet_grid(Compartment~Tests, scales="free", space="free") +
  xlab(label="") +
  scale_fill_gradient(name="Index",low="white", high="brown", na.value="#e3dddc") +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=16, colour="black", face = "italic")) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 16, colour = "black")) +
  theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90, hjust=0)) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 0)) +
    theme(strip.placement = "top", plot.title = element_text(hjust=0.5), axis.title.y = element_blank(), strip.background=element_rect(fill="#EEEEEE", color="#424242")) +
  ggtitle(label="")

ggsave("heatmap_selected9.png", path = "~/Documents/R_analysis/jki_seq2_seq4/isol_jki_seq2_seq4/output/", width = 25, height = 14, units = "cm",dpi = 300)

