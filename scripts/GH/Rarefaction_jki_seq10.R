#Create rarefied data from phyloseq object using vegan package or microbiome package
#Loading package
library("phyloseq")
library("vegan")
library("microbiome")
library("ggpubr")
library("ggplot2")
library("tidyr")

##Creating rarefaction curve of non-rarefied samples
#transpose S4 otu-table
mat_psO_jki_seq10 <- t(otu_table(psO_jki_seq10))

#Transform S4 objet to matrix (a warning message is normal)
class(mat_psO_jki_seq10) <- "matrix"

#Create curve using rarecurve()
system.time(rarecurve(mat_psO_jki_seq10, step = 1000, lwd=2, ylab="OTU", col = "black", label = FALSE))

##Rarefying samples to the minimum number of reads among samples
#Using rarefy_even_depth() rarefy to the lower number of total sequences in a sample
psO_jki_seq10_rarefied<-rarefy_even_depth(psO_jki_seq10, rngseed=2022, sample.size = min(sample_sums(psO_jki_seq10)),trimOTUs=TRUE)
psO_jki_seq10_rarefied
sample_sums(psO_jki_seq10_rarefied)

##Creating rarefaction curve of non-rarefied samples
#transpose S4 otu-table
mat_psO_jki_seq10_rarefied <- t(otu_table(psO_jki_seq10_rarefied))

#Transform S4 objet to matrix
class(mat_psO_jki_seq10_rarefied) <- "matrix"

#Create curve using rarecurve()
system.time(rarecurve(mat_psO_jki_seq10_rarefied, step = 1000, lwd=2, ylab="OTU", col = "blue", label = FALSE))

#Prepare file with meta data using meta function
psO_jki_seq10_rarefied_rare.meta <- meta(psO_jki_seq10_rarefied)
head(psO_jki_seq10_rarefied_rare.meta)
psO_jki_seq10_rarefied_rare.meta.Site_rot <- subset(psO_jki_seq10_rarefied_rare.meta, select = -c(Sample_id, Site, Rotation, Replicate, Microhabitat, Treatment_id, Site_rot_tre, Treatment))
head(psO_jki_seq10_rarefied_rare.meta.Site_rot)

#Create and transpose matrix to make x axis for samples
df_psO_jki_seq10_rarefied <- data.frame(otu_table(psO_jki_seq10_rarefied))
#write.csv(df_psO_jki_seq10_rarefied, "~/Documents/R_analysis/jki_seq10/output_jki_seq10/Tables_jki_seq10/df_psO_jki_seq10_rarefied.csv")

df_psO_jki_seq10_rarefied_t <- t(df_psO_jki_seq10_rarefied)

#Combine matrix in one
ASVs_rarefied_ed<-cbind(psO_jki_seq10_rarefied_rare.meta.Site_rot, df_psO_jki_seq10_rarefied_t)
ASVs_rarefied_ed[1:30,1]

##check for the distribution of any ASV counts in the sample group using "hist" function to plot, 
#"shapiro.test" function to test the Null hypothesis, and "qqnorm" function to qq plot
# If alpha result is lower than the alpha value chosen (p<0.05) 
#then the null hypotesis (the population is normally distributed) is rejected
hist(ASVs_rarefied_ed$sp1)
shapiro.test(ASVs_rarefied_ed$sp1)
qqnorm(ASVs_rarefied_ed$sp1)

#Go = "#1d5b65"   Go 2020 = "#216974"  Go 2021 = "#6baea5"
#Ki = "#A34828"   Ki 2020 = "#D1711F"  Ki 2021 = "#ebaf7a"

#Define colours
colour = rep(NA, length=length(ASVs_rarefied_ed[,1]))
colour[which(ASVs_rarefied_ed$Site_rot=="Go_W1")] = "#216974"
colour[which(ASVs_rarefied_ed$Site_rot=="Go_WM")] = "#6baea5"
colour[which(ASVs_rarefied_ed$Site_rot=="Ki_W1")] = "#D1711F"
colour[which(ASVs_rarefied_ed$Site_rot=="Ki_W3")] = "#ebaf7a"

dim(ASVs_rarefied_ed)

#plot
tiff("~/Documents/R_analysis/jki_seq10/output_jki_seq10/Alpha_div_jki_seq10/rarecurve_jki_seq10_rarefied.tiff", units="cm", width=12, height=12, res=300)
rarecurve(ASVs_rarefied_ed [,2:18617], step=1000, label=FALSE, col=colour, main="", xlab = "Number of sequences", ylab = "ASVs")
legend(legend=c("Go_W1", "Go_WM", "Ki_W1", "Ki_W3"),
       "bottomright", 
       bty="n",
       col=c("#216974", "#6baea5", "#D1711F", "#ebaf7a"),
       pch=15,
       pt.cex=1.5,
       cex=0.75,
       ncol=3)
#arrows(x0=30464, y0=300, y1=1800, angle=90, length=0, lty = 2)
dev.off()


## The end! Have fun! Enjoy the day!  : )