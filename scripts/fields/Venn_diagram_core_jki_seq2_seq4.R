#Create a venndiagram 

#loading packages
library("limma")
library("VennDiagram")
library("ggplot2")

#Venndiagram Go 2020 BS (triple)

venn_Go_2020_BS <- draw.triple.venn(
  area1 = 248,
  area2 = 247,
  area3 = 205,
  n12 = 190,
  n13 = 173,
  n23 = 179,
  n123 = 162,
  category = c("W1", "W2", "WM"),
  fill = c("#e7e5cc", "#9cc184", "#1e3d14"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#1e3d14", "#1e3d14", "#1e3d14")
);

tiff(filename = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Venn_jki_seq2_seq4/venn_Go_2020_BS.tiff", width = 14, height = 14, units = "cm", res = 200, compression = "lzw");
grid.draw(venn_Go_2020_BS);
dev.off()


#Venndiagram Go 2020 RH (triple)

venn_Go_2020_RH <- draw.triple.venn(
  area1 = 239,
  area2 = 229,
  area3 = 216,
  n12 = 184,
  n13 = 184,
  n23 = 179,
  n123 = 167,
  category = c("W1", "W2", "WM"),
  fill = c("#e7e5cc", "#9cc184", "#1e3d14"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#1e3d14", "#1e3d14", "#1e3d14")
);

tiff(filename = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Venn_jki_seq2_seq4/venn_Go_2020_RH.tiff", width = 14, height = 14, units = "cm", res = 200, compression = "lzw");
grid.draw(venn_Go_2020_RH);
dev.off()


#Venndiagram Go 2020 RP (triple)

venn_Go_2020_RP <- draw.triple.venn(
  area1 = 233,
  area2 = 235,
  area3 = 219,
  n12 = 191, #tem W1 e W2 juntos somente
  n13 = 186,
  n23 = 187,
  n123 = 171,
  category = c("W1", "W2", "WM"),
  fill = c("#e7e5cc", "#9cc184", "#1e3d14"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#1e3d14", "#1e3d14", "#1e3d14")
);

tiff(filename = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Venn_jki_seq2_seq4/venn_Go_2020_RP.tiff", width = 14, height = 14, units = "cm", res = 200, compression = "lzw");
grid.draw(venn_Go_2020_RP);
dev.off()


# draw.pairwise.venn
venn_Ki_2020_BS  <- draw.pairwise.venn(
  area1 = 199,
  area2 = 204,
  cross.area = 168,
  category = c("W1", "W3"),
  fill = c("#e7e5cc", "#447243"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#1e3d14", "#1e3d14")
);

tiff(filename = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Venn_jki_seq2_seq4/venn_Ki_2020_BS.tiff", width = 14, height = 14, units = "cm", res = 200, compression = "lzw");
grid.draw(venn_Ki_2020_BS);
dev.off()


venn_Ki_2020_RH  <- draw.pairwise.venn(
  area1 = 210,
  area2 = 229,
  cross.area = 185,
  category = c("W1", "W3"),
  fill = c("#e7e5cc", "#447243"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#1e3d14", "#1e3d14")
);

tiff(filename = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Venn_jki_seq2_seq4/venn_Ki_2020_RH.tiff", width = 14, height = 14, units = "cm", res = 200, compression = "lzw");
grid.draw(venn_Ki_2020_RH);
dev.off()


venn_Ki_2020_RP  <- draw.pairwise.venn(
  area1 = 207,
  area2 = 208,
  cross.area = 175,
  category = c("W1", "W3"),
  fill = c("#e7e5cc", "#447243"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("#1e3d14", "#1e3d14")
);

tiff(filename = "~/Documents/R_analysis/jki_seq2_seq4/output_jki_seq2_seq4/Venn_jki_seq2_seq4/venn_Ki_2020_RP.tiff", width = 14, height = 14, units = "cm", res = 200, compression = "lzw");
grid.draw(venn_Ki_2020_RP);
dev.off()

