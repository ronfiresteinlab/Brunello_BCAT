require(VennDiagram);
require(dplyr)
require(magrittr)
require("openxlsx")

df <- read.xlsx("./Brunello_lowGFP.xlsx")

df_epick <- read.xlsx("./Table S2.xlsx", sheet = 2, startRow = 2)

# Figure 3A
venn.diagram(
  x = list(
    I = df_epick %>% 
      filter(DLD1_MYC_low_score > 0 & DLD1_MYC_low_pvalue < 0.1) %>% 
      .$Symbol,
    IV = df_epick %>% 
      filter(RKO_MYC_low_score >0 & RKO_MYC_low_pvalue < 0.1) %>% 
      .$Symbol,
    III = df_epick %>% 
      filter(DLD1_AXIN2_low_score > 0 & DLD1_AXIN2_low_pvalue < 0.1) %>% 
      .$Symbol,
    II = df_epick %>% 
      filter(DLD1_TOPGC_low_score > 0 & DLD1_TOPGC_low_pvalue < 0.1) %>% 
      .$Symbol
  ),
  filename = "3A.tiff",
  col = "black",
  lty = "dotted",
  lwd = 4,
  fill = c("cornflowerblue", "green", "pink", "darkorchid1"),
  alpha = 0.50,
  label.col = c("blue", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 2.5,
  cat.fontfamily = "serif"
);

# Figure 3B
venn.diagram(
  x = list(
    I = df %>% 
      filter(DLD1_MYC_low_score > 0 & DLD1_MYC_low_pvalue < 0.1) %>% 
      .$Symbol,
    IV = df %>% 
      filter(RKO_MYC_low_score >0 & RKO_MYC_low_pvalue < 0.1) %>% 
      .$Symbol,
    III = df %>% 
      filter(DLD1_AXIN2_low_score > 0 & DLD1_AXIN2_low_pvalue < 0.1) %>% 
      .$Symbol,
    II = df %>% 
      filter(DLD1_TOPGC_low_score > 0 & DLD1_TOPGC_low_pvalue < 0.1) %>% 
      .$Symbol
  ),
  filename = "3B.tiff",
  col = "black",
  lty = "dotted",
  lwd = 4,
  fill = c("cornflowerblue", "green", "pink", "darkorchid1"),
  alpha = 0.50,
  label.col = c("blue", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 2.5,
  cat.fontfamily = "serif"
)

