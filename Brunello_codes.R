require("openxlsx")
require(tidyverse)
require(ggrepel)

##load data
#Brunello sgRNA fold-change data
brunello <- read.xlsx("E:/R/brunello_library_sum.xlsx")


#Brunello MAGEcK datasets, lowGFP, highGFP and proliferation
Enrichment_list <- read.xlsx("E:/R/Brunello_lowGFP.xlsx")

Enrichment_highGFP_list <- read.xlsx("E:/R/Brunello_HiGFP.xlsx")

brunello_proliferation <- read.xlsx("E:/R/Proliferation_D21_plasmid.xlsx")

#depmap data downloaded on 24-July-2019
Essential_gene <- read.csv("E:/R/Achilles_common_essentials.csv", header = T, stringsAsFactors = F)
DLD1_RKO <- read.xlsx("E:/R/DLD1_RKO_expression.xlsx")

# for Table S1
df1 <- Enrichment_list %>% 
  left_join(Enrichment_highGFP_list, by = "Symbol") %>% 
  left_join(brunello_proliferation, by = "Symbol") %>% 
  #left_join(DLD1_RKO, by = "Symbol") %>% 
  left_join(Essential_gene, by = "Symbol")

write.csv(df1, "E:/R/Table S1.csv")

##make enrichment plots
#median of LFC for all sgRNAs per gene
brunello1 <- brunello %>% 
  ##calculate log2(Fold_change) value
  mutate(TOPGC_FC = log2(DLD1_TOPGC),
         AXIN2_FC = log2(DLD1_Axin2),
         MYC_FC = log2(DLD1_MYC),
         RKO_FC = log2(RKO_MYC),
         TOPGC_highFC = log2(DLD1_TOPGC_HighGFP),
         AXIN2_highFC = log2(DLD1_Axin2_HighGFP),
         MYC_highFC = log2(DLD1_MYC_HighGFP),
         RKO_highFC = log2(RKO_MYC_HighGFP),
         DLD1_Proliferation = log2(DLD1_proliferation_FC),
         RKO_proliferation = log2(RKO_proliferaiton_FC)
  ) %>% 
  filter_all(all_vars(!is.infinite(.))) %>% 
  group_by(Symbol) %>% 
  ##calculate median log2(Fold_change) per gene
  transmute(TOPGC_median_FC = median(TOPGC_FC),
         AXIN2_median_FC = median(AXIN2_FC),
         MYC_median_FC = median(MYC_FC),
         RKO_median_FC = median(RKO_FC),
         TOPGC_median_highFC = median(TOPGC_highFC),
         AXIN2_median_highFC = median(AXIN2_highFC),
         MYC_median_highFC = median(MYC_highFC),
         RKO_median_highFC = median(RKO_highFC),
         DLD1_median_Proliferation = median(DLD1_Proliferation),
         RKO_median_proliferation = median(RKO_proliferation)
         ) %>% 
  unique(.) %>% 
##calculate Zscore for every gene
  transmute(
         TOPGC_Zscore = (TOPGC_median_FC - mean(.$TOPGC_median_FC))/sd(.$TOPGC_median_FC),
         AXIN2_Zscore = (AXIN2_median_FC - mean(.$AXIN2_median_FC))/sd(.$AXIN2_median_FC),
         MYC_Zscore = (MYC_median_FC - mean(.$MYC_median_FC))/sd(.$MYC_median_FC),
         RKO_Zscore = (RKO_median_FC - mean(.$RKO_median_FC))/sd(.$RKO_median_FC),
         TOPGC_highZscore = (TOPGC_median_highFC - mean(.$TOPGC_median_highFC))/sd(.$TOPGC_median_highFC),
         AXIN2_highZscore = (AXIN2_median_highFC - mean(.$AXIN2_median_highFC))/sd(.$AXIN2_median_highFC),
         MYC_highZscore = (MYC_median_highFC - mean(.$MYC_median_highFC))/sd(.$MYC_median_highFC),
         RKO_highZscore = (RKO_median_highFC - mean(.$RKO_median_highFC))/sd(.$RKO_median_highFC),
         DLD1_Proliferation_Zscore = (DLD1_median_Proliferation - mean(.$DLD1_median_Proliferation))/sd(.$DLD1_median_Proliferation),
         RKO_Proliferation_Zscore = (RKO_median_proliferation - mean(.$RKO_median_proliferation))/sd(.$RKO_median_proliferation)) 


#combine two sets of data, transform mageck p-values

Brunello_mageck1  <- brunello1 %>% 
  left_join(Enrichment_list, by = "Symbol") %>% 
  #define indifinte value as 5
  mutate(TOPGC_pscore = ifelse(is.finite(log10(DLD1_TOPGC_low_pvalue)),
                               -log10(DLD1_TOPGC_low_pvalue),5),
         AXIN2_pscore = ifelse(is.finite(log10(DLD1_AXIN2_low_pvalue)),
                               -log10(DLD1_AXIN2_low_pvalue),5),
         MYC_pscore = ifelse(is.finite(log10(DLD1_MYC_low_pvalue)),
                             -log10(DLD1_MYC_low_pvalue),5),
         RKO_pscore = ifelse(is.finite(log10(RKO_MYC_low_pvalue)),
                             -log10(RKO_MYC_low_pvalue),5)
  ) %>% 
  #remove NA values from dataframe
  drop_na() %>% 
  mutate(TOPGC = ifelse(DLD1_TOPGC_low_score >= 0, "N", "T"),
         AXIN2 = ifelse(DLD1_AXIN2_low_score >= 0, "N", "T"),
         MYC = ifelse(DLD1_MYC_low_score >= 0, "N", "T"),
         RKO = ifelse(RKO_MYC_low_score >= 0, "N", "T"))  

str(Brunello_mageck1)

##make plots
#TOPGC plot
ggplot()+
  geom_point(data = Brunello_mageck1,
             aes(x = TOPGC_Zscore, y = TOPGC_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck1, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")),
             aes(x = TOPGC_Zscore, y = TOPGC_pscore, color = TOPGC),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck1, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1","RNF43")), 
                  aes(x = TOPGC_Zscore, y = TOPGC_pscore,color = TOPGC,label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0
                  )+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 1D.pdf")

#MYC plot
ggplot()+
  geom_point(data = Brunello_mageck1,
             aes(x = MYC_Zscore, y = MYC_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck1, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L")), 
             aes(x = MYC_Zscore, y = MYC_pscore, color = MYC),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck1, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L")), 
                  aes(x = MYC_Zscore, y = MYC_pscore,color = MYC,
                      label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 2C.pdf")

#Axin2 plot
ggplot()+
  geom_point(data = Brunello_mageck1,
             aes(x = AXIN2_Zscore, y = AXIN2_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck1, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L")), 
             aes(x = AXIN2_Zscore, y = AXIN2_pscore, color = AXIN2),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck1, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L")), 
                  aes(x = AXIN2_Zscore, y = AXIN2_pscore,color = AXIN2,
                      label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0
  )+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 2D.pdf")

#RKO plot
ggplot()+
  geom_point(data = Brunello_mageck1,
             aes(x = RKO_Zscore, y = RKO_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck1, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L")), 
             aes(x = RKO_Zscore, y = RKO_pscore, color = RKO),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck1, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L")), 
                  aes(x = RKO_Zscore, y = RKO_pscore,color = RKO,
                      label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 2E.pdf")

###highGFP plot

str(Enrichment_highGFP_list)
Brunello_mageck_highGFP <- Enrichment_highGFP_list %>%
  #define indifinte value as 5
  transmute(Symbol = Symbol,
         TOPGC_pscore = ifelse(is.finite(-log10(DLD1_TOPGC_high_pvalue)),
                               -log10(DLD1_TOPGC_high_pvalue),5),
         AXIN2_pscore = ifelse(is.finite(-log10(DLD1_AXIN2_high_pvalue)),
                               -log10(DLD1_AXIN2_high_pvalue),5),
         MYC_pscore = ifelse(is.finite(-log10(DLD1_MYC_high_pvalue)),
                             -log10(DLD1_MYC_high_pvalue),5),
         RKO_pscore = ifelse(is.finite(-log10(RKO_MYC_high_pvalue)),
                             -log10(RKO_MYC_high_pvalue),5)
  ) %>% 
  left_join(brunello1, by = "Symbol") %>% 
  #remove NA values from dataframe
  drop_na() %>% 
  mutate(TOPGC = ifelse(TOPGC_highZscore >= 0, "N", "T"),
         AXIN2 = ifelse(AXIN2_highZscore >= 0, "N", "T"),
         MYC = ifelse(MYC_highZscore >= 0, "N", "T"),
         RKO = ifelse(RKO_highZscore >= 0, "N", "T"))  
#TOPGC-highGFP plot
ggplot()+
  geom_point(data = Brunello_mageck_highGFP,
             aes(x = TOPGC_highZscore, y = TOPGC_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck_highGFP, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")),
             aes(x = TOPGC_highZscore, y = TOPGC_pscore, color = TOPGC),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck_highGFP, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")), 
                  aes(x = TOPGC_highZscore, y = TOPGC_pscore,color = TOPGC,label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure S1B.pdf")

#MYC-highGFP plot
ggplot()+
  geom_point(data = Brunello_mageck_highGFP,
             aes(x = MYC_highZscore, y = MYC_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck_highGFP, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")),
             aes(x = MYC_highZscore, y = MYC_pscore, color = MYC),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck_highGFP, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")), 
                  aes(x = MYC_highZscore, y = MYC_pscore,color = MYC,label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("mageck_MYChighGFP_log10pvalue_vs_Zscore.pdf")

#AXIN2-highGFP plot
ggplot()+
  geom_point(data = Brunello_mageck_highGFP,
             aes(x = AXIN2_highZscore, y = AXIN2_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck_highGFP, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")),
             aes(x = AXIN2_highZscore, y = AXIN2_pscore, color = AXIN2),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck_highGFP, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")), 
                  aes(x = AXIN2_highZscore, y = AXIN2_pscore,color = AXIN2,label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("mageck_MYChighGFP_log10pvalue_vs_Zscore.pdf")


##make proliferation plots
str(brunello_proliferation)

Brunello_mageck_Proliferation <- brunello_proliferation %>% 
  transmute(Symbol= Symbol,
         DLD1_pscore = ifelse(is.finite(-log10(DLD1_Proliferation_pvalue)),
                               -log10(DLD1_Proliferation_pvalue),5),
         RKO_pscore = ifelse(is.finite(-log10(RKO_Proliferation_pvalue)),
                               -log10(RKO_Proliferation_pvalue),5)
  ) %>% 
  left_join(brunello1, by = "Symbol") %>% 
  #remove NA values from dataframe
  drop_na() %>% 
  mutate(DLD1 = ifelse(DLD1_Proliferation_Zscore >= 0, "T", "N"),
         RKO = ifelse(RKO_Proliferation_Zscore >= 0, "T", "N")) 

str(Brunello_mageck_Proliferation)

#DLD1 Proliferation
ggplot()+
  geom_point(data = Brunello_mageck_Proliferation,
             aes(x = DLD1_Proliferation_Zscore, y = DLD1_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck_Proliferation, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")),
             aes(x = DLD1_Proliferation_Zscore, y = DLD1_pscore, color = DLD1),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck_Proliferation, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")), 
                  aes(x = DLD1_Proliferation_Zscore, y = DLD1_pscore,color = DLD1,label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure S1E.pdf")

#RKO Proliferation
ggplot()+
  geom_point(data = Brunello_mageck_Proliferation,
             aes(x = RKO_Proliferation_Zscore, y = RKO_pscore),color = "grey",size = 1)+
  geom_point(data = subset(Brunello_mageck_Proliferation, 
                           Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")),
             aes(x = RKO_Proliferation_Zscore, y = RKO_pscore, color = RKO),size = 2)+
  geom_text_repel(data = subset(Brunello_mageck_Proliferation, 
                                Symbol %in% c("CTNNB1","MED12","LEF1","TCF7L2","BCL9L","AXIN2","APC","CSNK1A1")), 
                  aes(x = RKO_Proliferation_Zscore, y = RKO_pscore, color = RKO,label = Symbol),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_y = 0.1,
                  segment.color = NA,
                  segment.size = 0)+
  theme_classic()+
  ylim(0,6)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure S2A.pdf")


dev.off()


