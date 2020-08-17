#author Chunhua Wan (Chunhua.wan@monash.edu)

require("openxlsx")
require(tidyverse)
require(ggrepel)

##load data
#Brunello sgRNA sequencing raw reads
brunello_rawread <- read.xlsx("./brunello_rawreads.xlsx")

#Brunello MAGEcK-VISPR datasets, lowGFP, highGFP and proliferation
Enrichment_list <- read.xlsx("./Brunello_lowGFP.xlsx")

Enrichment_highGFP_list <- read.xlsx("./Brunello_HiGFP.xlsx")

brunello_proliferation <- read.xlsx("./Proliferation_D21_plasmid.xlsx")

##depmap data downloaded on 24-July-2019
#Common essential gene list
Essential_gene <- read.csv("./Achilles_common_essentials.csv", header = T, stringsAsFactors = F)
#CCLE expression
DLD1_RKO <- read.xlsx("./DLD1_RKO_expression.xlsx")

# for Table S1
df1 <- Enrichment_list %>% 
  left_join(Enrichment_highGFP_list, by = "Symbol") %>% 
  left_join(brunello_proliferation, by = "Symbol") %>% 
  left_join(DLD1_RKO, by = "Symbol") %>% 
  left_join(Essential_gene, by = "Symbol")

write.csv(df1, "Table S1.csv")

##process raw data for volcanoplot
#normalize the reads and calculate the fold-change values per sgRNA
brunello <- brunello_rawread %>% 
  transmute(Symbol = Symbol,
            Brunelo_plasmid = Brunello_Plasmid/sum(.$Brunello_Plasmid)*1000000,
            DLD1_axin_1st_DAY21 = DLD1_axin_1st_DAY21_S7/sum(.$DLD1_axin_1st_DAY21_S7)*1000000,
            DLD1_axin_2nd_DAY21 = DLD1_axin_2nd_DAY21_S8/sum(.$DLD1_axin_2nd_DAY21_S8)*1000000,
            DLD1_axin_1st_DAY7 = DLD1_axin_1st_DAY7_S5/sum(.$DLD1_axin_1st_DAY7_S5)*1000000,
            DLD1_axin_2nd_DAY7 = DLD1_axin_2nd_DAY7_S6/sum(.$DLD1_axin_2nd_DAY7_S6)*1000000,
            DLD1_axin_1st_LowGFP = DLD1_axin_1st_LowGFP_S9/sum(.$DLD1_axin_1st_LowGFP_S9)*1000000,
            DLD1_axin_2nd_LowGFP = DLD1_axin_2nd_LowGFP_S10/sum(.$DLD1_axin_2nd_LowGFP_S10)*1000000,
            TOPGC_1st_DAY7 = TOPGC_1st_DAY7_S1/sum(.$TOPGC_1st_DAY7_S1)*1000000,
            TOPGC_2nd_DAY7 = TOPGC_2nd_DAY7_S2/sum(.$TOPGC_2nd_DAY7_S2)*1000000,
            TOPGC_1st_LowGFP = TOPGC_1st_LowGFP_S3/sum(.$TOPGC_1st_LowGFP_S3)*1000000,
            TOPGC_2nd_LowGFP = TOPGC_2nd_LowGFP_S4/sum(.$TOPGC_2nd_LowGFP_S4)*1000000,
            DLD1_cmyc1st_Day7 = DLD1_cmyc1st_Day7_S1/sum(.$DLD1_cmyc1st_Day7_S1)*1000000,
            DLD1_cmyc2nd_Day7 = DLD1_cmyc2nd_Day7_S2/sum(.$DLD1_cmyc2nd_Day7_S2)*1000000,
            DLD1_cmyc1st_Day21 = DLD1_cmyc1st_Day21_S3/sum(.$DLD1_cmyc1st_Day21_S3)*1000000,
            DLD1_cmyc2nd_Day21 = DLD1_cmyc2nd_Day21_S4/sum(.$DLD1_cmyc2nd_Day21_S4)*1000000,
            DLD1_cmyc1st_LowGFP = DLD1_cmyc1st_LowGFP_S5/sum(.$DLD1_cmyc1st_LowGFP_S5)*1000000,
            DLD1_cmyc2nd_LowGFP = DLD1_cmyc2nd_LowGFP_S6/sum(.$DLD1_cmyc2nd_LowGFP_S6)*1000000,
            RKO_1st_Day7 = RKO_1st_Day7_S7/sum(.$RKO_1st_Day7_S7)*1000000,
            RKO_2nd_Day7 =  RKO_2nd_Day7_S8/sum(.$RKO_2nd_Day7_S8)*1000000,
            RKO_1st_Day21 = RKO_1st_Day21_S9/sum(.$RKO_1st_Day21_S9)*1000000,
            RKO_2nd_Day21 = RKO_2nd_Day21_S10/sum(.$RKO_2nd_Day21_S10)*1000000,
            RKO_1st_LowGFP =  RKO_1st_LowGFP_S11/sum(.$RKO_1st_LowGFP_S11)*1000000,
            RKO_2nd_LowGFP =  RKO_2nd_LowGFP_S12/sum(.$RKO_2nd_LowGFP_S12)*1000000,
            TOP_GC_1st_HiGFP = TOP_GC_1st_HiGFP/sum(.$TOP_GC_1st_HiGFP)*1000000,
            TOP_GC_2nd_HiGFP = TOP_GC_2nd_HiGFP/sum(.$TOP_GC_2nd_HiGFP)*1000000,
            DLD_myc_1st_HiGFP = DLD_myc_1st_HiGFP/sum(.$DLD_myc_1st_HiGFP)*1000000,
            DLD_myc_2nd_HiGFP  =  DLD_myc_2nd_HiGFP/sum(.$DLD_myc_2nd_HiGFP)*1000000,
            DLD1_axin_1st_HiGFP = DLD1_axin_1st_HiGFP/sum(.$DLD1_axin_1st_HiGFP)*1000000,
            DLD1_axin_2nd_HiGFP = DLD1_axin_2nd_HiGFP/sum(.$DLD1_axin_2nd_HiGFP)*1000000,
            RKO_myc_1st_HiGFP =  RKO_myc_1st_HiGFP/sum(.$RKO_myc_1st_HiGFP)*1000000,
            RKO_myc_2nd_HiGFP =  RKO_myc_2nd_HiGFP/sum(.$RKO_myc_2nd_HiGFP)*1000000
            ) %>% 
  transmute(Symbol = Symbol,
            DLD1_Axin2 = (DLD1_axin_1st_LowGFP + DLD1_axin_2nd_LowGFP)/(DLD1_axin_1st_DAY7 + DLD1_axin_2nd_DAY7),
            DLD1_TOPGC = (TOPGC_1st_LowGFP + TOPGC_2nd_LowGFP)/(TOPGC_1st_DAY7 + TOPGC_2nd_DAY7),
            DLD1_MYC = (DLD1_cmyc1st_LowGFP + DLD1_cmyc2nd_LowGFP)/(DLD1_cmyc1st_Day7 + DLD1_cmyc2nd_Day7),
            RKO_MYC = (RKO_1st_LowGFP + RKO_2nd_LowGFP)/(RKO_1st_Day7 + RKO_2nd_Day7),
            DLD1_Axin2_HighGFP = (DLD1_axin_1st_HiGFP + DLD1_axin_2nd_HiGFP)/(DLD1_axin_1st_DAY7 + DLD1_axin_2nd_DAY7),
            DLD1_TOPGC_HighGFP = (TOP_GC_1st_HiGFP + TOP_GC_2nd_HiGFP)/(TOPGC_1st_DAY7 + TOPGC_2nd_DAY7),
            DLD1_MYC_HighGFP = (DLD_myc_1st_HiGFP + DLD_myc_2nd_HiGFP)/(DLD1_cmyc1st_Day7 + DLD1_cmyc2nd_Day7),
            RKO_MYC_HighGFP = (RKO_myc_1st_HiGFP + RKO_myc_2nd_HiGFP)/(RKO_1st_Day7 + RKO_2nd_Day7),
            DLD1_proliferation_FC = (DLD1_cmyc1st_Day21 + DLD1_cmyc2nd_Day21)/Brunelo_plasmid * 0.5,
            RKO_proliferaiton_FC = (RKO_1st_Day21 + RKO_2nd_Day21)/Brunelo_plasmid * 0.5
            ) %>% 
  drop_na()

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
  inner_join(Enrichment_list, by = "Symbol") %>% 
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
  mutate(TOPGC = ifelse(DLD1_TOPGC_low_score >= 0, "N", "T"),
         AXIN2 = ifelse(DLD1_AXIN2_low_score >= 0, "N", "T"),
         MYC = ifelse(DLD1_MYC_low_score >= 0, "N", "T"),
         RKO = ifelse(RKO_MYC_low_score >= 0, "N", "T"))  

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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave(filename = "Figure 1D.pdf", width = 4.05, height = 3.6)

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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 2C.pdf", width = 4.05, height = 3.6)

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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 2D.pdf", width = 4.05, height = 3.6)

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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure 2E.pdf", width = 4.05, height = 3.6)

###highGFP plot

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
  inner_join(brunello1, by = "Symbol") %>% 
  mutate(TOPGC = ifelse(TOPGC_highZscore >= 0, "T","N"),
         AXIN2 = ifelse(AXIN2_highZscore >= 0, "T", "N"),
         MYC = ifelse(MYC_highZscore >= 0, "T", "N"),
         RKO = ifelse(RKO_highZscore >= 0, "T", "N"))  

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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure S1B.pdf", width = 4.05, height = 3.6)


##make proliferation plots
str(brunello_proliferation)

Brunello_mageck_Proliferation <- brunello_proliferation %>% 
  transmute(Symbol= Symbol,
         DLD1_pscore = ifelse(is.finite(-log10(DLD1_Proliferation_pvalue)),
                               -log10(DLD1_Proliferation_pvalue),5),
         RKO_pscore = ifelse(is.finite(-log10(RKO_Proliferation_pvalue)),
                               -log10(RKO_Proliferation_pvalue),5)
  ) %>% 
  inner_join(brunello1, by = "Symbol") %>% 
  mutate(DLD1 = ifelse(DLD1_Proliferation_Zscore >= 0, "T", "N"),
         RKO = ifelse(RKO_Proliferation_Zscore >= 0, "T", "N")) 

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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
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
  theme(legend.position = "none", axis.text=element_text(size=12))+
  ylim(0,5.5)+
  xlab("Z-score")+ylab("-log10(P-value)")

ggsave("Figure S2A.pdf")


dev.off()
