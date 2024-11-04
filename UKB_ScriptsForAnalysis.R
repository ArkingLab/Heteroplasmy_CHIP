library(dplyr)
library(ggplot2)
library(foreign)
library(haven)
library(e1071)
library(stringr)
library(survival)
library(survminer)
library(ggrepel)
library(gghalves)
rm(list=ls())

GeneralTheme <- theme(axis.title = element_text(color = "black", size = 20),
                      axis.text = element_text(color = "black", size = 18),
                      strip.text = element_text(size = 18, color = "black"),
                      axis.ticks = element_line(color = "black"),
                      legend.title = element_text(color = "black", size = 20),
                      legend.text = element_text(color = "black", size = 18),
                      legend.position = "top")



# Figure1 UKB -------------------------------------------------------------

UKB_Freq_Tab <- fin_df%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  arrange(-Freq)
  
#Supplementary Figure 1A 
UKB_Freq_Tab%>%
  head(n=15)%>%
  ggplot(aes(x = reorder(Var1, Freq), y = Freq))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  xlab("Gene")+
  ylab("Number of individuals")+
  coord_flip()


#Figure 1A
UKB_Freq_Tab%>%
  filter(Var1 %in% c("DNMT3A", "TET2", "ASXL1", "PPM1D","TP53", "SRSF2", "SF3B1", "U2AF1"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("DNMT3A", "TET2", "ASXL1", "PPM1D","TP53", "SRSF2", "SF3B1", "U2AF1"))))%>%
  ggplot(aes(x = Var1, y = Freq))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  xlab("Gene")+
  ylab("Number of individuals")+
  coord_flip()



FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D","TP53", "SRSF2", "SF3B1", "U2AF1")


#Figure 1C
fin_df%>%
  filter(Gene.refGene %in% FreqGenesToAssess)%>%
  mutate(Gene.refGene = factor(Gene.refGene, levels = rev(FreqGenesToAssess)))%>%
  ggplot(aes(x = Gene.refGene, y = VAF))+
  geom_violin(scale = "width")+
  geom_boxplot(alpha = 0, width = 0.15)+
  theme_classic()+
  GeneralTheme+
  xlab("Gene")+
  ylab("VAF")+
  #coord_cartesian(ylim = c(0, 100))+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = paste0(c(0, 25, 50, 75, 100), "%"),
                     limits = c(0, 100))+
  coord_flip()
  
CHIP_int_df <- fin_df
CHIP_int_df$SpliceHit <- CHIP_int_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")

wilcox.test(CHIP_int_df$VAF ~ CHIP_int_df$SpliceHit)

#Figure 1E
table(fin_df$SampID)%>%
  data.frame()%>%
  mutate(FreqN = ifelse(Freq>3, 4, Freq))%>%
  pull(FreqN)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  mutate(Var1 = as.numeric(as.character(Var1)))%>%
  ggplot(aes(x = Freq, y = reorder(Var1, -Var1)))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Number of mutations")+
  xlab("Number of individuals")




# Figure2 UKB -------------------------------------------------------------

#Figure2A
HCCOMP_DF <- (PhenoDF[,c(1, which(grepl("^hetcount_complex_", colnames(PhenoDF))))])%>%
  melt(., id.vars = "id")

table(HCCOMP_DF$variable, HCCOMP_DF$value>0)%>%
  data.frame()%>%
  filter(Var2 != FALSE)%>%
  mutate(Var1 = str_remove(Var1, "hetcount_complex_"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                        labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  ggplot(aes(x = Freq, y = Var1))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Complex")+
  xlab("Number of individuals")


#Figure2C
KIV <- HCCOMP_DF%>%
  filter(value>0)%>%
  mutate(CombVar = paste0(id, "_", variable))%>%
  pull(CombVar)%>%
  as.character()%>%
  str_replace(., "hetcount", "mMSS")

MSSCOMP_DF <- (PhenoDF[,c(1, which(grepl("^mMSS_complex_", colnames(PhenoDF))))])%>%
  melt(., id.vars = "id")

MSSCOMP_DF <- MSSCOMP_DF[paste0(MSSCOMP_DF$id, "_", MSSCOMP_DF$variable) %in% KIV,]




MSSCOMP_DF%>%
  mutate(variable = str_remove(variable, "mMSS_complex_"))%>%
  mutate(variable = factor(variable, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                           labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  ggplot(aes(y = variable, x = value))+
  geom_violin(scale = "width")+
  geom_boxplot(alpha = 0, width = 0.15)+
  theme_classic()+
  GeneralTheme+
  xlab("mMSS")+
  ylab("Complex")+
  coord_cartesian(xlim = c(0, 2))+
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2))


MSSCOMP_DF$TRvar <- MSSCOMP_DF$variable %in% c("mMSS_complex_RRNA", "mMSS_complex_TRNA")
wilcox.test(MSSCOMP_DF$value ~ MSSCOMP_DF$TRvar)


#Figure2E
table(PhenoDF$count_het)%>%
  data.frame()%>%
  mutate(Var1 = as.numeric(as.character(Var1)))%>%
  filter(Var1>0)%>%
  ggplot(aes(x = Freq, y = reorder(Var1, -Var1)))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Number of heteroplasmies")+
  xlab("Number of individuals")


# Figure 3A UKB -----------------------------------------------------------

prop.table(table(PhenoDF$AnyCHIP_init_filt, PhenoDF$het_yn))*100

# Figure 3C UKB-------------------------------------------------------------

FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D","TP53", "SRSF2", "SF3B1", "U2AF1")

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()


for(i in c("AnyCHIP_init_filt", "AnyCHIP_init_filt_Small", "AnyCHIP_init_filt_VAF20", "NOMsingle", "NOMmulti")){
  
  int_df <- PhenoDF
  
  int_df$HCov0 <- int_df$count_het>0

  
  int_CHIP <- fin_df
  
  int_CHIP$Gene <- int_CHIP$Gene.refGene
  

  
  if(grepl("VAF20", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_VAF20 == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyCHIP_init_filt_Small"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_Small == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_CHIP > 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]

        ForLOGDF <- int_df%>%filter(count_het==0 | het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
          #filter(!is.na(smk_ever))
        
        si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
          data.frame()%>%
          filter(IsGOI == TRUE)%>%
          mutate(s = sum(Freq))%>%
          filter(het_VOI == TRUE)%>%
          mutate(P = Freq/s*100)
        
        
        
        Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
        si_logistic_sum <- summary(Logistic_multi)
        
        CHIPfiltch <- c(CHIPfiltch, i)
        Genech <- c(Genech, g)
        hvarch <- c(hvarch, hvar)
        
        if(length(si_ct$IsGOI)==1){
          
          HetCHch <- c(HetCHch, si_ct$Freq)
          TotCHch <- c(TotCHch, si_ct$s)
          PercHet <- c(PercHet, si_ct$P)
          
        }else{
          
          HetCHch <- c(HetCHch, NA)
          TotCHch <- c(TotCHch, NA)
          PercHet <- c(PercHet, NA)
          
        }
        

        pch <- c(pch, (si_logistic_sum$coefficients[2,4]))
        
        rm(Logistic_multi)
        rm(si_logistic_sum)
        

    }
  }
}




si_ct <- with(PhenoDF, table((AnyCHIP_init_filt==FALSE), (count_het>0)))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  mutate(s = sum(Freq))%>%
  filter(Var2 == TRUE)%>%
  mutate(P = Freq/s*100)



CHIPfiltch <- c(CHIPfiltch, "AnyCHIP_init_filt")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
pch <- c(pch, NA)

CHIPfiltch <- c(CHIPfiltch, "AnyCHIP_init_filt_Small")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
pch <- c(pch, NA)

CHIPfiltch <- c(CHIPfiltch, "AnyCHIP_init_filt_VAF20")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
pch <- c(pch, NA)


FreqGeneCHIPhet_df <- data.frame(CHIPfiltch, Genech, hvarch, HetCHch, TotCHch, PercHet, pch)

FreqGeneCHIPhet_df$P_adj <- p.adjust(FreqGeneCHIPhet_df$pch, method = "BH")


FreqGeneCHIPhet_df <- FreqGeneCHIPhet_df%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyCHIP_init_filt_Small No CHIP", "AnyCHIP_init_filt any", "AnyCHIP_init_filt_Small any", "AnyCHIP_init_filt_VAF20 any", "NOMsingle any", "NOMmulti any", paste0("AnyCHIP_init_filt ", FreqGenesToAssess))),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<20%", "VAF>20%", "Single mutation", "Multiple mutations", FreqGenesToAssess))))%>%
  filter(!(is.na(PlotGF)))

FreqGeneCHIPhet_df$SplVar <- case_when(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("2%<VAF<20%", "VAF>20%") ~ "xSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                                       TRUE ~ "zG")


FreqGeneCHIPhet_df$Significant <- ifelse(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP"), NA, "*")
FreqGeneCHIPhet_df$Significant <- ifelse(is.na(FreqGeneCHIPhet_df$Significant), "Not significant", "Significant")


FreqGeneCHIPhet_df%>%
  ggplot(aes(y = PlotGF, x = PercHet, fill = Significant))+
  geom_col(position = position_dodge())+
  theme_classic()+
  GeneralTheme+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  geom_text(aes(x = 1, y = PlotGF, label = HetCHch), color = "white", hjust = 0, size = 6)+
  ylab("Gene")+
  xlab("Heteroplasmy prevalence")



ForLOGDF <- PhenoDF
Logistic_multi <- glm(het_yn ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt == TRUE)

Logistic_multi <- glm(het_yn ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

Logistic_multi <- glm(het_yn ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)


ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt == TRUE)

ForLOGDF$IsGOI <- ForLOGDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])

Logistic_multi <- glm(het_yn ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)


# Figure 3E UKB-------------------------------------------------------------

FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D","TP53", "SRSF2", "SF3B1", "U2AF1")

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()


for(i in c("AnyCHIP_init_filt", "AnyCHIP_init_filt_Small", "AnyCHIP_init_filt_VAF20", "NOMsingle", "NOMmulti")){
  
  int_df <- PhenoDF%>%
    filter(count_het>0)
  
  int_df$HCov0 <- int_df$count_het>1

  
  int_CHIP <- fin_df

  int_CHIP$Gene <- int_CHIP$Gene.refGene
  

  
  if(grepl("VAF20", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_VAF20 == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyCHIP_init_filt_Small"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_Small == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_CHIP > 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]

      ForLOGDF <- int_df%>%filter(count_het==1 | het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
      #filter(!is.na(smk_ever))
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
      si_logistic_sum <- summary(Logistic_multi)
      
      CHIPfiltch <- c(CHIPfiltch, i)
      Genech <- c(Genech, g)
      hvarch <- c(hvarch, hvar)
      
      if(length(si_ct$IsGOI)==1){
        
        HetCHch <- c(HetCHch, si_ct$Freq)
        TotCHch <- c(TotCHch, si_ct$s)
        PercHet <- c(PercHet, si_ct$P)
        
      }else{
        
        HetCHch <- c(HetCHch, NA)
        TotCHch <- c(TotCHch, NA)
        PercHet <- c(PercHet, NA)
        
      }
      

      pch <- c(pch, (si_logistic_sum$coefficients[2,4]))
      
      rm(Logistic_multi)
      rm(si_logistic_sum)
      

    }
  }
}


si_ct <- with(PhenoDF%>%filter(count_het>0), table((AnyCHIP_init_filt==FALSE), (count_het>1)))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  mutate(s = sum(Freq))%>%
  filter(Var2 == TRUE)%>%
  mutate(P = Freq/s*100)



CHIPfiltch <- c(CHIPfiltch, "AnyCHIP_init_filt")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
pch <- c(pch, NA)

CHIPfiltch <- c(CHIPfiltch, "AnyCHIP_init_filt_Small")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
pch <- c(pch, NA)

CHIPfiltch <- c(CHIPfiltch, "AnyCHIP_init_filt_VAF20")
Genech <- c(Genech, "No CHIP")
hvarch <- c(hvarch, "HCov0")
HetCHch <- c(HetCHch, si_ct$Freq)
TotCHch <- c(TotCHch, si_ct$s)
PercHet <- c(PercHet, si_ct$P)
pch <- c(pch, NA)


FreqGeneCHIPhet_df <- data.frame(CHIPfiltch, Genech, hvarch, HetCHch, TotCHch, PercHet, pch)



FreqGeneCHIPhet_df <- FreqGeneCHIPhet_df%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyCHIP_init_filt_Small No CHIP", "AnyCHIP_init_filt any", "AnyCHIP_init_filt_Small any", "AnyCHIP_init_filt_VAF20 any", "NOMsingle any", "NOMmulti any", paste0("AnyCHIP_init_filt ", FreqGenesToAssess))),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<20%", "VAF>20%", "Single mutation", "Multiple mutations", FreqGenesToAssess))))%>%
  filter(!(is.na(PlotGF)))

FreqGeneCHIPhet_df$SplVar <- case_when(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("2%<VAF<20%", "VAF>20%") ~ "xSC",
                                       FreqGeneCHIPhet_df$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                                       TRUE ~ "zG")


FreqGeneCHIPhet_df$Significant <- ifelse(FreqGeneCHIPhet_df$PlotGF %in% c("No CHIP"), NA, "*")
FreqGeneCHIPhet_df$Significant <- ifelse(is.na(FreqGeneCHIPhet_df$Significant), "Not significant", "Significant")


FreqGeneCHIPhet_df%>%
  ggplot(aes(y = PlotGF, x = PercHet, fill = Significant))+
  geom_col(position = position_dodge())+
  theme_classic()+
  GeneralTheme+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  geom_text(aes(x = 1, y = PlotGF, label = HetCHch), color = "white", hjust = 0, size = 6)+
  ylab("Gene")+
  xlab("Prevalence of multiple heteroplasmies")




ForLOGDF <- PhenoDF%>%
  filter(count_het>0)
Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt == TRUE & count_het>0)

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt == TRUE & count_het>0)

ForLOGDF$IsGOI <- ForLOGDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])

Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)






# Figure 3G UKB------------------------------------------

FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D","TP53", "SRSF2", "SF3B1", "U2AF1")

UKB_LIN_DF <- data.frame()


for(i in c("AnyCHIP_init_filt", "AnyCHIP_init_filt_Small", "AnyCHIP_init_filt_VAF20", "NOMsingle", "NOMmulti")){
  
  int_df <- PhenoDF
  int_df$LogmMSS <- log10(int_df$mMSS+1)
  int_df$HCov0 <- int_df$count_het>0
  
  int_CHIP <- fin_df
  int_CHIP$Gene <- int_CHIP$Gene.refGene
  

  
  if(grepl("VAF20", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_VAF20 == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyCHIP_init_filt_Small"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_Small == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_CHIP > 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      
        
        ForLOGDF <- int_df%>%filter(het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
        
        
        ToBindDf <- ForLOGDF%>%
          filter(IsGOI)%>%
          select(mMSS)
        
        if(length(ToBindDf$mMSS)==0){
          next
          }
          
          Linear_multi <- lm(LogmMSS ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
          si_linear_sum <- summary(Linear_multi)

          ToBindDf$CHIPfiltch <- i
          ToBindDf$Genech <- g
          ToBindDf$hvarch <- hvar
          ToBindDf$pch <- (si_linear_sum$coefficients[2,4])
          
          UKB_LIN_DF <- rbind(UKB_LIN_DF, ToBindDf)

          rm(Linear_multi)
          rm(si_linear_sum)

    }
  }
}

UKB_LIN_DF <- PhenoDF%>%
  filter(AnyCHIP_init_filt==FALSE & count_het>0)%>%
  select(mMSS)%>%
  mutate(CHIPfiltch = "AnyCHIP_init_filt", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(UKB_LIN_DF, .)

UKB_LIN_DF <- PhenoDF%>%
  filter(AnyCHIP_init_filt==FALSE & count_het>0)%>%
  select(mMSS)%>%
  mutate(CHIPfiltch = "AnyCHIP_init_filt_Small", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(UKB_LIN_DF, .)

UKB_LIN_DF <- PhenoDF%>%
  filter(AnyCHIP_init_filt==FALSE & count_het>0)%>%
  select(mMSS)%>%
  mutate(CHIPfiltch = "AnyCHIP_init_filt_VAF20", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(UKB_LIN_DF, .)

UKB_LIN_DF <- UKB_LIN_DF%>%
  distinct(CHIPfiltch, Genech, hvarch, pch)%>%
  mutate(P_adj = p.adjust(pch, method = "BH"))%>%
  select(-pch)%>%
  left_join(UKB_LIN_DF, ., by = c("CHIPfiltch", "Genech", "hvarch"))





UKB_LIN_DF <- UKB_LIN_DF%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyCHIP_init_filt_Small No CHIP", "AnyCHIP_init_filt any", "AnyCHIP_init_filt_Small any", "AnyCHIP_init_filt_VAF20 any", "NOMsingle any", "NOMmulti any", paste0("AnyCHIP_init_filt ", FreqGenesToAssess))),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<20%", "VAF>20%", "Single mutation", "Multiple mutations", FreqGenesToAssess))))%>%
  filter(!(is.na(PlotGF)))

UKB_LIN_DF$SplVar <- case_when(UKB_LIN_DF$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                                       UKB_LIN_DF$PlotGF %in% c("2%<VAF<20%", "VAF>20%") ~ "xSC",
                                       UKB_LIN_DF$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                                       TRUE ~ "zG")


UKB_LIN_DF$Significant <- ifelse(UKB_LIN_DF$PlotGF %in% c("No CHIP"), NA, "*")
UKB_LIN_DF$Significant <- ifelse(is.na(UKB_LIN_DF$Significant), "Not significant", "Significant")




UKB_LIN_DF%>%
  ggplot(aes(x = PlotGF, y = log10(mMSS+1), fill = Significant))+
  geom_half_violin(side = "r", scale = "width")+
  geom_half_boxplot(side = "l", alpha = 0)+
  theme_classic()+
  GeneralTheme+
  facet_grid(SplVar~., scales = "free_y", space = "free_y")+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  ylab("Log10(mMSS+1)")+
  xlab("Gene")+
  coord_flip(ylim = c(0, log10(2+1)))+
  scale_y_continuous(breaks = c(log10(c(0, 0.5, 1, 1.5, 2)+1)),
                     labels = c(0, 0.5, 1, 1.5, 2))




by(UKB_LIN_DF$mMSS, paste0(UKB_LIN_DF$CHIPfiltch, UKB_LIN_DF$Genech), quantile)


ForLOGDF <- PhenoDF%>%filter(count_het>0)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
Linear_multi <- lm(LogmMSS ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4)+count_het + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$AnyCHIP_init_filt, quantile)

ForLOGDF <- PhenoDF%>%filter(count_het>0)%>%filter(AnyCHIP_init_filt == TRUE)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
Linear_multi <- lm(LogmMSS ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4)+count_het + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$AnyCHIP_init_filt_VAF20, quantile)


Linear_multi <- lm(LogmMSS ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) +count_het+ sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$NOM_CHIP>1, quantile)


ForLOGDF <- PhenoDF%>%filter(count_het>0)%>%filter(AnyCHIP_init_filt == TRUE)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
ForLOGDF$IsGOI <- ForLOGDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])
Linear_multi <- lm(LogmMSS ~ as.numeric(IsGOI) + rms::rcs(age, df = 4)+count_het + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)
by(ForLOGDF$mMSS, ForLOGDF$IsGOI, quantile)




# Figure 4A,C UKB -----------------------------------------------------------

int_df <- PhenoDF
m <- "AnyCHIP_init_filt"

int_df$BindingHetVAR <- int_df$count_het>0

int_df <- int_df[((int_df$AnyCHIP_init_filt==0)|(int_df[,which(colnames(int_df) == m)])),]
int_df <- int_df[((int_df$BindingHetVAR)|(int_df$count_het==0)),]

int_df$AnyCHIP_VAR <- (int_df[,which(colnames(int_df) == m)])
int_df$CHIP_VAR <- paste0((int_df[,which(colnames(int_df) == m)]),
                          (int_df$BindingHetVAR))

int_df$CHIP_VAR <- factor(int_df$CHIP_VAR, levels = c("FALSEFALSE", "FALSETRUE", "TRUEFALSE", "TRUETRUE"),
                          labels = c("NO_CHIP_NO_Het", "NO_CHIP_YES_Het", "YES_CHIP_NO_Het", "YES_CHIP_YES_Het"))

surv_object <- Surv(time = int_df$first_incidence_fin, event = int_df$HEME_CENSOR)

fit.surv <- survfit(surv_object ~ CHIP_VAR, data = int_df)
fit.surv <- survfit(surv_object ~ BindingHetVAR, data = int_df)

ggsurvplot(fit.surv, data = int_df, 
           censor = FALSE,
           fun = "event",
           palette= c("#595959FF", "#00A087FF", "#4DBBD5FF", "#E64B35FF"),
           risk.table.col="strata",
           risk.table.y.text=FALSE,
           break.time.by=5,
           surv.scale="percent",
           ylab="% Individuals",
           xlab="Years",
           risk.table = TRUE,
           ggtheme=theme(axis.text = element_text(color = "black", size = 18),
                         axis.title = element_text(color = "black", size = 20),
                         axis.line = element_line(color = "black"),
                         axis.ticks = element_blank(),
                         panel.background = element_blank(),
                         panel.grid = element_blank())
)



fit.coxph <- coxph(surv_object ~ as.numeric(BindingHetVAR) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, 
                   data = int_df)

fit.coxph <- coxph(surv_object ~ CHIP_VAR + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, 
                   data = int_df)


exp(fit.coxph[[1]])
exp(confint(fit.coxph))
summary(fit.coxph)


# Table 1 UKB -----------------------------------------------------------------

t.test(PhenoDF$age ~ PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$sex, PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$race_na=="white", PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$smk_ever, PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$anemia, PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$thrombocytopenia, PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$neutropenia, PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$cytopenia_yn, PhenoDF$AnyCHIP_init_filt)

t.test(PhenoDF$mcv ~ PhenoDF$AnyCHIP_init_filt)

t.test(PhenoDF$rbcdw ~ PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$prev_cancer_yn, PhenoDF$AnyCHIP_init_filt)

fisher.test(PhenoDF$het_yn, PhenoDF$AnyCHIP_init_filt)

# Table 2 UKB -----------------------------------------------------------------

t.test(PhenoDF$age ~ PhenoDF$het_yn)

fisher.test(PhenoDF$sex, PhenoDF$het_yn)

fisher.test(PhenoDF$race_na=="white", PhenoDF$het_yn)

fisher.test(PhenoDF$smk_ever, PhenoDF$het_yn)

fisher.test(PhenoDF$anemia, PhenoDF$het_yn)

fisher.test(PhenoDF$thrombocytopenia, PhenoDF$het_yn)

fisher.test(PhenoDF$neutropenia, PhenoDF$het_yn)

fisher.test(PhenoDF$cytopenia_yn, PhenoDF$het_yn)

t.test(PhenoDF$mcv ~ PhenoDF$het_yn)

t.test(PhenoDF$rbcdw ~ PhenoDF$het_yn)

fisher.test(PhenoDF$prev_cancer_yn, PhenoDF$het_yn)

fisher.test(PhenoDF$AnyCHIP_init_filt, PhenoDF$het_yn)

# Supplementary Figure 2A UKB ----------------------------------------------

FreqGenesToAssess <- fin_df%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  unique()%>%
  as.character()

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()

Medch <- c()


for(i in c("AnyCHIP_init_filt")){
  
  int_df <- PhenoDF
  
  int_df$HCov0 <- int_df$count_het>0
  
  int_CHIP <- fin_df
  int_CHIP$Gene <- int_CHIP$Gene.refGene
  

  
  if(grepl("VAF20", i)){
    int_CHIP <- int_CHIP[int_CHIP$VAF>=20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_VAF20 == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  
  if(i == "AnyCHIP_init_filt_Small"){
    int_CHIP <- int_CHIP[int_CHIP$VAF<20,]
    int_df <- int_df[(int_df$AnyCHIP_init_filt_Small == TRUE)|(int_df$AnyCHIP_init_filt == FALSE),]
  }
  
  if(i == "NOMmulti"){
    int_df <- int_df[(int_df$NOM_CHIP > 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$SampID %in% int_df$SampID,]
    
  }
  
  
  
  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$SampID %in% int_CHIP$SampID
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
    
      ForLOGDF <- int_df%>%filter(count_het==0 | het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
      
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
      si_logistic_sum <- summary(Logistic_multi)
      
      CHIPfiltch <- c(CHIPfiltch, i)
      Genech <- c(Genech, g)
      hvarch <- c(hvarch, hvar)
      
      if(length(si_ct$IsGOI)==1){
        
        HetCHch <- c(HetCHch, si_ct$Freq)
        TotCHch <- c(TotCHch, si_ct$s)
        PercHet <- c(PercHet, si_ct$P)
        
        if(si_ct$Freq == 0){
          Medch <- c(Medch, 0)
        }else{
          Medch <- c(Medch, (ForLOGDF%>%
                               filter(IsGOI & het_VOI)%>%
                               pull(mMSS)%>%
                               median()))
          
        }
        
      }else{
        
        HetCHch <- c(HetCHch, NA)
        TotCHch <- c(TotCHch, NA)
        PercHet <- c(PercHet, NA)
        
      }
      

      pch <- c(pch, (si_logistic_sum$coefficients[2,4]))
      
      rm(Logistic_multi)
      rm(si_logistic_sum)
      
    }
  }
}

FreqGeneCHIPhet_df <- data.frame(CHIPfiltch, Genech, hvarch, HetCHch, TotCHch, PercHet, pch, Medch)

FreqGeneCHIPhet_df%>%
  filter(TotCHch > 10)%>%
  pull(Genech)%>%
  write.csv(., file = "GOI_UKB_OutAlt3.csv", row.names = FALSE, quote = FALSE)

FreqGeneCHIPhet_df$TotCHch <- ifelse(FreqGeneCHIPhet_df$Genech == "any", 1, FreqGeneCHIPhet_df$TotCHch)#To keep figure proportions; edited to the correct size in the final figure
FreqGeneCHIPhet_df$SzCat <- cut(FreqGeneCHIPhet_df$TotCHch, c(0, 10, 100, 1000, 100000))

options(ggrepel.max.overlaps = Inf)
FreqGeneCHIPhet_df%>%
  filter(TotCHch > 10|Genech == "any")%>%
  ggplot(aes(x = PercHet, y = Medch))+
  geom_point(aes(size = SzCat), color = "#595959FF")+
  geom_text_repel(aes(label = Genech), size = 2)+
  theme_classic()+
  GeneralTheme+
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     labels = paste0(seq(0, 100, by = 10), "%"),
                     limits = c(15, 80))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 0.8))+
  xlab("Heteroplasmy prevalence")+
  ylab("Median mMSS")
  



# Supplementary Figure 3 --------------------------------------------------

#Supplementary Figure 3A
with(PhenoDF, table(ccus, het_yn))%>%
  data.frame()%>%
  group_by(ccus)%>%
  mutate(s = sum(Freq))%>%
  ungroup()%>%
  mutate(P = Freq/s*100)%>%
  filter(het_yn == 1)%>%
  mutate(ccus = factor(ccus, levels = rev(c("none", "chip", "ccus")),
                       labels = rev(c("No CHIP", "CHIP\nwithout\ncytopenia", "CCUS"))))%>%
  mutate(Significant = ccus != "No CHIP")%>%
  ggplot(aes(y = ccus, x = P, fill = Significant))+
  geom_col(position = position_dodge())+
  geom_text(aes(x = 1, y = ccus, label = Freq), color = "white", hjust = 0, size = 6)+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  xlab("Heteroplasmy prevalence")+
  ylab("")

ForLOGDF <- PhenoDF[as.character(PhenoDF$ccus)!="ccus",]
ForLOGDF$IsGOI <- as.numeric(as.character(ForLOGDF$ccus) == "chip")
Logistic_multi <- glm(het_yn ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt == TRUE)
ForLOGDF$IsGOI <- as.numeric(as.character(ForLOGDF$ccus) == "ccus")
Logistic_multi <- glm(het_yn ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)


#Supplementary Figure3C
with(PhenoDF%>%filter(het_yn==1), table(ccus, count_het>1))%>%
  data.frame()%>%
  group_by(ccus)%>%
  mutate(s = sum(Freq))%>%
  ungroup()%>%
  mutate(P = Freq/s*100)%>%
  filter(Var2 == TRUE)%>%
  mutate(ccus = factor(ccus, levels = rev(c("none", "chip", "ccus")),
                       labels = rev(c("No CHIP", "CHIP\nwithout\ncytopenia", "CCUS"))))%>%
  mutate(Significant = ccus != "No CHIP")%>%
  ggplot(aes(y = ccus, x = P, fill = Significant))+
  geom_col(position = position_dodge())+
  geom_text(aes(x = 1, y = ccus, label = Freq), color = "white", hjust = 0, size = 6)+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  xlab("Prevalence of multiple heteroplasmies")+
  ylab("")

ForLOGDF <- (PhenoDF[as.character(PhenoDF$ccus)!="ccus",])%>%
  filter(count_het>0)
ForLOGDF$IsGOI <- as.numeric(as.character(ForLOGDF$ccus) == "chip")
Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)


ForLOGDF <- PhenoDF%>%
  filter(AnyCHIP_init_filt == TRUE & count_het>0)
ForLOGDF$IsGOI <- as.numeric(as.character(ForLOGDF$ccus) == "ccus")
Logistic_multi <- glm(as.numeric(count_het>1) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



PhenoDF%>%
  filter(het_yn==1)%>%
  mutate(ccus = factor(ccus, levels = rev(c("none", "chip", "ccus")),
                       labels = rev(c("No CHIP", "CHIP\nwithout\ncytopenia", "CCUS"))))%>%
  mutate(Significant = ccus != "No CHIP")%>%
  ggplot(aes(x = ccus, y = log10(mMSS+1), fill = Significant))+
  geom_half_violin(side = "r", scale = "width")+
  geom_half_boxplot(side = "l", alpha = 0)+
  #geom_half_point(side = "l", alpha=0.5)+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  theme(legend.position = "none")+
  ylab("Log10(mMSS+1)")+
  xlab("")+
  coord_flip(ylim = c(0, log10(2+1)))+
  scale_y_continuous(breaks = c(log10(c(0, 0.5, 1, 1.5, 2)+1)),
                     labels = c(0, 0.5, 1, 1.5, 2))


ForLOGDF <- (PhenoDF[as.character(PhenoDF$ccus)!="ccus",])%>%
  filter(count_het>0)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
ForLOGDF$IsGOI <- as.numeric(as.character(ForLOGDF$ccus) == "chip")
Linear_multi <- lm(LogmMSS ~ as.numeric(IsGOI) + rms::rcs(age, df = 4)+count_het + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)


ForLOGDF <- (PhenoDF[as.character(PhenoDF$ccus)!="none",])%>%
  filter(count_het>0)
ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
ForLOGDF$IsGOI <- as.numeric(as.character(ForLOGDF$ccus) == "ccus")
Linear_multi <- lm(LogmMSS ~ as.numeric(IsGOI) + rms::rcs(age, df = 4)+count_het + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
summary(Linear_multi)

# Supplementary Figure 4 ----------------------------------------------

HCCOMP_DF <- (PhenoDF[,c(1, which(grepl("^hetcount_complex_", colnames(PhenoDF))))])%>%
  melt(., id.vars = "id")

TBDF <- PhenoDF%>%
  select(SampID, AnyCHIP_init_filt, AnyCHIP_init_filt_VAF20, NOM_CHIP)%>%
  mutate(MultiHit = NOM_CHIP>1)%>%
  select(-NOM_CHIP)

TBDF$IsGOI <- TBDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])


HCCOMP_DF_TP <- left_join(HCCOMP_DF, TBDF, by = c("id" = "SampID"))

HCCOMP_DF_TP <- HCCOMP_DF_TP[HCCOMP_DF_TP$value>0,]

int_df <- HCCOMP_DF_TP

TabIntDF <- table(int_df$variable, int_df$AnyCHIP_init_filt)%>%
  data.frame()


TabIntDF <- table(PhenoDF$count_het>0, PhenoDF$AnyCHIP_init_filt)%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  select(-Var1)%>%
  left_join(TabIntDF, ., by = "Var2")%>%
  mutate(P = Freq.x/Freq.y*100)


#Supplementary Figure4A
TabIntDF%>%
  mutate(Var1 = str_remove(Var1, "hetcount_complex_"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                       labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(Var2 = factor(Var2, levels = rev(c(FALSE, TRUE)),
                       labels = rev(c("No CHIP", "CHIP"))))%>%
  ggplot(aes(x = P, y = Var1, fill = Var2))+
  geom_col(position = position_dodge(width = 0.9))+
  geom_text(aes(x = 1, y = Var1, label = Freq.x), color = "white", hjust = 0, size = 6, position = position_dodge(width = 0.9))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     labels = paste0(seq(0, 100, by = 10), "%"))+
  ylab("Complex")+
  xlab("Percent of all individuals with heteroplasmy")



ForLOGDF <- PhenoDF%>%
  filter(count_het>0)


nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
ForLOGDF$IsGOI <- as.numeric((ForLOGDF[,i])>0)

Logistic_multi <- glm(as.numeric(IsGOI) ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
si_tmp <- summary(Logistic_multi)

nch <- c(nch, colnames(ForLOGDF)[i])
pch <- c(pch, si_tmp$coefficients[,4][[2]])

print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)







int_df <- HCCOMP_DF_TP%>%
  filter(AnyCHIP_init_filt == TRUE)

TabIntDF <- table(int_df$variable, int_df$AnyCHIP_init_filt_VAF20)%>%
  data.frame()



TabIntDF <- with(PhenoDF%>%filter(AnyCHIP_init_filt == TRUE), table(count_het>0, AnyCHIP_init_filt_VAF20==TRUE))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  select(-Var1)%>%
  left_join(TabIntDF, ., by = "Var2")%>%
  mutate(P = Freq.x/Freq.y*100)

#Supplementary Figure4B
TabIntDF%>%
  mutate(Var1 = str_remove(Var1, "hetcount_complex_"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                       labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(Var2 = factor(Var2, levels = rev(c(FALSE, TRUE)),
                       labels = rev(c("2%<VAF<20%", "VAF>20%"))))%>%
  ggplot(aes(x = P, y = Var1, fill = Var2))+
  geom_col(position = position_dodge(width = 0.9))+
  geom_text(aes(x = 1, y = Var1, label = Freq.x), color = "white", hjust = 0, size = 6, position = position_dodge(width = 0.9))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     labels = paste0(seq(0, 100, by = 10), "%"))+
  ylab("Complex")+
  xlab("Percent of all individuals with heteroplasmy")


ForLOGDF <- PhenoDF%>%
  filter(count_het>0)%>%
  filter(AnyCHIP_init_filt == TRUE)


nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF$IsGOI <- as.numeric((ForLOGDF[,i])>0)
  
  Logistic_multi <- glm(as.numeric(IsGOI) ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
  si_tmp <- summary(Logistic_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)




int_df <- HCCOMP_DF_TP%>%
  filter(AnyCHIP_init_filt == TRUE)

TabIntDF <- table(int_df$variable, int_df$MultiHit)%>%
  data.frame()



TabIntDF <- with(PhenoDF%>%filter(AnyCHIP_init_filt == TRUE), table(count_het>0, NOM_CHIP>1))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  select(-Var1)%>%
  left_join(TabIntDF, ., by = "Var2")%>%
  mutate(P = Freq.x/Freq.y*100)


#Supplementary Figure4C
TabIntDF%>%
  mutate(Var1 = str_remove(Var1, "hetcount_complex_"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                       labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(Var2 = factor(Var2, levels = rev(c(FALSE, TRUE)),
                       labels = rev(c("Single mutation", "Multiple mutations"))))%>%
  ggplot(aes(x = P, y = Var1, fill = Var2))+
  geom_col(position = position_dodge(width = 0.9))+
  geom_text(aes(x = 1, y = Var1, label = Freq.x), color = "white", hjust = 0, size = 6, position = position_dodge(width = 0.9))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     labels = paste0(seq(0, 100, by = 10), "%"))+
  ylab("Complex")+
  xlab("Percent of all individuals with heteroplasmy")



ForLOGDF <- PhenoDF%>%
  filter(count_het>0)%>%
  filter(AnyCHIP_init_filt == TRUE)


nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF$IsGOI <- as.numeric((ForLOGDF[,i])>0)
  
  Logistic_multi <- glm(as.numeric(IsGOI) ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
  si_tmp <- summary(Logistic_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)






int_df <- HCCOMP_DF_TP%>%
  filter(AnyCHIP_init_filt == TRUE)

TabIntDF <- table(int_df$variable, int_df$IsGOI)%>%
  data.frame()



TabIntDF <- with(PhenoDF%>%filter(AnyCHIP_init_filt == TRUE), table(count_het>0, (SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")]))))%>%
  data.frame()%>%
  filter(Var1 == TRUE)%>%
  select(-Var1)%>%
  left_join(TabIntDF, ., by = "Var2")%>%
  mutate(P = Freq.x/Freq.y*100)


#Supplementary Figure4D
TabIntDF%>%
  mutate(Var1 = str_remove(Var1, "hetcount_complex_"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                       labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(Var2 = factor(Var2, levels = rev(c(FALSE, TRUE)),
                       labels = rev(c("Other CHIP", "Spliceosome"))))%>%
  ggplot(aes(x = P, y = Var1, fill = Var2))+
  geom_col(position = position_dodge(width = 0.9))+
  geom_text(aes(x = 1, y = Var1, label = Freq.x), color = "white", hjust = 0, size = 6, position = position_dodge(width = 0.9))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     labels = paste0(seq(0, 100, by = 10), "%"))+
  ylab("Complex")+
  xlab("Percent of all individuals with heteroplasmy")



ForLOGDF <- PhenoDF%>%
  filter(count_het>0)%>%
  filter(AnyCHIP_init_filt == TRUE)


ForLOGDF$IsGOISplice <- ForLOGDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])

nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF$IsGOI <- as.numeric((ForLOGDF[,i])>0)
  
  Logistic_multi <- glm(as.numeric(IsGOI) ~ as.numeric(IsGOISplice) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF, family = binomial)
  si_tmp <- summary(Logistic_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)







KIV <- HCCOMP_DF%>%
  filter(value>0)%>%
  mutate(CombVar = paste0(id, "_", variable))%>%
  pull(CombVar)%>%
  as.character()%>%
  str_replace(., "hetcount", "mMSS")

MSSCOMP_DF <- (PhenoDF[,c(1, which(grepl("^mMSS_complex_", colnames(PhenoDF))))])%>%
  melt(., id.vars = "id")

MSSCOMP_DF <- MSSCOMP_DF[paste0(MSSCOMP_DF$id, "_", MSSCOMP_DF$variable) %in% KIV,]



TBDF <- PhenoDF%>%
  select(SampID, AnyCHIP_init_filt, AnyCHIP_init_filt_VAF20, NOM_CHIP)%>%
  mutate(MultiHit = NOM_CHIP>1)%>%
  select(-NOM_CHIP)

TBDF$IsGOI <- TBDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])


MSSCOMP_DF_TP <- left_join(MSSCOMP_DF, TBDF, by = c("id" = "SampID"))

#Supplementary Figure4E
MSSCOMP_DF_TP%>%
  mutate(variable = str_remove(variable, "mMSS_complex_"))%>%
  mutate(variable = factor(variable, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                           labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(AnyCHIP_init_filt = factor(AnyCHIP_init_filt, levels = rev(c(FALSE, TRUE)),
                       labels = rev(c("No CHIP", "CHIP"))))%>%
  ggplot(aes(x = variable, y = value, fill = AnyCHIP_init_filt))+
  geom_half_violin(side = "r", scale = "width", position = position_dodge(width = 0.85))+
  geom_half_boxplot(side = "l", alpha = 0, position = position_dodge(width = 0.85))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  ylab("mMSS")+
  xlab("Complex")+
  coord_flip(ylim = c(0, 2))



nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF <- PhenoDF
  ForLOGDF <- ForLOGDF[((ForLOGDF[,i])>0),]
  
  ForLOGDF$mMSS <- ForLOGDF[,which(colnames(ForLOGDF) == str_replace(colnames(ForLOGDF)[i], "hetcount", "mMSS"))]
  ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
  
  Linear_multi <- lm(LogmMSS ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
  si_tmp <- summary(Linear_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)






#Supplementary Figure4F
MSSCOMP_DF_TP%>%
  filter(AnyCHIP_init_filt == TRUE)%>%
  mutate(variable = str_remove(variable, "mMSS_complex_"))%>%
  mutate(variable = factor(variable, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                           labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(AnyCHIP_init_filt_VAF20 = factor(AnyCHIP_init_filt_VAF20, levels = rev(c(FALSE, TRUE)),
                                    labels = rev(c("2%<VAF<20%", "VAF>20%"))))%>%
  ggplot(aes(x = variable, y = value, fill = AnyCHIP_init_filt_VAF20))+
  geom_half_violin(side = "r", scale = "width", position = position_dodge(width = 0.85))+
  geom_half_boxplot(side = "l", alpha = 0, position = position_dodge(width = 0.85))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  ylab("mMSS")+
  xlab("Complex")+
  coord_flip(ylim = c(0, 2))







nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF <- PhenoDF%>%
    filter(AnyCHIP_init_filt == TRUE)
  ForLOGDF <- ForLOGDF[((ForLOGDF[,i])>0),]
  
  ForLOGDF$mMSS <- ForLOGDF[,which(colnames(ForLOGDF) == str_replace(colnames(ForLOGDF)[i], "hetcount", "mMSS"))]
  ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
  
  Linear_multi <- lm(LogmMSS ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
  si_tmp <- summary(Linear_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)






#Supplementary Figure4G
MSSCOMP_DF_TP%>%
  filter(AnyCHIP_init_filt == TRUE)%>%
  mutate(variable = str_remove(variable, "mMSS_complex_"))%>%
  mutate(variable = factor(variable, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                           labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(MultiHit = factor(MultiHit, levels = rev(c(FALSE, TRUE)),
                                          labels = rev(c("Single mutation", "Multiple mutations"))))%>%
  ggplot(aes(x = variable, y = value, fill = MultiHit))+
  geom_half_violin(side = "r", scale = "width", position = position_dodge(width = 0.85))+
  geom_half_boxplot(side = "l", alpha = 0, position = position_dodge(width = 0.85))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  ylab("mMSS")+
  xlab("Complex")+
  coord_flip(ylim = c(0, 2))








nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF <- PhenoDF%>%
    filter(AnyCHIP_init_filt == TRUE)
  ForLOGDF <- ForLOGDF[((ForLOGDF[,i])>0),]
  
  ForLOGDF$mMSS <- ForLOGDF[,which(colnames(ForLOGDF) == str_replace(colnames(ForLOGDF)[i], "hetcount", "mMSS"))]
  ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
  
  Linear_multi <- lm(LogmMSS ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
  si_tmp <- summary(Linear_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)







#Supplementary Figure4H
MSSCOMP_DF_TP%>%
  filter(AnyCHIP_init_filt == TRUE)%>%
  mutate(variable = str_remove(variable, "mMSS_complex_"))%>%
  mutate(variable = factor(variable, levels = rev(c("I", "III", "IV", "V", "DLOOP", "RRNA", "TRNA")),
                           labels = rev(c("I", "III", "IV", "V", "D-loop", "rRNA", "tRNA"))))%>%
  mutate(IsGOI = factor(IsGOI, levels = rev(c(FALSE, TRUE)),
                                          labels = rev(c("Other CHIP", "Spliceosome"))))%>%
  ggplot(aes(x = variable, y = value, fill = IsGOI))+
  geom_half_violin(side = "r", scale = "width", position = position_dodge(width = 0.85))+
  geom_half_boxplot(side = "l", alpha = 0, position = position_dodge(width = 0.85))+
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#E64B35FF", "#595959FF"))+
  ylab("mMSS")+
  xlab("Complex")+
  coord_flip(ylim = c(0, 2))







nch <- c()
pch <- c()

for(i in which(grepl("^hetcount_complex_", colnames(PhenoDF)))){
  
  ForLOGDF <- PhenoDF%>%
    filter(AnyCHIP_init_filt == TRUE)
  ForLOGDF$IsGOISplice <- ForLOGDF$SampID %in% (fin_df$SampID[fin_df$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])
  ForLOGDF <- ForLOGDF[((ForLOGDF[,i])>0),]
  
  ForLOGDF$mMSS <- ForLOGDF[,which(colnames(ForLOGDF) == str_replace(colnames(ForLOGDF)[i], "hetcount", "mMSS"))]
  ForLOGDF$LogmMSS <- log10(ForLOGDF$mMSS + 1)
  
  Linear_multi <- lm(LogmMSS ~ as.numeric(IsGOISplice) + rms::rcs(age, df = 4) + sex + smk_ever + prev_cancer_yn, data = ForLOGDF)
  si_tmp <- summary(Linear_multi)
  
  nch <- c(nch, colnames(ForLOGDF)[i])
  pch <- c(pch, si_tmp$coefficients[,4][[2]])
  
  print(colnames(ForLOGDF)[i])
}

data.frame(nch, pch)




