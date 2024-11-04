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

GeneralTheme <- theme(axis.title = element_text(color = "black", size = 16),
                      axis.text = element_text(color = "black", size = 14),
                      strip.text = element_text(size = 14, color = "black"),
                      axis.ticks = element_line(color = "black"),
                      legend.title = element_text(color = "black", size = 16),
                      legend.text = element_text(color = "black", size = 14),
                      legend.position = "top")


# Figure 1 ARIC -------------------------------------------------------------

ARIC_Freq_Tab <- Base_CHIP_vars%>%
  distinct(SampID, Gene.refGene)%>%
  pull(Gene.refGene)%>%
  table()%>%
  data.frame()%>%
  rename("Var1" = ".")%>%
  arrange(-Freq)

#Supplementary Figure 1B
ARIC_Freq_Tab%>%
  head(n=15)%>%
  ggplot(aes(x = reorder(Var1, Freq), y = Freq))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  xlab("Gene")+
  ylab("Number of individuals")+
  coord_flip()


#Figure 1B
ARIC_Freq_Tab%>%#additional genes like ZBTB33 are just to keep graph proportions; edited with SRSF2 and removed the bar
  filter(Var1 %in% c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "ZBTB33","SRSF2", "SF3B1", "U2AF1"))%>%
  mutate(Var1 = factor(Var1, levels = rev(c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "ZBTB33", "SRSF2", "SF3B1", "U2AF1"))))%>%
  ggplot(aes(x = Var1, y = Freq))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  xlab("Gene")+
  ylab("Number of individuals")+
  coord_flip()

FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "ZBTB33", "SRSF2", "SF3B1", "U2AF1")

#Figure 1D
Base_CHIP_vars%>%
  filter(Gene.refGene %in% FreqGenesToAssess)%>%
  mutate(Gene.refGene = factor(Gene.refGene, levels = rev(FreqGenesToAssess)))%>%
  ggplot(aes(x = Gene.refGene, y = VAF))+
  geom_violin(scale = "width")+
  geom_boxplot(alpha = 0, width = 0.15)+
  theme_classic()+
  GeneralTheme+
  xlab("Gene")+
  ylab("VAF")+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = paste0(c(0, 25, 50, 75, 100), "%"),
                     limits = c(0, 100))+
  coord_flip()


CHIP_var_int <- Base_CHIP_vars
CHIP_var_int$GOI <- CHIP_var_int$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")
wilcox.test(CHIP_var_int$VAF ~ CHIP_var_int$GOI)



#Figure 1F
table(Base_CHIP_vars$SampID)%>%
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


# Figure 2 ARIC-------------------------------------------------------------

HCCOMP_DF <- (ARIC_CANCER[,c(2, which(grepl("^hetcount_complex_", colnames(ARIC_CANCER))))])%>%
  melt(., id.vars = "id")

#Figure2B
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


KIV <- HCCOMP_DF%>%
  filter(value>0)%>%
  mutate(CombVar = paste0(id, "_", variable))%>%
  pull(CombVar)%>%
  as.character()%>%
  str_replace(., "hetcount", "mMSS_ukb")

MSSCOMP_DF <- (ARIC_CANCER[,c(2, which(grepl("^mMSS_ukb_complex_", colnames(ARIC_CANCER))))])%>%
  melt(., id.vars = "id")

MSSCOMP_DF <- MSSCOMP_DF[paste0(MSSCOMP_DF$id, "_", MSSCOMP_DF$variable) %in% KIV,]

#Figure2D
MSSCOMP_DF%>%
  mutate(variable = str_remove(variable, "mMSS_ukb_complex_"))%>%
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


MSSCOMP_DF$TRvar <- MSSCOMP_DF$variable %in% c("mMSS_ukb_complex_RRNA", "mMSS_ukb_complex_TRNA")
wilcox.test(MSSCOMP_DF$value ~ MSSCOMP_DF$TRvar)

#Figure2F
table(ARIC_CANCER$count_het)%>%
  data.frame()%>%
  mutate(Var1 = as.numeric(as.character(Var1)))%>%
  filter(Var1>0)%>%
  ggplot(aes(x = Freq, y = reorder(Var1, -Var1)))+
  geom_col()+
  theme_classic()+
  GeneralTheme+
  ylab("Number of heteroplasmies")+
  xlab("Number of individuals")


# Figure 3B ARIC ----------------------------------------------------------

prop.table(table(ARIC_CANCER$count_het>0, ARIC_CANCER$AnyCHIP_init_filt))*100

# Figure 3D ARIC-------------------------------------------------------------

FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "SRSF2", "SF3B1", "U2AF1")

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
ORch <- c()
SEch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()

for(i in c("AnyCHIP_init_filt", "AnyCHIP_init_filt_Small", "AnyCHIP_init_filt_VAF20", "NOMsingle", "NOMmulti")){
  
  int_df <- ARIC_CANCER
  
  int_df$HCov0 <- int_df$het_count>0
  
  int_CHIP <- Base_CHIP_vars
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
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  #getting the gene section

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid[int_CHIP$Gene == g]
    
    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid
    }
    

    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      

      
      ForLOGDF <- int_df%>%filter(het_count==0 | het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)

      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
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




si_ct <- with(ARIC_CANCER, table((AnyCHIP_init_filt==FALSE), (het_count>0)))%>%
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
  xlab("Heteroplasmy prevalence")


ForLOGDF <- ARIC_CANCER
Logistic_multi <- glm(as.numeric(het_count>0) ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)




ForLOGDF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt == TRUE)

Logistic_multi <- glm(as.numeric(het_count>0) ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



Logistic_multi <- glm(as.numeric(het_count>0) ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)




ForLOGDF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt == TRUE)

ForLOGDF$IsGOI <- ForLOGDF$aricid %in% (Base_CHIP_vars$aricid[Base_CHIP_vars$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])

Logistic_multi <- glm(as.numeric(het_count>0) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



# Figure 3F ARIC-------------------------------------------------------------
FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "SRSF2", "SF3B1", "U2AF1")

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()

for(i in c("AnyCHIP_init_filt", "AnyCHIP_init_filt_Small", "AnyCHIP_init_filt_VAF20", "NOMsingle", "NOMmulti")){
  
  int_df <- ARIC_CANCER%>%
    filter(het_count>0)
  
  int_df$HCov0 <- int_df$het_count>1
  
  int_CHIP <- Base_CHIP_vars
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
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  #getting the gene section

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      
      ForLOGDF <- int_df%>%filter(het_count==1 | het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
      
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
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



si_ct <- with(ARIC_CANCER%>%filter(het_count>0), table((AnyCHIP_init_filt==FALSE), (het_count>1)))%>%
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







ForLOGDF <- ARIC_CANCER%>%
  filter(het_count>0)
Logistic_multi <- glm(as.numeric(het_count>1) ~ as.numeric(AnyCHIP_init_filt)+ rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

ForLOGDF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt == TRUE & het_count>0)

Logistic_multi <- glm(as.numeric(het_count>1) ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



Logistic_multi <- glm(as.numeric(het_count>1) ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



ForLOGDF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt == TRUE & het_count>0)

ForLOGDF$IsGOI <- ForLOGDF$aricid %in% (Base_CHIP_vars$aricid[Base_CHIP_vars$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])

Logistic_multi <- glm(as.numeric(het_count>1) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



# Figure 3H ARIC------------------------------------------



FreqGenesToAssess <- c("DNMT3A", "TET2", "ASXL1", "ZBTB33", "TP53", "YLPM1", "SF3B1", "U2AF1")#ZBTB33 and YLPM1 were used to keep figure proportions in cases with no heteroplasmy, their specific plots were further removed

ARIC_LIN_DF <- data.frame()



for(i in c("AnyCHIP_init_filt", "AnyCHIP_init_filt_Small", "AnyCHIP_init_filt_VAF20", "NOMsingle", "NOMmulti")){
  
  int_df <- ARIC_CANCER
  
  int_df$LogmMSS_ukb <- log10(int_df$mMSS_ukb+1)
  int_df$HCov0 <- int_df$het_count>0

  int_CHIP <- Base_CHIP_vars
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
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  #getting the gene section
  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid[int_CHIP$Gene == g]

    
    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid
    }
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]
      
      
      
      ForLOGDF <- int_df%>%filter(het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
      
      
      ToBindDf <- ForLOGDF%>%
        filter(IsGOI)%>%
        select(mMSS_ukb)
      
      if(length(ToBindDf$mMSS_ukb)==0){
        next
      }
      
      Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF)
      si_linear_sum <- summary(Linear_multi)
      
      ToBindDf$CHIPfiltch <- i
      ToBindDf$Genech <- g
      ToBindDf$hvarch <- hvar
      ToBindDf$pch <- (si_linear_sum$coefficients[2,4])
      
      ARIC_LIN_DF <- rbind(ARIC_LIN_DF, ToBindDf)
      
      rm(Linear_multi)
      rm(si_linear_sum)
      
    }
  }
}

ARIC_LIN_DF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt==FALSE & het_count>0)%>%
  select(mMSS_ukb)%>%
  mutate(CHIPfiltch = "AnyCHIP_init_filt", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(ARIC_LIN_DF, .)

ARIC_LIN_DF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt==FALSE & het_count>0)%>%
  select(mMSS_ukb)%>%
  mutate(CHIPfiltch = "AnyCHIP_init_filt_Small", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(ARIC_LIN_DF, .)

ARIC_LIN_DF <- ARIC_CANCER%>%
  filter(AnyCHIP_init_filt==FALSE & het_count>0)%>%
  select(mMSS_ukb)%>%
  mutate(CHIPfiltch = "AnyCHIP_init_filt_VAF20", Genech = "No CHIP", hvarch = "HCov0", pch = NA)%>%
  rbind(ARIC_LIN_DF, .)

ARIC_LIN_DF <- ARIC_LIN_DF%>%
  distinct(CHIPfiltch, Genech, hvarch, pch)%>%
  mutate(P_adj = p.adjust(pch, method = "BH"))%>%
  select(-pch)%>%
  left_join(ARIC_LIN_DF, ., by = c("CHIPfiltch", "Genech", "hvarch"))






ARIC_LIN_DF <- ARIC_LIN_DF%>%
  mutate(PlotGF = paste(CHIPfiltch, Genech))%>%
  mutate(PlotGF = factor(PlotGF, levels = rev(c("AnyCHIP_init_filt_Small No CHIP", "AnyCHIP_init_filt any", "AnyCHIP_init_filt_Small any", "AnyCHIP_init_filt_VAF20 any", "NOMsingle any", "NOMmulti any", paste0("AnyCHIP_init_filt ", FreqGenesToAssess))),
                         labels = rev(c("No CHIP", "CHIP", "2%<VAF<20%", "VAF>20%", "Single mutation", "Multiple mutations", FreqGenesToAssess))))%>%
  filter(!(is.na(PlotGF)))

ARIC_LIN_DF$SplVar <- case_when(ARIC_LIN_DF$PlotGF %in% c("No CHIP", "CHIP") ~ "wSC",
                               ARIC_LIN_DF$PlotGF %in% c("2%<VAF<20%", "VAF>20%") ~ "xSC",
                               ARIC_LIN_DF$PlotGF %in% c("Single mutation", "Multiple mutations") ~ "ySC",
                               TRUE ~ "zG")


ARIC_LIN_DF$Significant <- ifelse(ARIC_LIN_DF$PlotGF %in% c("No CHIP"), NA, "*")
ARIC_LIN_DF$Significant <- ifelse(is.na(ARIC_LIN_DF$Significant), "Not significant", "Significant")





ARIC_LIN_DF%>%
  ggplot(aes(x = PlotGF, y = log10(mMSS_ukb+1), fill = Significant))+
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




ForLOGDF <- ARIC_CANCER%>%filter(het_count>0)
ForLOGDF$LogmMSS_ukb <- log10(ForLOGDF$mMSS_ukb + 1)
Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(AnyCHIP_init_filt) + rms::rcs(age, df = 4)+het_count + gender + EVR_SMK_fin, data = ForLOGDF)
summary(Linear_multi)




ForLOGDF <- ARIC_CANCER%>%filter(het_count>0)%>%filter(AnyCHIP_init_filt == TRUE)
ForLOGDF$LogmMSS_ukb <- log10(ForLOGDF$mMSS_ukb + 1)
Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(AnyCHIP_init_filt_VAF20) + rms::rcs(age, df = 4)+het_count + gender + EVR_SMK_fin, data = ForLOGDF)
summary(Linear_multi)



Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(NOM_CHIP>1) + rms::rcs(age, df = 4)+het_count + gender + EVR_SMK_fin, data = ForLOGDF)
summary(Linear_multi)





ForLOGDF <- ARIC_CANCER%>%filter(het_count>0)%>%filter(AnyCHIP_init_filt == TRUE)
ForLOGDF$LogmMSS_ukb <- log10(ForLOGDF$mMSS_ukb + 1)
ForLOGDF$IsGOI <- ForLOGDF$aricid %in% (Base_CHIP_vars$aricid[Base_CHIP_vars$Gene.refGene %in% c("SRSF2", "SF3B1", "U2AF1")])
Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(IsGOI) + rms::rcs(age, df = 4)+het_count + gender + EVR_SMK_fin, data = ForLOGDF)
summary(Linear_multi)





# Figure 4B,D ARIC -----------------------------------------------------------

int_df <- ARIC_CANCER

m <- "AnyCHIP_init_filt"

int_df$BindingHetVAR <- int_df$het_count>0

int_df <- int_df[((int_df$AnyCHIP_init_filt==0)|(int_df[,which(colnames(int_df) == m)])),]
int_df <- int_df[((int_df$BindingHetVAR)|(int_df$het_count==0)),]



int_df$AnyCHIP_VAR <- (int_df[,which(colnames(int_df) == m)])
int_df$CHIP_VAR <- paste0((int_df$AnyCHIP_VAR),
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



fit.coxph <- coxph(surv_object ~ as.numeric(BindingHetVAR) + rms::rcs(age, df = 4)  + gender + EVR_SMK_fin, 
                   data = int_df)

fit.coxph <- coxph(surv_object ~ CHIP_VAR + rms::rcs(age, df = 4)  + gender + EVR_SMK_fin, 
                   data = int_df)



exp(fit.coxph[[1]])
exp(confint(fit.coxph))
summary(fit.coxph)

# Table1 ARIC -------------------------------------------------------------

t.test(ARIC_CANCER$age ~ ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$gender, ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$RACEGRP=="W", ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$EVR_SMK_fin, ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$Anemia_fin, ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$Thrombocytopenia_fin, ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$Neutropenia_fin, ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$Cytopenia_fin, ARIC_CANCER$AnyCHIP_init_filt)

t.test(ARIC_CANCER$MCV_fin ~ ARIC_CANCER$AnyCHIP_init_filt)

fisher.test(ARIC_CANCER$het_count>0, ARIC_CANCER$AnyCHIP_init_filt)


# Table2 ARIC -------------------------------------------------------------

t.test(ARIC_CANCER$age ~ ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$gender, ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$RACEGRP=="W", ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$EVR_SMK_fin, ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$Anemia_fin, ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$Thrombocytopenia_fin, ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$Neutropenia_fin, ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$Cytopenia_fin, ARIC_CANCER$het_count>0)

t.test(ARIC_CANCER$MCV_fin ~ ARIC_CANCER$het_count>0)

fisher.test(ARIC_CANCER$AnyCHIP_init_filt, ARIC_CANCER$het_count>0)



# Supplementary Figure2B ARIC ----------------------------------------------



FreqGenesToAssess <- as.character(unique(Base_CHIP_vars$Gene.refGene))

CHIPfiltch <- c()
Genech <- c()
hvarch <- c()
HetCHch <- c()
TotCHch <- c()
PercHet <- c()
pch <- c()

Medch <- c()



for(i in c("AnyCHIP_init_filt")){
  
  int_df <- ARIC_CANCER
  
  int_df$HCov0 <- int_df$het_count>0

  
  int_CHIP <- Base_CHIP_vars
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
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  
  
  if(i == "NOMsingle"){
    int_df <- int_df[(int_df$NOM_CHIP == 1)|(int_df$AnyCHIP_init_filt == FALSE),]
    int_CHIP <- int_CHIP[int_CHIP$aricid %in% int_df$aricid,]
    
  }
  

  for (g in c("any", as.character(unique(FreqGenesToAssess)))){
    print(i)
    print(g)
    
    int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid[int_CHIP$Gene == g]

    
    
    if(g == "any"){
      int_df$IsGOI <- int_df$aricid %in% int_CHIP$aricid
    }
    
    
    for (hvar in c("HCov0")){
      
      int_df$het_VOI <- int_df[,which(colnames(int_df) == hvar)]

      
      ForLOGDF <- int_df%>%filter(het_count==0 | het_VOI)%>%filter(AnyCHIP_init_filt == 0| IsGOI)
      
      
      si_ct <- with(ForLOGDF, table(IsGOI, het_VOI))%>%
        data.frame()%>%
        filter(IsGOI == TRUE)%>%
        mutate(s = sum(Freq))%>%
        filter(het_VOI == TRUE)%>%
        mutate(P = Freq/s*100)
      
      
      
      Logistic_multi <- glm(het_VOI ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
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
                               pull(mMSS_ukb)%>%
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

FreqGeneCHIPhet_df$UKBalso <- FreqGeneCHIPhet_df$Genech %in% read.csv("GOI_UKB_OutAlt3.csv")$x
FreqGeneCHIPhet_df$TotCHch <- ifelse(FreqGeneCHIPhet_df$Genech == "any", 10000, FreqGeneCHIPhet_df$TotCHch)#for consistency between cohorts, edited to the correct cathegory in the final form
FreqGeneCHIPhet_df$SzCat <- cut(FreqGeneCHIPhet_df$TotCHch, c(0, 10, 100, 1000, 100000))
options(ggrepel.max.overlaps = Inf)
FreqGeneCHIPhet_df%>%
  filter(UKBalso)%>%
  filter(TotCHch>3)%>%
  filter(HetCHch>1)%>%
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






# Supplementary Figure3 ARIC --------------------------------------------

PhenoDF <- ARIC_CANCER
PhenoDF$het_yn <- as.numeric(PhenoDF$het_count>0)
PhenoDF <- PhenoDF[!is.na(PhenoDF$Cytopenia_fin),]
PhenoDF$ccus <- case_when((PhenoDF$Cytopenia_fin&PhenoDF$AnyCHIP_init_filt) ~ "ccus",
                          ((PhenoDF$Cytopenia_fin==FALSE)&PhenoDF$AnyCHIP_init_filt) ~ "chip",
                          TRUE ~ "none")


#Supplementary Figure 3B
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
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  xlab("Heteroplasmy prevalence")+
  ylab("")+
  geom_text(aes(x = 1, y = ccus, label = Freq), color = "white", hjust = 0, size = 6)


ForLOGDF <- PhenoDF%>%
  filter(ccus != "ccus")


ForLOGDF$IsGOI <- as.numeric(ForLOGDF$ccus=="chip")
Logistic_multi <- glm(as.numeric(het_count>0) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)


ForLOGDF <- PhenoDF%>%
  filter(ccus != "none")


ForLOGDF$IsGOI <- as.numeric(ForLOGDF$ccus=="ccus")
Logistic_multi <- glm(as.numeric(het_count>0) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)



#Supplementary Figure 3D
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
  theme_classic()+
  GeneralTheme+
  scale_fill_manual(values = c("#595959FF", "#E64B35FF"))+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = 0:10*10,
                     labels = paste0(0:10*10, "%"))+
  xlab("Prevalence of multiple heteroplasmies")+
  ylab("")+
  geom_text(aes(x = 1, y = ccus, label = Freq), color = "white", hjust = 0, size = 6)



ForLOGDF <- PhenoDF%>%
  filter(ccus != "ccus")%>%
  filter(het_count>0)


ForLOGDF$IsGOI <- as.numeric(ForLOGDF$ccus=="chip")
Logistic_multi <- glm(as.numeric(het_count>1) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)

ForLOGDF <- PhenoDF%>%
  filter(ccus != "none")%>%
  filter(het_count>0)


ForLOGDF$IsGOI <- as.numeric(ForLOGDF$ccus=="ccus")
Logistic_multi <- glm(as.numeric(het_count>1) ~ as.numeric(IsGOI) + rms::rcs(age, df = 4) + gender + EVR_SMK_fin, data = ForLOGDF, family = binomial)
summary(Logistic_multi)


PhenoDF%>%
  filter(het_yn==1)%>%
  filter(ccus == "none")%>%
  pull(mMSS_ukb)%>%
  quantile()

#Supplementary Figure 3F
PhenoDF%>%
  filter(het_yn==1)%>%
  mutate(ccus = factor(ccus, levels = rev(c("none", "chip", "ccus")),
                       labels = rev(c("No CHIP", "CHIP\nwithout\ncytopenia", "CCUS"))))%>%
  mutate(Significant = ccus != "No CHIP")%>%
  ggplot(aes(x = ccus, y = log10(mMSS_ukb+1), fill = Significant))+
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

ForLOGDF <- PhenoDF%>%
  filter(ccus != "ccus")%>%
  filter(het_count>0)

ForLOGDF$IsGOI <- as.numeric(ForLOGDF$ccus=="chip")

ForLOGDF$LogmMSS_ukb <- log10(ForLOGDF$mMSS_ukb + 1)
Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(IsGOI) + rms::rcs(age, df = 4)+het_count + gender + EVR_SMK_fin, data = ForLOGDF)
summary(Linear_multi)


ForLOGDF <- PhenoDF%>%
  filter(ccus != "none")%>%
  filter(het_count>0)

ForLOGDF$IsGOI <- as.numeric(ForLOGDF$ccus=="ccus")

ForLOGDF$LogmMSS_ukb <- log10(ForLOGDF$mMSS_ukb + 1)
Linear_multi <- lm(LogmMSS_ukb ~ as.numeric(IsGOI) + rms::rcs(age, df = 4)+het_count + gender + EVR_SMK_fin, data = ForLOGDF)
summary(Linear_multi)




