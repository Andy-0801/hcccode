####Surivival analysis in CPTAC data

##Ribosome-genes
library(survminer)
library(survival)
library(dplyr)
library(tibble)
rb.list <- openxlsx::read.xlsx("../data/NMF_Ribosome b_genes.xlsx", sheet = 2,rowNames = F)%>%as.data.frame()
rb.list <- as.data.frame(rb.list[,1])
names(rb.list) <- "Symbol" 

rb.list <- rb.list%>% distinct(Symbol, .keep_all = T)

cptac.data <- openxlsx::read.xlsx("../data/S049_Liver_Cancer_Gao2019/hcc_proteome_TTexp.xlsx", sheet = 1, rowNames = F)

rb.data <- merge(rb.list, cptac.data, by="Symbol")
rb.data <- as.data.frame(t(column_to_rownames(rb.data, var = "Symbol")))
rb.data <- rownames_to_column(rb.data, var = "TumorID")

CPTAC.clinical <- openxlsx::read.xlsx("../data/S049_Liver_Cancer_Gao2019/Zhou_Liver_Cancer_ProtemoeMeta.xlsx", sheet = 2, rowNames = F)

surv.data <- merge(CPTAC.clinical[,c(1,24:27)], rb.data, by="TumorID")

cutoff <- surv_cutpoint(surv.data, time = 'Survival_time', 
                        event = 'Survival_test', 
                        variables = c("XPO1","RPRD1B","RAN","CSNK2A1"))

summary(cutoff)

ncut <- cut(surv.data$XPO1, breaks = c(0, 1.798, Inf), 
            labels = c('Low', 'High'))
scut <- cut(surv.data$RPRD1B, breaks = c(0, 1.372, Inf), 
            labels = c('Low', 'High'))
pcut <- cut(surv.data$RAN, breaks = c(0, 2.185, Inf), 
            labels = c('Low', 'High'))
ccut <- cut(surv.data$CSNK2A1, breaks = c(0, 1.68, Inf), 
            labels = c('Low', 'High'))

fit <- survfit(Surv(Survival_time, Survival_test)~scut, data = surv.data)

p <- ggsurvplot(fit, pval=TRUE, surv.median.line = "hv", 
                palette=c("blue", "red"),
                legend.title="RPRD1B", legend.lab=c('Low','High'),
                legend=c(0.85,0.15),
                ggtheme = theme_bw(),
                title="Overall Survival in CPTAC")
## ND-genes
ne.list <- openxlsx::read.xlsx("../data/NMF_C1-2-3_DEPs.xlsx", sheet = 3,rowNames = F)
ne.list <- as.data.frame(ne.list[,1])
names(ne.list) <- "Symbol" 

ne.list <- ne.list%>% distinct(Symbol, .keep_all = T)

cptac.data <- openxlsx::read.xlsx("../data/S049_Liver_Cancer_Gao2019/hcc_proteome_TTexp.xlsx", sheet = 1, rowNames = F)

ne.data <- merge(ne.list, cptac.data, by="Symbol")
ne.data <- as.data.frame(t(column_to_rownames(ne.data, var = "Symbol")))
ne.data <- rownames_to_column(ne.data, var = "TumorID")

CPTAC.clinical <- openxlsx::read.xlsx("../data/S049_Liver_Cancer_Gao2019/Zhou_Liver_Cancer_ProtemoeMeta.xlsx", sheet = 2, rowNames = F)

surv.data <- merge(CPTAC.clinical[,c(1,24:27)], ne.data, by="TumorID")

cutoff <- surv_cutpoint(surv.data, time = 'Survival_time', 
                        event = 'Survival_test', 
                        variables = c("ICAM1","CEACAM8","ITGAM","S100A8",
                                      "ITGB2","LCN2","MPO","TIMP1","VCAM1",
                                      "MMP9","CD63","S100A9","ALOX5"))

summary(cutoff)

mcut <- cut(surv.data$MMP9, breaks = c(0, 1.900, Inf), 
            labels = c('Low', 'High'))
ccut <- cut(surv.data$CD63, breaks = c(0, 3.652, Inf), 
            labels = c('Low', 'High'))
scut <- cut(surv.data$S100A9, breaks = c(0, 0.222, Inf), 
            labels = c('Low', 'High'))
acut <- cut(surv.data$ALOX5, breaks = c(0, 5.196, Inf), 
            labels = c('Low', 'High'))

cut1 <- cut(surv.data$ICAM1, breaks = c(0, 0.938, Inf), 
            labels = c('Low', 'High'))
cut2 <- cut(surv.data$CEACAM8, breaks = c(0, 2.165, Inf), 
            labels = c('Low', 'High'))
cut3 <- cut(surv.data$ITGAM, breaks = c(0, 6.538, Inf), 
            labels = c('Low', 'High'))
cut4 <- cut(surv.data$S100A8, breaks = c(0, 0.439, Inf), 
            labels = c('Low', 'High'))
cut5 <- cut(surv.data$ITGB2, breaks = c(0, 2.480, Inf), 
            labels = c('Low', 'High'))
cut6 <- cut(surv.data$LCN2, breaks = c(0, 0.751, Inf), 
            labels = c('Low', 'High'))
cut7 <- cut(surv.data$MPO, breaks = c(0, 2.626, Inf), 
            labels = c('Low', 'High'))
cut8 <- cut(surv.data$TIMP1, breaks = c(0, 14.636, Inf), 
            labels = c('Low', 'High'))
cut9 <- cut(surv.data$VCAM1, breaks = c(0, 0.097, Inf), 
            labels = c('Low', 'High'))

fit <- survfit(Surv(Survival_time, Survival_test)~cut9, data = surv.data)

p <- ggsurvplot(fit, pval=TRUE, surv.median.line = "hv", 
                palette=c("blue", "red"),
                legend.title="VCAM1", legend.lab=c('Low','High'),
                legend=c(0.85,0.15),
                title="Overall Survival in CPTAC")

a <- ggsurvplot(fit, pval=TRUE, surv.median.line = "hv", 
                palette=c("blue", "red"),
                legend.title="ALOX5", legend.lab=c('Low','High'),
                legend=c(0.85,0.15),
                ggtheme = theme_bw(),
                title="Overall Survival in CPTAC")
