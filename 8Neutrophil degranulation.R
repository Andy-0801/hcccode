#### Neutrophil degranulation

## xCell for immune cells
library(xCell)
library(dplyr)
library(ComplexHeatmap)
exprMatrix <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1, rowNames = T)%>%as.data.frame()

xcell.data <- as.data.frame(xCellAnalysis(exprMatrix))

## NMF clusters heatmap in Metabolome and lipidome
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

nmf_list <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/NMF_C1-2-3_DALs.xlsx", sheet = 2, rowNames = F)
nmf_list <- as.data.frame(nmf_list[,1])
names(nmf_list) <- "Name" 

nmf_list <- nmf_list%>% distinct(Name, .keep_all = T)

prodata <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 3, rowNames = F)%>%as.data.frame()

dep.data <- merge(nmf_list, prodata, by = "Name", sort = FALSE)
dep.data <- column_to_rownames(dep.data, var = "Name")

plotdata <- as.data.frame(t(scale(log2(1+t(dep.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 2, rowNames = T)%>%as.data.frame()

sub = HeatmapAnnotation(Subtypes = clinicaldata$Subtypes,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Subtypes = c("C1" = "#2874C5",
                                                "C2" = "#EABF00",
                                                "C3" = "#C6524A")))

col_fun = circlize::colorRamp2(c(-4,-2, 0, 2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(123)
p <- Heatmap(plotdata, 
             col = col_fun, column_split = clinicaldata$Subtypes,
             top_annotation = sub,
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F, 
             cluster_rows = F, 
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscores"))

##  Neutrophil degranulation corrplot
library(dplyr)
library(ggplot2)
library(ggpubr)
data <- openxlsx::read.xlsx("../data/NMF_C3_neu-pufa.xlsx", sheet = 1, rowNames = T)%>%as.data.frame()

plotdata <- as.data.frame(scale(log2(data)))

p <- ggplot(plotdata, aes(y=ALOX5, x=MMP9)) +
  geom_point(color="black") + stat_smooth(method="lm")+
  stat_cor(data=plotdata, method = "pearson")

m <- ggplot(plotdata, aes(y=ALOX5, x=MPO)) +
  geom_point(color="black") + stat_smooth(method="lm")+
  stat_cor(data=plotdata, method = "pearson")

i <- ggplot(plotdata, aes(y=ALOX5, x=ITGAM)) +
  geom_point(color="black") + stat_smooth(method="lm")+
  stat_cor(data=plotdata, method = "pearson")

g <- ggarrange(ggarrange(p, m, i,
                         labels = c("", "", ""),
                         ncol = 3, nrow = 1))
## PUFA metabolism
library(dplyr)
library(tibble)
library(circlize)
library(ComplexHeatmap)
pufa.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolome_PUFA.xlsx", sheet = 1, rowNames = T)%>%as.data.frame()

pufa.data <- as.data.frame(t(pufa.data[,-c(1,13:19)]))

plotdata <- as.data.frame(t(scale(log2(1+t(pufa.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3, rowNames = T)%>%as.data.frame()

sub = HeatmapAnnotation(Subtypes = clinicaldata$Subtypes,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Subtypes = c("NAT"= "grey",
                                                "C1" = "#2874C5",
                                                "C2" = "#EABF00",
                                                "C3" = "#C6524A")))

col_fun = circlize::colorRamp2(c(-4,-2, 0, 2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(111)
p <- Heatmap(plotdata,
             col = col_fun,column_split = clinicaldata$Subtypes,
             top_annotation = sub,
             show_row_names = T,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F, 
             cluster_rows = F, 
             row_names_gp = gpar(fontsize = 7.5),
             heatmap_legend_param = list(title = "Zscores"))

## Metaboanalysis anotation
library(ggplot2)
library(dplyr)
plotdata <- read.csv("../data/2022ProteomeMetabolism/NMF_C3_pathway_results.csv", header = T)
plotdata <- plotdata[1:6,]

p <- ggplot(plotdata, aes(y = reorder(factor(Name),logP),x=logP)) + 
  geom_bar(stat = 'identity', position = 'dodge', fill = "#C6524A") +
  theme_classic()+scale_x_reverse()+
  scale_y_discrete(position = "right")+
  xlab("-log10 (p-value)")+ylab("")

