## metabolite-interacting proteins (MIPros)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggpubr)
prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1,rowNames = F)

proteinlist <- openxlsx::read.xlsx("../data/2022proteome/MIPros_list.xlsx", sheet = 2,rowNames = F)

MIPros_exp <- merge(proteinlist, prodata, by="Name")
MIPros_exp <- column_to_rownames(MIPros_exp, var = "Name")

plotdata <- as.data.frame(t(log2(MIPros_exp+1)))
groups <- as.data.frame(c(rep("NAT",39), rep("TT",39))) 
names(groups) <- "Groups"

plotdata <- cbind(groups, plotdata)

data <- melt(plotdata, id.vars = "Groups")

p <- ggplot(data, aes(x=Groups, y=value, color=Groups))+
  geom_boxplot(outlier.size = 0.5,alpha=0.5)+
  scale_color_manual(values=c(NAT = "#2c83b7", TT = "#c2395a"))+
  stat_compare_means()+
  theme_classic()+ 
  ylab("Log2(Normalized expression value)")+
  xlab("Proteins assigned to metabolic pathways")+
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

## DEP in Tumor vs NAT
library(dplyr)
library(limma)
prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1,rowNames = T)

groups <- c(rep("NAT",39), rep("TT",39))

design<- model.matrix(~0+factor(groups))
colnames(design)<- levels(factor(groups))
#voom transform
v <- voom(prodata, design, plot=TRUE)

fit <- lmFit(v, design)

cont.matrix <- makeContrasts('TT-NAT', levels = design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
write.csv(output, file = "../results/Proteome_DEPs_limmavoom.csv",row.names = T)


#metabolome anaysis
library(dplyr)
library(limma)
lipid_data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites2023.xlsx", sheet = 1,rowNames = T)%>% as.data.frame()

group <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3, rowNames = T)
design<- model.matrix(~0+factor(group$Groups))
colnames(design)<- levels(factor(group$Groups))  

v <- voom(lipid_data[,-1], design, plot=TRUE)
v <- voom(lipid_data[,-c(1:12)], design, plot=TRUE)

fit <- lmFit(v, design)

cont.matrix <- makeContrasts('TT-NAT', levels = design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
write.csv(output, file = "../results/Metabolite_limmavoom.csv",row.names = T)

##heatmap of DEPs in TT vs NAT
library(ComplexHeatmap)
library(tibble)
library(dplyr)
DEP_list <- openxlsx::read.xlsx("../data/Proteome_DEGs_limmavoom.xlsx", sheet = 3)
DEP_list <- DEP_list[,1]%>%as.data.frame()
names(DEP_list) <- "Name"

prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1,rowNames = F)

#prodata <- rownames_to_column(prodata, var = "Name")

DEP_exp <- merge(DEP_list, prodata, by="Name", sort = F)
DEP_exp <- column_to_rownames(DEP_exp, var = "Name")

plotdata <- as.data.frame(t(scale(log2(t(DEP_exp)))))

groups <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3)

mark_gene <- openxlsx::read.xlsx("../data/Proteome_DEGs_limmavoom.xlsx", sheet = 4)%>%as.data.frame()

row_anno <-  rowAnnotation(link = anno_mark(at = which(rownames(plotdata) 
                                                       %in% mark_gene$Name), 
                                            labels = mark_gene$Name,
                                            labels_gp = gpar(fontsize = 6)))

ht = HeatmapAnnotation(Cluster = groups$Groups,
                       show_legend = F,
                       col = list(Cluster = c("NAT" = "#2c83b7",
                                              "TT"   = "#c2395a")))
col_fun = circlize::colorRamp2(c(-4, -2, 0, 2, 4),
                               c("#2c83b7","deepskyblue1", "white",
                                 "red", "#c2395a"))

p <- Heatmap(plotdata, 
             col = col_fun,
             top_annotation = ht,
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F,
             cluster_rows = F,
             heatmap_legend_param = list(title = "Zscore"))
p+row_anno

##heatmap of DAMs in TT vs NAT

library(ComplexHeatmap)
library(tibble)
metabolites_exp <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites2.xlsx", sheet = 2)
metabolites_list <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites_limmavoom.xlsx", sheet = 2)

metabolites_list <- as.data.frame(metabolites_list[,1])
names(metabolites_list) <- "Name"

DAM_exp <- merge(metabolites_list, metabolites_exp, sort = F, by="Name")
DAM_exp <- column_to_rownames(DAM_exp, var = "Name")

plotdata <- as.data.frame(t(scale(t(log2(DAM_exp+1)))))

groups <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3)

ht = HeatmapAnnotation(Cluster = groups$Groups,
                       show_legend = F,
                       col = list(Cluster = c("NAT" = "#2c83b7",
                                              "TT" = "#c2395a")))

col_fun = circlize::colorRamp2(c(-4, -2, 0, 2, 4),
                               c("#2c83b7","deepskyblue1", "white",
                                 "red", "#c2395a"))
marker <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites_limmavoom.xlsx", sheet = 3)%>%as.data.frame()

row_anno <-  rowAnnotation(link = anno_mark(at = which(rownames(plotdata) 
                                                       %in% marker$Name), 
                                            labels = marker$Name,
                                            labels_gp = gpar(fontsize = 6)))

p <- Heatmap(plotdata, 
             col = col_fun,
             top_annotation = ht, 
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F,
             cluster_rows = F,
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscore", 
                                         direction = "horizontal"))

p+row_anno

## heatmap of DALs in TT vs NAT
library(ComplexHeatmap)
library(tibble)
lipid_exp <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 1)
lipid_list <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipids_DAL_VIP.xlsx", sheet = 2)

lipid_list <- as.data.frame(lipid_list[,1])
names(lipid_list) <- "Name"

DAL_exp <- merge(lipid_list, lipid_exp, sort = F, by="Name")
DAL_exp <- column_to_rownames(DAL_exp, var = "Name")

plotdata <- as.data.frame(t(scale(t(log2(DAL_exp[,-1]+1)))))

groups <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3)

ht = HeatmapAnnotation(Cluster = groups$Groups,
                       show_legend = F,
                       col = list(Cluster = c("NAT" = "#2c83b7",
                                              "TT" = "#c2395a")))

col_fun = circlize::colorRamp2(c(-4, -2, 0, 2, 4),
                               c("#2c83b7","deepskyblue1", "white",
                                 "red", "#c2395a"))
marker <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipids_DAL_VIP.xlsx", sheet = 3)%>%as.data.frame()

row_anno <-  rowAnnotation(link = anno_mark(at = which(rownames(plotdata) 
                                                       %in% marker$Name), 
                                            labels = marker$Name,
                                            labels_gp = gpar(fontsize = 6)))

p <- Heatmap(plotdata, 
             col = col_fun,
             top_annotation = ht, 
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F,
             cluster_rows = F,
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscore", 
                                         direction = "horizontal"))

p+row_anno

## Metascale anotation

library(ggplot2)

Tumor.data <- openxlsx::read.xlsx("../data/2022proteome/TTvsNAT-metascape.xlsx", sheet = 1)

tt <- ggplot(Tumor.data, aes(y = reorder(factor(Description), `-Log(q-value)`),
                             x = `-Log(q-value)`)) + 
  geom_bar(stat = 'identity', fill = '#c2395a', width = 0.7)+
  xlab("")+ylab("-log10(q-value)")+
  theme_classic() +
  theme(legend.position ="",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

NAT.data <- openxlsx::read.xlsx("../data/2022proteome/TTvsNAT-metascape.xlsx", sheet = 2)

nat <- ggplot(NAT.data, aes(x = reorder(factor(Description), `-Log(q-value)`),
                            y = `-Log(q-value)`)) + 
  geom_bar(stat = 'identity', fill = "#2c83b7", width = 0.7) +
  xlab("")+ylab("-log10(p-value)")+
  coord_flip()+
  theme_classic() +
  theme(legend.position ="",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

c2 <- ggplot(cluster2, aes(x = reorder(factor(Description), LogP), y = LogP)) + 
  geom_bar(stat = 'identity', fill = '#0f4392', width = 0.7)+
  xlab("")+ylab("-log10(p-value)")+
  coord_flip()+
  theme_classic() +
  theme(legend.position ="",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))


#### Bubble plot of MSEA in TTvsNAT

library(ggplot2)

nat.data <- read.csv("../data/2022ProteomeMetabolism/02NAT_msea_ora_result.csv", header = T)

plotdata <- nat.data[c(1:13),c(1,8,9)]

p <- ggplot(plotdata,aes(x = LogP, 
                         y = EnrichmentRatio, 
                         size = EnrichmentRatio,
                         colour=LogP)) +
  geom_point(shape = 16) +               
  labs(x = "LogP", y = "EnrichmentRatio")+          
  scale_colour_gradient2(name="-log10(P-value)",
                         low = "#ffbd39",
                         mid = "#e61c5d",
                         high = "#930077",
                         midpoint = 2)+
  theme_bw()          


tt.data <- read.csv("../data/2022ProteomeMetabolism/01TT_msea_ora_result.csv", header = T)

plotdata <- tt.data[c(1:10),c(1,8,9)]

```
