#### Multi-omics landscape of NAT samples

library(dplyr)

prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1, rowNames = T)%>%as.data.frame()
kw_data <- as.data.frame(t(prodata[,c(1:39)]))

group <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 4, rowNames = T)%>%as.data.frame()

kw_data <- cbind(group[6], kw_data)

data1 <- kruskal.test(as.matrix(kw_data[2]) ~Segments, data = kw_data)

pval=c()

for (i in 2:ncol(kw_data)) {
  p = kruskal.test(as.matrix(kw_data[i]) ~Segments, data = kw_data)$p.value
  pval=c(pval,p)
  
}

table_out <- as.data.frame(cbind(row.names(prodata),pval))

write.csv(table_out, file = "../results/NAT-KW-test_in_proteome.csv")

## heatmap of NAT_DEP
library(ComplexHeatmap)
library(tibble)
prodata <- openxlsx::read.xlsx("./2022proteome/2022Proteome.xlsx", sheet = 1, rowNames = F)

NAT_DEP_list <- openxlsx::read.xlsx("./2022proteome/NAT-KW-test_in_proteome.xlsx", sheet = 2, rowNames = F)

NAT_DEP_exp <- merge(NAT_DEP_list[1], prodata[,1:40], by="Name")
NAT_DEP_exp <- column_to_rownames(NAT_DEP_exp, var = "Name")

plotdata <- as.data.frame(t(scale(t(log2(NAT_DEP_exp+1)))))

groups <- openxlsx::read.xlsx("./2022proteome/Clinicaldata2.xlsx", sheet = 4)

ht = HeatmapAnnotation(Segments = groups$Segments,
                       show_legend = T,
                       col = list(Segments = c("S1" = "#48466d",
                                               "S2"   = "#3d84a8",
                                               "S3"   = "#46cdcf",
                                               "S4"   = "#abedd8",
                                               "S5"   = "#f9ed69",
                                               "S6"   = "#f08a5d",
                                               "S7"   = "#b83b5e",
                                               "S8"   = "#6a2c70")))

col_fun = circlize::colorRamp2(c(-4, -2, 0, 2, 4),
                               c("#2c83b7","deepskyblue1", "white",
                                 "red", "#c2395a"))

p <- Heatmap(plotdata, 
             col = col_fun,
             top_annotation = ht, 
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = T,
             show_column_names = F,
             cluster_columns = T,
             clustering_method_columns = "ward.D",
             clustering_method_rows = "ward.D",
             cluster_rows = T,
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscore"))

## corrplot of NAT_DEP

res <- rcorr(as.matrix(log2(1+NAT_DEP_exp)))

col1 <- colorRampPalette(c( "black","#00007f","blue","#007fff","cyan",
                            "#7fff7f","#fff5a5","#ffaa64","#ff8264","#ff6464","#930077")) 

corrplot(res$r, order = "hclust",col=col1(300), col.lim=c(0.8,1),
         hclust.method = "ward.D",method = "circle",
         tl.col = "black", tl.pos = "lt", addrect = 2)


## NAT_LDA
library(MASS)
library(ggplot2)
library(MOVICS)
prodata <- openxlsx::read.xlsx("./2022proteome/2022Proteome.xlsx", sheet = 1, rowNames = T)
pt.elite <- getElites(dat=prodata[,1:39],
                      method="mad",
                      na.action="rm",
                      elite.pct = 0.3)
pt.data <- as.data.frame(scale(t(log2(pt.elite$elite.dat+1))))

lda.data <- as.data.frame(t(na.omit(t(pt.data))))

groups <- as.factor(c(rep("S1",3),rep("S2",5),rep("S3",5),
                      rep("S4",5),rep("S5",5),rep("S6",5),
                      rep("S7",5),rep("S8",6)))

lda.pt = lda(lda.data, groups)

result = predict(lda.pt, lda.data)

table(groups, result$class)

P = lda.pt$scaling 
means = lda.pt$means %*% P

total_means = as.vector(lda.pt$prior %*% means)
n_samples = nrow(lda.data)

x<-as.matrix(lda.data) %*% P - (rep(1, n_samples) %o% total_means)

subtybes <- as.data.frame(c(rep("Subtybe1",13),rep("Subtybe2",26)))
names(groups) <- "Groups"
names(subtybes) <- "Subtybes"
lda_data <- cbind(x, groups, subtybes)

p <- ggplot(data = lda_data, aes(x = LD1, y = LD2)) +
  geom_point(aes(color = groups, shape=groups), size = 2) +
  scale_shape_manual(values = c(0:2,5,6,3,4,8))+
  scale_color_manual(values = c("#48466d","#3d84a8","#46cdcf","#abedd8",
                                "#f9ed69","#f08a5d","#b83b5e","#6a2c70")) + 
  theme_bw() +
  theme(legend.position = "",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(x =  paste('LD1'), y = paste('LD2'), color = '')

g <- p + stat_ellipse(aes(fill = Subtybes), geom = 'polygon', level = 0.95, 
                      alpha = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = c("#0f4392","#ff5151"))

### NAT_Lipds_OPLS
library(MASS)
library(ggplot2)
library(ropls)
lipid.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 1, rowNames = T)

nat.data <- as.data.frame(log2(lipid.data[,c(2:40)]+1))

pls.data <- as.data.frame(t(nat.data))

groups <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 4, rowNames = T)

pdf("NAT_Lipid.oplsda.pdf")
lipid.oplsda <- opls(as.matrix(pls.data), groups$Subtypes2, predI = 1, orthoI = NA)
dev.off()

plot(lipid.oplsda,
     typeVc = "x-score",
     parAsColFcVn = groups$Subtypes2,
     parPaletteVc = c("#00b8a9","#f6416c"))

## NAT_Metabolites_OPLS
library(MASS)
library(ggplot2)
library(ropls)
lipid.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites2.xlsx", sheet = 1, rowNames = T)

pt.data <- as.data.frame(log2(lipid.data[,c(1:39)]+1))

pls.data <- as.data.frame(t(pt.data))

groups <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 4, rowNames = T)

pdf("NAT_Metabolites.oplsda2.pdf")
lipid.oplsda <- opls(as.matrix(pls.data), groups$Subtypes2, predI = 1, orthoI = NA)
dev.off()

plot(lipid.oplsda,
     typeVc = "x-score",
     parAsColFcVn = groups$Subtypes2,
     parPaletteVc = c("#00b8a9","#f6416c"))

## DEP in NAT Clusters

library(dplyr)
library(limma)
prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1,rowNames = T)

groups <- c(rep("Cluster1",13), rep("Cluster2",26))

design<- model.matrix(~0+factor(groups))
colnames(design)<- levels(factor(groups))

v <- voom(prodata[,1:39], design, plot=TRUE)

fit <- lmFit(v, design)

cont.matrix <- makeContrasts('Cluster1-Cluster2', levels = design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
write.csv(output, file = "NAT_DEPs_limmavoom.csv",row.names = T)

##DAL in NAT Clusters

library(dplyr)
library(limma)
lipid.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 1,rowNames = T)%>% as.data.frame()
nat.data <- as.data.frame(lipid.data[,c(2:40)]+1)

groups <- c(rep("Cluster1",18), rep("Cluster2",21))

design<- model.matrix(~0+factor(groups))
colnames(design)<- levels(factor(groups))
v <- voom(nat.data, design, plot=TRUE)

fit <- lmFit(v, design)

cont.matrix <- makeContrasts('Cluster1-Cluster2', levels = design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
write.csv(output, file = "NAT_DALs_limmavoom.csv",row.names = T)

##DAM in NAT Clusters

metabolites.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites2.xlsx", sheet = 2,rowNames = T)

groups <- c(rep("Cluster1",18), rep("Cluster2",21))

design<- model.matrix(~0+factor(groups))
colnames(design)<- levels(factor(groups))
v <- voom(metabolites.data[,1:39], design, plot=TRUE)

fit <- lmFit(v, design)

cont.matrix <- makeContrasts('Cluster1-Cluster2', levels = design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

output <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
write.csv(output, file = "NAT_DAMs_limmavoom.csv",row.names = T)

## Heatmap_NAT_hCluster2

library(tibble)
library(dplyr)
library(ComplexHeatmap)

PT_list <- openxlsx::read.xlsx("../data/NAT_DEPs_limmavoom2.xlsx", sheet = 2, rowNames = F)
PT_list <- as.data.frame(PT_list[,1])
names(PT_list) <- "Name" 

PT_list <- PT_list%>% distinct(Name, .keep_all = T)

prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1, rowNames = F)
## Read into KW test_expression matrix
dep.data <- merge(PT_list, prodata[,1:40], by = "Name")
dep.data <- column_to_rownames(dep.data, var = "Name")

plotdata <- as.data.frame(t(scale(log2(0.5+t(dep.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 4, rowNames = T)

ht = HeatmapAnnotation(Segments = clinicaldata$Segments,
                       show_legend = T,annotation_name_side = "left",
                       col = list(Segments = c("S1" = "#48466d",
                                               "S2"   = "#3d84a8",
                                               "S3"   = "#46cdcf",
                                               "S4"   = "#abedd8",
                                               "S5"   = "#f9ed69",
                                               "S6"   = "#f08a5d",
                                               "S7"   = "#b83b5e",
                                               "S8"   = "#6a2c70")))
ht2 = HeatmapAnnotation(Clusters = clinicaldata$Subtypes2,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Clusters = c("Cluster1" = "#0f4392",
                                                "" = "#ff5151"))) 
col_fun = circlize::colorRamp2(c(-4,-2, 0,2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(123)
p <- Heatmap(plotdata, column_split=clinicaldata$Subtypes2,
             col = col_fun,
             top_annotation = c(ht2,ht),
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,clustering_distance_columns = "euclidean",
             cluster_columns = T, clustering_method_columns = "ward.D",
             cluster_rows = T, 
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscores"))


cls <- plyr::ldply (row_order(p), data.frame)
names(cls) <- c("id", "rowid") 
cls <- dplyr::arrange(cls, rowid)
plotdata$cluster <- cls$id
write.csv(plotdata, file = "Heatmap_PT_DEP_cluster.csv")

### Heatmap of DAMs and DALs in NAT

library(tibble)
library(dplyr)
library(ComplexHeatmap)

dam_list <- openxlsx::read.xlsx("../data/NAT_DAMs_limmavoom2.xlsx", sheet = 2, rowNames = F)
dam_list <- as.data.frame(dam_list[,1])
names(dam_list) <- "Name" 

dam_list <- dam_list%>% distinct(Name, .keep_all = T)

metabolites.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites2.xlsx", sheet = 1, rowNames = F)

dep.data <- merge(dam_list, metabolites.data[,1:40], by = "Name")
dep.data <- column_to_rownames(dep.data, var = "Name")

plotdata <- as.data.frame(t(scale(log2(t(dep.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 4, rowNames = T)

ht = HeatmapAnnotation(Segments = clinicaldata$Segments,
                       show_legend = T,annotation_name_side = "left",
                       col = list(Segments = c("S1" = "#48466d",
                                               "S2"   = "#3d84a8",
                                               "S3"   = "#46cdcf",
                                               "S4"   = "#abedd8",
                                               "S5"   = "#f9ed69",
                                               "S6"   = "#f08a5d",
                                               "S7"   = "#b83b5e",
                                               "S8"   = "#6a2c70")))
ht2 = HeatmapAnnotation(Clusters = clinicaldata$Subtypes2,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Clusters = c("Cluster1" = "#0f4392",
                                                "Cluster2" = "#ff5151"))) 
col_fun = circlize::colorRamp2(c(-4,-2, 0,2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(123)
p <- Heatmap(plotdata, 
             col = col_fun,
             top_annotation = c(ht2,ht),
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F, 
             cluster_rows = T, 
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscores"))


dal_list <- openxlsx::read.xlsx("../data/NAT/NAT_DALs_limmavoom.xlsx", sheet = 2, rowNames = F)
dal_list <- as.data.frame(dal_list[,1])
names(dal_list) <- "Name" 

dal_list <- dal_list%>% distinct(Name, .keep_all = T)

lipids.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 1, rowNames = F)

dep.data <- merge(dal_list, lipids.data[,c(1:41)], by = "Name")
dep.data <- column_to_rownames(dep.data[,-2], var = "Name")

plotdata <- as.data.frame(t(scale(log2(t(dep.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 4, rowNames = T)

ht = HeatmapAnnotation(Segments = clinicaldata$Segments,
                       show_legend = T,annotation_name_side = "left",
                       col = list(Segments = c("S1" = "#48466d",
                                               "S2"   = "#3d84a8",
                                               "S3"   = "#46cdcf",
                                               "S4"   = "#abedd8",
                                               "S5"   = "#f9ed69",
                                               "S6"   = "#f08a5d",
                                               "S7"   = "#b83b5e",
                                               "S8"   = "#6a2c70")))
ht2 = HeatmapAnnotation(Clusters = clinicaldata$Subtypes2,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Clusters = c("Cluster1" = "#0f4392",
                                                "Cluster2" = "#ff5151"))) 
col_fun = circlize::colorRamp2(c(-4,-2, 0,2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(123)
p <- Heatmap(plotdata, 
             col = col_fun,
             top_annotation = c(ht2,ht),
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,clustering_distance_columns = "euclidean",
             cluster_columns = F, clustering_method_columns = "ward.D",
             cluster_rows = T, 
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscores"))

## Bubble plot of MSEA in NAT-cluster

library(ggplot2)

ll.data <- read.csv("../data/2022ProteomeMetabolism/NAT-Cluster1-msea_ora_result.csv", header = T)

plotdata <- ll.data[c(1:8),c(1,8,9)]

p <- ggplot(plotdata,aes(x = EnrichmentRatio, 
                         y = reorder(Name,EnrichmentRatio,sum), 
                         size = EnrichmentRatio,
                         colour=LogP)) +
  geom_point(shape = 16) +                
  labs(x = "EnrichmentRatio", y = "")+          
  scale_colour_gradient2(name="-log10(P-value)",
                         low = "#ffbd39",
                         mid = "#e61c5d",
                         high = "#930077",
                         midpoint = 3)+
  theme_bw()          


RL.data <- read.csv("../data/2022ProteomeMetabolism/NAT-cluster2-msea_ora_result.csv", header = T)

plotdata <- RL.data[c(1:8),c(1,8,9)]

## LineChart_NAT_DAL

library(ggplot2)
library(reshape2)
library(tibble)
library(dplyr)
hex.data <- openxlsx::read.xlsx("../data/NAT/NAT_lipid_mean2023.xlsx", sheet = 9, rowNames = T)

hex.data <- as.data.frame(t(scale(t(hex.data))))
hex.data <- rownames_to_column(hex.data, var = "Name")

plotdata <- melt(hex.data, by="Name")

mean <- plotdata %>% 
  group_by(variable) %>% 
  summarise(value = mean(value))  
mean$variable <- as.numeric(mean$variable)

p <- ggplot(plotdata, aes(x=variable, y=value)) +
  geom_line(aes(group=Name), color="#e84545", alpha=0.03, size=2)+
  geom_line(data = mean, color="#e84545", alpha = 0.7, size=2) + 
  geom_point(data = mean)+
  xlab("")+ylab("Relative Abundance")+
  theme_bw()


## NAT_annotation

library(ggplot2)
cluster1 <- openxlsx::read.xlsx("../data/NAT/NAT_hcluster_FINAL_GO.xlsx", sheet = 2)

c1 <- ggplot(cluster1, aes(x = reorder(factor(Description), `-LogP`), 
                           y = `-LogP`)) + 
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

cluster2 <- openxlsx::read.xlsx("../data/NAT/NAT_hcluster_FINAL_GO.xlsx", sheet = 3)
c2 <- ggplot(cluster2, aes(x = reorder(factor(Description), `-LogP`), 
                           y = `-LogP`)) + 
  geom_bar(stat = 'identity', fill = '#ff5151', width = 0.7)+
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


