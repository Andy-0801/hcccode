#### HCC molecular subtypes
## proteome NMF cluster

library(NMF)
library(doMPI)
library(dplyr)
prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 2, rowNames = T)%>%as.data.frame()

prodata <- log2(1+prodata)

mads <- apply(prodata,1,mad) 
pro.exp <- prodata[rev(order(mads))[1:4000],] 

res <- nmf(pro.exp,2:8,nrun=10,seed=123,method = "lee")
plot(res)

res <- nmf(pro.exp,3,nrun=10,seed=123)

group <- predict(res)
table(group)
write.csv(as.data.frame(group), file = "Proteome_subtypes4000.csv")

coefmap(res)
consensusmap(res)

index <- extractFeatures(res,0.95) 
sig.order <- unlist(index)
NMF.Exp <- prodata[sig.order,]
NMF.Exp <- na.omit(NMF.Exp)

library(tinyarray)
dp = NMF.Exp[,order(group)]
draw_heatmap(dp,sort(group), split_column = 3,
             color_an = c("#2874C5","#EABF00","#C6524A","#868686"),
             annotation_legend = T,
             cluster_cols = F,
             show_rownames = T)


jco <- c("#2874C5","#EABF00","#C6524A")
res2 <- nmf(NMF.Exp,3,nrun=10,seed=123)

coefmap(res2)
consensusmap(res2)
group2 <- predict(res2)

basismap(res, distfun="euclidean", hclustfun="ward",
         cexCol = 1,
         cexRow = 0.3,
         annColors=list(c("1"=jco[1],"2"=jco[2],"3"=jco[3])))

group <- predict(res) 
table(group)
consensusmap(res,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(res)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3])))

## proteome 3 clusters

## Different expression analysis

library(limma)
library(dplyr)
##Proteome
prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 2, rowNames = T)%>%as.data.frame()
prodata <- as.data.frame(log2(1+prodata))

groups <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 3, rowNames = TRUE)
groups <- groups[1:39,]
design<- model.matrix(~0+factor(groups$Subtypes))
colnames(design)<- levels(factor(groups$Subtypes))

cont.wt <- makeContrasts("C1-C2", "C1-C3","C2-C3", "C2-C1",
                         levels=design) 

fit <- lmFit(prodata, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTableF(fit2, adjust="BH",sort.by="F",n=Inf)

write.csv(tT, file = "NMF_C1-2-3_DEPs.csv",row.names = T)

##metabolome 
metabo.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Metabolites2023.xlsx", sheet = 5, rowNames = T)%>%as.data.frame()
metabo.data <- as.data.frame(log2(1+metabo.data))

groups <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 3, rowNames = TRUE)
groups <- groups[1:39,]
design<- model.matrix(~0+factor(groups$Subtypes))
colnames(design)<- levels(factor(groups$Subtypes))

cont.wt <- makeContrasts("C1-C2", "C1-C3","C2-C3", "C2-C1",
                         levels=design) 

fit <- lmFit(metabo.data, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTableF(fit2, adjust="BH",sort.by="F",n=Inf)

write.csv(tT, file = "NMF_C1-2-3_DAMs.csv",row.names = T)

## lipidome
lipid.data <- openxlsx::read.xlsx("../data/2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 3, rowNames = T)%>%as.data.frame()
lipid.data <- as.data.frame(log2(1+lipid.data))

groups <- openxlsx::read.xlsx("../data/2022proteome/Clinicaldata2.xlsx", sheet = 3, rowNames = TRUE)
design<- model.matrix(~0+factor(groups$Subtypes))
colnames(design)<- levels(factor(groups$Subtypes))

cont.wt <- makeContrasts("C1-C2", "C1-C3","C2-C3", "C2-C1",
                         levels=design) 

fit <- lmFit(lipid.data, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTableF(fit2, adjust="BH",sort.by="F",n=Inf)

write.csv(tT, file = "../results/NMF_C1-2-3_DALs.csv",row.names = T)


## NMF clusters heatmap
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
nmf_list <- openxlsx::read.xlsx("../data/NMF_C1-2-3_DEPs.xlsx", sheet = 2, rowNames = F)
nmf_list <- as.data.frame(nmf_list[,1])
names(nmf_list) <- "Name" 

nmf_list <- nmf_list%>% distinct(Name, .keep_all = T)

prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 2, rowNames = F)%>%as.data.frame()

dep.data <- merge(nmf_list, prodata, by = "Name", sort = FALSE)
dep.data <- column_to_rownames(dep.data, var = "Name")

plotdata <- as.data.frame(t(scale(log2(1+t(dep.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 2, rowNames = T)%>%as.data.frame()

sub = HeatmapAnnotation(Subtypes = clinicaldata$Subtypes,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Subtypes = c("C1" = "#2874C5",
                                                "C2" = "#EABF00",
                                                "C3" = "#C6524A")))
seg = HeatmapAnnotation(Segments = clinicaldata$Segments,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Segments  = c("S1" = "#48466d", "S2"   = "#3d84a8",
                                                 "S3" = "#46cdcf", "S4"   = "#abedd8",
                                                 "S5" = "#f9ed69", "S6"   = "#f08a5d",
                                                 "S7" = "#b83b5e", "S8"   = "#6a2c70")))

col_fun = circlize::colorRamp2(c(-4,-2, 0, 2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(123)
p <- Heatmap(plotdata, km=3,
             col = col_fun, column_split = clinicaldata$Subtypes,
             top_annotation = c(sub, seg),
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = F, 
             cluster_rows = F, 
             row_names_gp = gpar(fontsize = 6),
             heatmap_legend_param = list(title = "Zscores"))


## NMF clusters anotations

library(pheatmap)
plotdata <- openxlsx::read.xlsx("../data/NMF3_FINAL_GO.xlsx", sheet=3, rowNames = T)
p <- pheatmap(plotdata, cluster_cols = FALSE, cluster_rows = FALSE,
              border_color = "black",
              color = colorRampPalette(c("white","red","firebrick3"))(20))

