#### Ribosome biogenesis

library(dplyr)
library(tibble)
prodata <- openxlsx::read.xlsx("../data/2022proteome/2022Proteome.xlsx", sheet = 1, rowNames = F)%>%as.data.frame()
rb.list <- openxlsx::read.xlsx("../data/NMF_Ribosome b_genes.xlsx", sheet = 6, rowNames = F, startRow = 5)
rb.list <- as.data.frame(rb.list[,1])
names(rb.list) <- "Name" 

rb.list <- rb.list%>% distinct(Name, .keep_all = T)

rb.data <- merge(rb.list, prodata, by="Name",sort = F)
write.csv(rb.data, file = "Ribosome b_genes.csv", row.names = F)

rb.data <- column_to_rownames(rb.data, var = "Name")
rb.data <- as.data.frame(t(rb.data))

clinicaldata <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3, rowNames = T)%>%as.data.frame()

rb.aov <- cbind(clinicaldata$Subtypes, rb.data)
names(rb.aov)[1] <- "Subtypes"
rb.aov <- rb.aov[-c(1:39),]

anova1 <- summary(aov(rb.aov[,2]~Subtypes, rb.aov))

output <- as.data.frame(anova1[[1]])


for (i in 2:ncol(rb.aov)) {
  anova = summary(aov(rb.aov[,i]~Subtypes, rb.aov))
  output1 <- as.data.frame(anova[[1]])
  output <- rbind(output1,output)
}

write.csv(output, file="../results/Annova.results.csv")

#heatmap

library(dplyr)
library(ComplexHeatmap)
rb.mean <- openxlsx::read.xlsx("../data/NMF_Ribosome b_genes.xlsx", sheet = 7, rowNames = T)

plotdata <- as.data.frame(scale(t(rb.mean)))

col_fun = circlize::colorRamp2(c(-1.5,-0.5, 0, 0.5, 1.5),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(plotdata[i, j], 1), 
    x, y,
    gp = gpar(
      fontsize = 10
    ))
}

set.seed(123)
p <- Heatmap(plotdata,
             rect_gp = gpar(col = "white", lwd = 2), 
             col = col_fun, cell_fun = cell_fun,
             show_row_names = T,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = T,
             cluster_columns = F, 
             cluster_rows = F, 
             row_names_gp = gpar(fontsize = 10),
             heatmap_legend_param = list(title = "Zscores"))

# Ribosome proteins

library(dplyr)
library(ComplexHeatmap)
RP.data <- openxlsx::read.xlsx("../data/NMF_Ribosome b_genes.xlsx", sheet = 3, rowNames = T)%>%as.data.frame()

plotdata <- as.data.frame(t(scale(log2(1+t(RP.data)))))

clinicaldata <- openxlsx::read.xlsx("../data/2022MBR/Clinicaldata.xlsx", sheet = 3, rowNames = T)%>%as.data.frame()

sub = HeatmapAnnotation(Subtypes = clinicaldata$Subtypes,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Subtypes = c("NAT"= "grey",
                                                "C1" = "#2874C5",
                                                "C2" = "#EABF00",
                                                "C3" = "#C6524A")))
seg = HeatmapAnnotation(Segments = clinicaldata$Segments,
                        show_legend = T,annotation_name_side = "left",
                        col = list(Segments  = c("S1" = "#48466d", "S2"   = "#3d84a8",
                                                 "S3" = "#46cdcf", "S4"   = "#abedd8",
                                                 "S5" = "#f9ed69", "S6"   = "#f08a5d",
                                                 "S7" = "#b83b5e", "S8"   = "#6a2c70")))
hlob = HeatmapAnnotation(Hepatic_lobes = clinicaldata$Hepatic_lobes,
                         show_legend = T,annotation_name_side = "left",
                         col = list(Hepatic_lobes=c("Left"= "#0f4392","Right"= "#ff5151")))

group = HeatmapAnnotation(Groups = clinicaldata$Groups,
                          show_legend = T,annotation_name_side = "left",
                          col = list(Groups=c("NAT"= "#0f4392","TT"= "#ff5151")))

mark_gene <- openxlsx::read.xlsx("../data/NMF_Ribosome b_genes.xlsx", sheet = 4)%>%
  as.data.frame()

row_anno <-  rowAnnotation(link = anno_mark(at = which(rownames(plotdata) 
                                                       %in% mark_gene$Name), 
                                            labels = mark_gene$Name,
                                            labels_gp = gpar(fontsize = 7.5)))


col_fun = circlize::colorRamp2(c(-4,-2, 0, 2, 4),
                               c("#49a09d","#00AFBB","white",
                                 "#FC4E07", "#b20a2c"))
set.seed(111)
p <- Heatmap(plotdata,
             col = col_fun, column_split = clinicaldata$Subtypes,
             top_annotation = c(sub,seg,hlob),
             show_row_names = F,
             show_row_dend = F,
             show_column_dend = F,
             show_column_names = F,
             cluster_columns = T, 
             cluster_rows = F, 
             row_names_gp = gpar(fontsize = 7.5),
             heatmap_legend_param = list(title = "Zscores"))
p+row_anno

# Ribosome biogenesis and ribosome proteins

library(ggplot2)
library(ggrepel)
output <- openxlsx::read.xlsx("../data/Proteome_DEGs_limmavoom.xlsx", sheet = 3, rowNames = T)
names(output) <- c('Groups','log2FC', 'AveExpr', 'Pvalue', 'BH-Pvalue')

yanse <- ifelse(output$`BH-Pvalue`<0.01 & abs(output$log2FC)>= 1, ifelse(output$log2FC> 1,'red','grey'),'green')

p <- ggplot(output, aes(output$log2FC, -log10(output$`BH-Pvalue`), col=Groups))+geom_point(size=1)

g <- p+theme_bw()+ 
  scale_color_manual(values = c('#a82ffc','#ef7e56'))+ 
  labs(x="log2 (Fold Change)",y="-log10 (BH-Pvalue)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

label <- ifelse(-log10(output$`BH-Pvalue`) > 5 & abs(output$log2FC) >= 1, rownames(output),"")

l <- g+geom_text_repel(data=output, aes(x = output$log2FC, 
                                        y = -log10(output$`BH-Pvalue`), 
                                        label = label),
                       size = 3,box.padding = unit(0.2, "lines"),
                       segment.color = "black", 
                       show.legend = FALSE)




