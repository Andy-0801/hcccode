##Proteome_PCA (R 3.6.3)
library(FactoMineR)
library(ggplot2)
library(dplyr)
prodata <- openxlsx::read.xlsx("./2022proteome/2022Proteome.xlsx", sheet = 1,rowNames = T)

plotdata <- as.data.frame(t(log2(prodata)))

prodata.pca <- PCA(plotdata, ncp = 2, scale.unit = TRUE, graph = FALSE)

plot(prodata.pca)

pca_data <- data.frame(prodata.pca$ind$coord[,1:2])

pca_eig1 <- round(prodata.pca$eig[1,2], 2)
pca_eig2 <- round(prodata.pca$eig[2,2], 2)

groups <- as.data.frame(c(rep("NAT",39), rep("TT",39)))
names(groups) <- "Groups"

pca_data <- cbind(pca_data, groups)
pca_data$SampleID <- rownames(pca_data)

p <- ggplot(data = pca_data, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Groups), size = 1) +
  scale_color_manual(values = c('#2c83b7', '#c2395a')) + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  labs(x =  paste('PCA1:', pca_eig1, '%'), 
       y = paste('PCA2:', pca_eig2, '%'), color = '')

g <- p + stat_ellipse(aes(fill = Groups), geom = 'polygon', level = 0.95, 
                      alpha = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = c('#2c83b7', '#c2395a'))

##Corplot_QC for MultiOmics data
library(Hmisc)
library(corrplot)
QC_data <- read.csv("../data/2022ProteomeMetabolism/QC for lipidome.csv", header = T, row.names = 1)

plotdata <- as.data.frame(log2(QC_data+1))

res <- rcorr(as.matrix(plotdata))

corrplot(res$r, order = "hclust", 
         hclust.method = "ward.D2",method = "circle",
         tl.col = "black", type = "upper", tl.pos = "d")

corrplot(res$r, add = T,order = "hclust", 
         hclust.method = "ward.D2",
         method = "number",type = "lower", tl.col = "black",diag=FALSE,
         tl.pos = "n", cl.pos="n")              


