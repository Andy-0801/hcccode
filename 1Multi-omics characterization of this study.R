## CPTAC data
library(dplyr)
prodata <- data.table::fread("./S049_Liver_Cancer_Gao2019/Liver_Cancer_Proteome_CDAP_Protein_Report.r1/Zhou_Liver_Cancer_Proteome.tmt11.tsv", header = T)%>% as.data.frame()
rt <- prodata[!duplicated(prodata$Gene),]
rownames(rt)<-rt[,1]
rt<-rt[,-1]
rt=as.matrix(rt)

selectCol=seq(2,ncol(rt),2) 
rt=rt[,selectCol]

hcc_data <- as.data.frame(rt[-(1:3),])
openxlsx::write.xlsx(hcc_data,file="CPTAC_proteome_data.xlsx", rowNames = T, colName = T)

## Venn graph
library(VennDiagram)
library(grDevices)
library(dplyr)
filedata <- openxlsx::read.xlsx("./2022MBR/Venn/Venngene_list.xlsx", sheet = 1)%>%as.data.frame()
p <- venn.diagram(
  list(CPTAC=filedata$CPTAC, This_Study=filedata$This.Study), 
  na = "remove", filename = NULL,  col = "black", lty = 1, lwd = 3,
  fill = c("cornflowerblue", "darkorchid1"), alpha = 0.50, cex = 2.0,
  fontfamily = "serif", fontface = "bold",
  cat.col = c("cornflowerblue", "darkorchid4"), cat.cex = 1.8,
  cat.fontface = "bold",  cat.fontfamily = "serif")

pdf(file="CPTACvsThis.study_Venn.pdf")
grid.draw(p)
dev.off()

pt <- venn.diagram(
  list(Tumor=filedata$Tumor, Peritumor=filedata$Peritumor), 
  na = "remove", filename = NULL,  lty = 1, lwd = 3,
  col = c('#c2395a', '#2c83b7'), alpha = 0.50, cex = 2.0,
  fontfamily = "serif", fontface = "bold",
  cat.col = c('#c2395a', '#2c83b7'), cat.cex = 1.8,
  cat.fontface = "bold",  cat.fontfamily = "serif")

pdf(file="2NATvsTT_Venn.pdf")
grid.draw(pt)
dev.off()

## Lollipop chart
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(ggpubr)
no.data <- openxlsx::read.xlsx("./2022MBR/Proteome_MBR.xlsx", sheet = 3)
lollipop <- no.data[,c(4:6)]

g <- ggplot(lollipop) +
  geom_segment(aes(x=SampleID, xend=SampleID, y=Peritumor, yend=Tumor),
               color="grey") +
  geom_point(aes(x=SampleID, y=Peritumor), color="#2c83b7", size=2 ) +
  geom_point(aes(x=SampleID, y=Tumor), color="#c2395a", size=2 ) +
  theme(legend.position = "none") +
  theme_classic() + ylim(5000,8000)+
  xlab("Paired tumor and peritumor samples (n = 39)") +
  ylab("Number of protein identifications")

## Boxplot：the counts of identified protein in TT and NAT
plotdata <- no.data[,c(1:3)]

b <- ggplot(plotdata, aes(x = Groups, y = No.)) + 
  geom_jitter(aes(color = Groups), size=1, width = 0.3)+
  geom_boxplot(aes(color = Groups), size=1) + 
  scale_color_manual(values=c(Peritumor = "#2c83b7", Tumor = "#c2395a"))+
  stat_compare_means()+
  theme_classic() + ylim(5000,8000)+
  theme(legend.position = "none") +
  xlab("") +
  ylab("Number of protein identifications")

## Boxplot：the relative expression of Liver sepecific protein in TT and NAT
library(tibble)
library(reshape2)
library(ggplot2)
library(ggpubr)
prodata <- openxlsx::read.xlsx("./2022proteome/2022Proteome.xlsx", sheet = 1,rowNames = F)

proteinlist <- openxlsx::read.xlsx("./2022proteome/Liver_specific_genes.xlsx", sheet = 1,rowNames = F)
proteinlist <- as.data.frame(proteinlist[,1])
names(proteinlist) <- "Name"
specific_protein_exp <- merge(proteinlist, prodata, by="Name")
specific_protein_exp <- column_to_rownames(specific_protein_exp, var = "Name")

##0.00001 replaced by NA
plotdata <- as.data.frame(t(log2(specific_protein_exp+1)))
groups <- as.data.frame(c(rep("NAT",39), rep("TT",39))) 
names(groups) <- "Groups"

plotdata <- cbind(groups, plotdata)

data <- melt(plotdata, id.vars = "Groups")

p <- ggplot(data, aes(x=Groups, y=value, color=Groups))+geom_boxplot()+
  scale_color_manual(values=c(NAT = "#2c83b7", TT = "#c2395a"))+
  stat_compare_means()+
  theme_classic()+ 
  ylab("Log2(Normalized expression value)")+
  xlab("Liver specific proteins")+
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


## Metabolites Count Pieplot
library(ggplot2)
library(RColorBrewer)
metabolite_data <- openxlsx::read.xlsx("./2022ProteomeMetabolism/Metabolites2023.xlsx", sheet = 1, rowNames = F)

count <- as.data.frame(table(metabolite_data$Super_Class))
write.csv(count, file = "Metablites_class.csv", row.names = F)

count <- read.csv("./2022ProteomeMetabolism/Metablites_class.csv")
count <- count[order(count$Var1, decreasing = TRUE),]

p <- ggplot(data = count, mapping = aes(x = '', y = Freq, fill = Var1)) +
  geom_bar(stat = 'identity', position = 'stack', width = 1)+
  coord_polar(theta = 'y')

label_value <- paste('(', round(count$Freq/sum(count$Freq) * 100, 2), '%)', 
                     sep = '')
label_value

label <- paste(count$Var1, label_value, sep = '')

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=5, face="bold")
  )
g <- p + scale_fill_brewer(palette = "Paired")+
  labs(x = '', y = '', title = '') + blank_theme +
  theme(legend.position = "none") + 
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = count$Freq/5 +
                  c(0,cumsum(count$Freq)[-length(count$Freq)]),
                x = sum(count$Freq)/600, label = label)) 


## Venn graph of metabolites
library(VennDiagram)
library(grDevices)
library(dplyr)
filedata <- openxlsx::read.xlsx("./2022ProteomeMetabolism/Metabolites2023.xlsx", sheet = 4)%>%as.data.frame()
p <- venn.diagram(
  list(Negative=filedata$Negative, Positive=filedata$Positive), 
  na = "remove", filename = NULL,  lty = 1, lwd = 3,
  col = c('#c2395a', '#2c83b7'), alpha = 0.50, cex = 2.0,
  fontfamily = "serif", fontface = "bold",
  cat.col = c('#c2395a', '#2c83b7'), cat.cex = 1.8,
  cat.fontface = "bold",  cat.fontfamily = "serif")

pdf(file="2quadruple_metabolites.pdf")
grid.draw(p)
dev.off()

##Lipid Count Pieplot
library(ggplot2)
library(RColorBrewer)
lipid_data <- openxlsx::read.xlsx("./2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 1, rowNames = F)

count <- as.data.frame(table(lipid_data$Category))
write.csv(count, file = "Lipid_class.csv", row.names = F)

count <- read.csv("./2022ProteomeMetabolism/Lipid_class.csv")
count <- count[order(count$Var1, decreasing = TRUE),]

p <- ggplot(data = count, mapping = aes(x = '', y = Freq, fill = Var1)) +
  geom_bar(stat = 'identity', position = 'stack', width = 1)+
  coord_polar(theta = 'y')

label_value <- paste('(', round(count$Freq/sum(count$Freq) * 100, 2), '%)', 
                     sep = '')
label_value

label <- paste(count$Var1, label_value, sep = '')
##bank theme
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=5, face="bold")
  )
g <- p + scale_fill_brewer(palette = "Set1")+
  labs(x = '', y = '', title = '') + blank_theme +
  theme(legend.position = "none") + 
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = count$Freq/10 +
                  c(0,cumsum(count$Freq)[-length(count$Freq)]),
                x = sum(count$Freq)/2000, label = label)) 

## Venn graph of lipids
library(VennDiagram)
library(grDevices)
library(dplyr)
filedata <- openxlsx::read.xlsx("./2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 4)%>%as.data.frame()
p <- venn.diagram(
  list(Negative=filedata$Negative, Positive=filedata$Positive), 
  na = "remove", filename = NULL,  lty = 1, lwd = 3,
  col = c('#c2395a', '#2c83b7'), alpha = 0.50, cex = 2.0,
  fontfamily = "serif", fontface = "bold",
  cat.col = c('#c2395a', '#2c83b7'), cat.cex = 1.8,
  cat.fontface = "bold",  cat.fontfamily = "serif")

pdf(file="2quadruple_lipids.pdf")
grid.draw(p)
dev.off()

```
