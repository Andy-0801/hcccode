##OPLS-DA in Metabolome

library(ropls)
Metabolite_data <- openxlsx::read.xlsx("./2022ProteomeMetabolism/Metabolites2023.xlsx", sheet = 1,rowNames = T)

group <- openxlsx::read.xlsx("./2022proteome/Clinicaldata2.xlsx", sheet = 2,rowNames = F)

# OPLS-DA
pdf(../results/"Metabolite_oplsda.pdf")
Metabolite.oplsda <- opls(t(log2(Metabolite_data[,-c(1:12)])), group$Groups, predI = 1, orthoI = NA)
dev.off()

plot(Metabolite.oplsda,
     typeVc = "x-score",
     parAsColFcVn = group$Groups,
     parPaletteVc = c("#2c83b7", "#c2395a"))

#Plot in ggplot2
library(ggplot2)
library(ggsci)
library(tidyverse)

sample.score = Metabolite.oplsda@scoreMN %>%  
  as.data.frame() %>%
  mutate(group = group$Groups,
         o1 = Metabolite.oplsda@orthoScoreMN[,1]) 
head(sample.score)

p <- ggplot(sample.score, aes(p1, o1, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  #geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'OPLS-DA:t1(6.0%)',y = 'OPLS-DA:to1') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) + 
  scale_color_manual(values = c('#2c83b7','#c2395a')) +
  theme_bw() +
  theme(legend.position = c(0.92,0.92),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

Metabolite.vip <- as.data.frame(getVipVn(Metabolite.oplsda))
write.csv(lipid.vip, file = "../results/Metabolite_VIP.csv", row.names = T)

## OPLS-DA in Lipidome

Lipid_data <- openxlsx::read.xlsx("./2022ProteomeMetabolism/Lipid2023.xlsx", sheet = 1,rowNames = T)

group <- openxlsx::read.xlsx("./2022proteome/Clinicaldata2.xlsx", sheet = 2,rowNames = F)

## OPLS-DA
pdf("Lipid_data_oplsda.pdf")
Lipid.oplsda <- opls(t(log2(Lipid_data[,-1])), group$Groups, predI = 1, orthoI = NA)
dev.off()

plot(Lipid.oplsda,
     typeVc = "x-score",
     parAsColFcVn = group$Groups,
     parPaletteVc = c("#2c83b7", "#c2395a"))
#Plot in ggplot2 
sample.score = Lipid.oplsda@scoreMN %>% 
  as.data.frame() %>%
  mutate(group = group$Groups,
         o1 = Lipid.oplsda@orthoScoreMN[,1]) 
head(sample.score)

p <- ggplot(sample.score, aes(p1, o1, color = group)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) + 
  geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
  geom_point() +
  #geom_point(aes(-10,-10), color = 'white') +
  labs(x = 'OPLS-DA:t1(6.0%)',y = 'OPLS-DA:to1') +
  stat_ellipse(level = 0.95, linetype = 'solid', 
               size = 1, show.legend = FALSE) + 
  scale_color_manual(values = c('#2c83b7','#c2395a')) +
  theme_bw() +
  theme(legend.position = c(0.92,0.92),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

lipid.vip <- as.data.frame(getVipVn(Lipid.oplsda))
write.csv(lipid.vip, file = "../results/Lipid_VIP.csv", row.names = T)


