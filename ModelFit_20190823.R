setwd("~/Documents/ANALYSES/TVB_post/results_TVBii_VS")
ModelFit <- read.csv("ModelFit_VS.csv")
str(ModelFit)
ModelFit$Condition <- factor(ModelFit$Condition, ordered = T, 
                             levels = c('PSE_Fcsim-Fcemp','PSE_SC-Fcemp','VS_Fcsim-Fcemp','VS_SC-Fcemp'),
                             labels = c("PSE: FCsim", "PSE: SC", "VS: FCsim", "VS: SC"))

library(ggplot2)
png('ModelFit2g.png', width=1000, height=600)
ggplot(ModelFit, aes(x=subID, y=cor)) + 
  geom_point(aes(fill=Condition,shape=Condition), size=10) +
  scale_shape_manual(values=c(23,21,24,22)) +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 30),
        legend.key.size=unit(3,"line"), 
        legend.title=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(size=6)), size=FALSE) +
  labs(x = "", y = "Prediction accuracy") +
  scale_fill_brewer(palette="Paired")
dev.off()

#tryout with barchart
ggplot(ModelFit, aes(x=subID, y=cor, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge2()) + 
  labs(x = "", y = "Model fit (correlation with FCemp)") +
  scale_fill_brewer(palette="Paired",  labels = c("PSE: SC", "PSE: FCsim", "VS: SC", "VS: FCsim"))



