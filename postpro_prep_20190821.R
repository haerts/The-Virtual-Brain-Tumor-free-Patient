###################################################################################################################
###                                                                                                             ###
###                   MODELING BRAIN DYNAMICS AFTER TUMOR RESECTION USING THE VIRTUAL BRAIN                     ###
###               ==============================================================================                ###
###                                       PART 1: PREPARATION & DESCRIPTIVES                                    ###
###                                                                                                             ###
### Created by Hannelore Aerts                                                                                  ###
### Date last update: 21/08/2019                                                                                ###
###################################################################################################################

# Analysis remarks:
# --> use Ji corrected for ROI size only, not ROI size + in-strength   
#     UPDATE: use uncorrected Ji!
#     UPDATE 2: also use corrected Ji to check robustness of results
# --> use J=JiT=JiNT for controls
#     08/04/2019: also use Jres=JiTres=JiNTres for controls
# --> use second batch of TVBii analyses with normalization factor=pre-op normalization factor
#     = 75K (instead of 45K in first batch)
# --> 21/03/2019: don't regress out motivation from cognitive performance scores + new plots with ggplot
# --> 11/04/2019: resnorm use !!pre-operative!! mean and sd 


library('PMCMR')
library('car')
library('viridis')
library('corrplot')
library('colorspace')
library('ggplot2')
library('gridExtra')


### Read in data and prepare for analyses ------------------------------------------------------------------------#

# Wide data format

setwd("/home/hannelore/Documents/ANALYSES/TVB_post")
results=read.table(file="RESULTS_ALL.csv", header=TRUE, sep=",")
str(results)

results=within(results, {
  subID=as.character(subID)
  date_t1 = as.Date(date_t1, "%d/%m/%y")
  date_t2 = as.Date(date_t2, "%d/%m/%y")
  group = factor(group, ordered=TRUE, levels=c('CON', 'MEN', 'GLI'))
  fmri_TR_t1 = as.factor(fmri_TR_t1)
  fmri_TR_t2 = as.factor(fmri_TR_t2)
})

# Select only subjects with post-op data, including MRI
idx=!is.na(results$intnorm_scaling_t2)
results=results[idx,]
rm(idx)

attach(results)


###################################################################################################################
###                                   PART 1: DEMOGRAPHICS DESCRIPTIVES                                         ###
###################################################################################################################

# 1) follow-up time 
results$date_diff = difftime(date_t2, date_t1, units="days") / 30
detach(results); attach(results)
summary(as.numeric(date_diff))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.200   6.583   8.167   7.906   8.783  10.733 
plot(as.numeric(date_diff)~group, ylab="follow-up time (months)", col="gray")
Anova(aov(as.numeric(date_diff)~group)) #F(2,24)=0.18, p=0.84


# 2) age per group (baseline age in post-op participants)
plot(age ~ group, col="gray")
aggregate(age, list(group), mean)
#  Group.1        x
#1     CON 59.60000
#2     MEN 57.90000
#3     GLI 50.71429
vars=aggregate(age, list(group), var)
sqrt(vars$x)
#[1] 10.31935 10.94887 11.65782
rm(vars)
Anova(aov(age~group)) #F(2,24)=1.47, p=0.25


# 3) sex per group (in post-op participants)
plot(sex ~ group)
table(sex,group)
#group
#   CON MEN GLI
#F   4   8   3
#M   6   2   4


# 4) handedness post-op + difference wrt pre-op per group
table(hand_t2, group)
#hand_t2 CON MEN GLI
#-0.9    1   0   0
#-0.79   0   1   0
#-0.6    0   0   1
#0.68    0   0   1
#0.85    1   0   0
#1       8   9   5
results$handdiff=hand_t2-hand_t1; detach(results); attach(results)
boxplot(handdiff~group, pch=21, bg="black", ylab="handedness difference (t2 - t1)", xlab="group")
Anova(aov(handdiff~group)) #F(2,24)=0.40, p=0.68


# 5) lesion volume post-op + difference wrt pre-op per group
plot(lesion_vol_cm3_t2 ~ group, col="gray", ylab="lesion volume (cm³)", pch=21, bg="black")
results$lesiondiff=lesion_vol_cm3_t2-lesion_vol_cm3; detach(results); attach(results)
plot(lesiondiff~group, col="gray", ylab="lesion volume difference (t2 - t1)", pch=21, bg="black")


# 6) intensity normalization scaling factor post-op + difference wrt pre-op per group
plot(intnorm_scaling_t2 ~ group, col="gray", ylab="intensity normalization scaling factor", pch=21, bg="black")
Anova(aov(intnorm_scaling_t2~group)) #F(2,24)=0.35, p=0.71
results$intnormdiff=intnorm_scaling_t2-intnorm_scaling_t1; detach(results); attach(results)
plot(intnormdiff~group, col="gray", ylab="intensity normalization scaling factor difference (t2 - t1)", 
     pch=21, bg="black")
Anova(aov(intnormdiff~group)) #F(2,24)=1.79, p=0.19


# 7) motion post-op + difference wrt pre-op per group
plot(fmri_MD_t2 ~ group, col="gray", ylab="fMRI mean displacement", pch=21, bg="black")
Anova(aov(fmri_MD_t2~group)) #F(2,24)=0.29, p=0.75
results$fmri_MDdiff=fmri_MD_t2-fmri_MD_t1; detach(results); attach(results)
plot(fmri_MDdiff~group, col="gray", ylab="fMRI mean displacement difference (t2 - t1)", pch=21, bg="black")



###################################################################################################################
###                                   PART 2: MODEL PARAMETERS DESCRIPTIVES                                     ###
###################################################################################################################

setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii_post")

## G ##

# Pre
Anova(aov(G_t1~group)) #F(2,24)=2.23, p=0.13        
shapiro.test(residuals(aov(G_t1~group))) #p=0.31

# Post
Anova(aov(G_t2~group)) #F(2,24)=0.75, p=0.48        
shapiro.test(residuals(aov(G_t2~group))) #p=0.76

# Difference post-pre
results$G_diff=G_t2-G_t1; detach(results); attach(results)
Anova(aov(G_diff~group)) #F(2,24)=3.41, p=0.0496
plot(G_diff~group)
shapiro.test(residuals(aov(G_diff~group))) # p=0.0178
kruskal.test(G_diff~group) #X²(2)=3.65, p=0.16
leveneTest(G_diff~group) #F(2,24)=1.07,p=0.36

# Difference score =/= 0?
t.test(G_diff) #t(26)=1.24, p=0.23


#--------------------------------------------------------------------------------------------------------------------#

## JiT (uncorrected!)

## Tests

# Pre
Anova(aov(JiT_t1~group)) #F(2,24)=5.50, p=0.0108
shapiro.test(residuals(aov(JiT_t1~group))) #p=0.03
kruskal.test(JiT_t1~group) #X²(2)=7.78, p=0.0204
leveneTest(JiT_t1~group) #F(2,24)=4.12, p=0.0291

# Post
Anova(aov(JiT_t2~group)) #F(2,24)=5.76, p=0.0091
shapiro.test(residuals(aov(JiT_t2~group))) #p=0.10
TukeyHSD(aov(JiT_t2~group))
#  Tukey multiple comparisons of means
#95% family-wise confidence level
#Fit: aov(formula = JiT_t2 ~ group)
#$group
#               diff         lwr        upr     p adj
#MEN-CON  0.38740612  0.08281592  0.691996309 0.0109336
#GLI-CON  0.05072155 -0.28492050  0.386363594 0.9247434
#GLI-MEN -0.33668457 -0.67232661 -0.001042521 0.0491887
leveneTest(JiT_t2~group) #F(2,24)=9.69, p=0.0008

# Difference post-pre
results$JiT_diff=JiT_t2-JiT_t1; detach(results); attach(results)
Anova(aov(JiT_diff~group)) #F(2,24)=0.79, p=0.47
shapiro.test(residuals(aov(JiT_diff~group))) #p=0.0140
kruskal.test(JiT_diff~group) #X²(2)=2.71, p=0.26
leveneTest(JiT_diff~group) #F(2,24)=2.76, p=0.08

# Difference score =/= 0?
t.test(JiT_diff) #t(26)=0.27, p=0.79


#--------------------------------------------------------------------------------------------------------------------#

## JiNT (uncorrected!)

## Tests

# Pre
Anova(aov(JiNT_t1~group)) #F(2,24)=1.9, p=0.17
shapiro.test(residuals(aov(JiNT_t1~group))) #p=0.72
leveneTest(aov(JiNT_t1~group)) #F(2,24)=0.55, p=0.58

# Post
Anova(aov(JiNT_t2~group)) #F(2,24)=1.15, p=0.33
shapiro.test(residuals(aov(JiNT_t2~group))) #p=0.89
leveneTest(aov(JiNT_t2~group)) #F(2,24)=1.10, p=0.35

# Difference score different by group
results$JiNT_diff=JiNT_t2-JiNT_t1; detach(results); attach(results)
Anova(aov(JiNT_diff~group)) #F(2,24)=0.18, p=0.84
shapiro.test(residuals(aov(JiNT_diff~group))) #p=0.61
leveneTest(aov(JiNT_diff~group)) #F(2,24)=0.22, p=0.81

# Difference =/= 0?
t.test(JiNT_diff) #t(26)=1.21, p=0.24


#--------------------------------------------------------------------------------------------------------------------#

## JiT (corrected!)

## Tests

# Pre
Anova(aov(JiT_res_t1~group)) #F(2,24)=12.16, p=0.0002
shapiro.test(residuals(aov(JiT_res_t1~group))) #p=0.02
kruskal.test(JiT_res_t1~group) #X²(2)=18.95, p<.0001
leveneTest(JiT_res_t1~group) #F(2,24)=153.3, p<.0001

# Post
Anova(aov(JiT_res_t2~group)) #F(2,24)=10.74, p=0.0005
shapiro.test(residuals(aov(JiT_res_t2~group))) #p=0.04
kruskal.test(JiT_res_t2~group)#X²(2)=15.28, p=0.0005
leveneTest(JiT_res_t2~group) #F(2,24)=62.07, p<.0001

# Difference post-pre
results$JiT_res_diff=JiT_res_t2-JiT_res_t1; detach(results); attach(results)
Anova(aov(JiT_res_diff~group)) #F(2,24)=1.67, p=0.21
shapiro.test(residuals(aov(JiT_res_diff~group))) #p=0.04
kruskal.test(JiT_res_diff~group) #X²(2)=2.89,p=0.24
leveneTest(JiT_res_diff~group) #F(2,24)=2.00, p=0.16

# Difference score =/= 0?
t.test(JiT_res_diff) #t(26)=1.58, p=0.13



#--------------------------------------------------------------------------------------------------------------------#

## JiNT (corrected!)

## Tests

# Pre
Anova(aov(JiNT_res_t1~group)) #F(2,24)=4.58, p=0.0207
shapiro.test(residuals(aov(JiNT_res_t1~group))) #p=0.14
leveneTest(aov(JiNT_res_t1~group)) #F(2,24)=1.53, p=0.24

# Post
Anova(aov(JiNT_res_t2~group)) #F(2,24)=7.82, p=0.0024
shapiro.test(residuals(aov(JiNT_res_t2~group))) #p=0.64
leveneTest(aov(JiNT_res_t2~group)) #F(2,24)=1.57, p=0.23

# Difference score different by group
results$JiNT_res_diff=JiNT_res_t2-JiNT_res_t1; detach(results); attach(results)
Anova(aov(JiNT_res_diff~group)) #F(2,24)=0.66, p=0.53
shapiro.test(residuals(aov(JiNT_res_diff~group))) #p=0.23
leveneTest(aov(JiNT_res_diff~group)) #F(2,24)=0.25, p=0.78

# Difference =/= 0?
t.test(JiNT_res_diff) #t(26)=2.84, p=0.0087

#-----------------------------------------------------------------------------

## Plot 

## G

# 1) Scatter plot pre vs. post with main diagonal
results$group=factor(results$group, ordered=TRUE, levels=c('MEN', 'GLI', 'CON'))
detach(results); attach(results)

range(G_t1)
range(G_t2)

p1a<-ggplot(results, aes(x=G_t1,y=G_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(1.25,2.5) + ylim(1.25,2.5) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = "G post", x = "G pre") + 
  scale_fill_brewer(palette="Paired") + #, breaks=c("CON", "MEN", "GLI")) 
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p1b<-ggplot(results, aes(x=group, y=G_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group)) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8) + 
  #scale_x_discrete(labels = "")  
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = "G difference") 


## JiT res

# 1) Scatter plot pre vs. post with main diagonal
range(JiT_res_t1)
range(JiT_res_t2)

p2a<-ggplot(results, aes(x=JiT_res_t1,y=JiT_res_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-1.5,0.8) + ylim(-1.5,0.8) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = expression('J'[tumor]*' post'), x = expression('J'[tumor]*' pre')) + 
  scale_fill_brewer(palette="Paired") +#, breaks=c("CON", "MEN", "GLI")) 
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p2b<-ggplot(results, aes(x=group, y=JiT_res_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = expression('J'[tumor]*' difference')) 


## JiNT res

# 1) Scatter plot pre vs. post with main diagonal
range(JiNT_res_t1)
range(JiNT_res_t2)

p3a<-ggplot(results, aes(x=JiNT_res_t1,y=JiNT_res_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(0.48,0.72) + ylim(0.48,0.72) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = expression('J'[non-tumor]*' post'), x = expression('J'[non-tumor]*' pre')) + 
  scale_fill_brewer(palette="Paired") +#, breaks=c("CON", "MEN", "GLI")) 
  scale_shape_manual(values=c(21,24,22))


# 2) Box plot post-pre differences per group
p3b<-ggplot(results, aes(x=group, y=JiNT_res_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = expression('J'[non-tumor]*' difference')) 


## JiT (uncorrected)

# 1) Scatter plot pre vs. post with main diagonal
range(JiT_t1)
range(JiT_t2)

p4a<-ggplot(results, aes(x=JiT_t1,y=JiT_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(1.2,2.5) + ylim(1.2,2.5) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = expression('J'[tumor]*' post'), x = expression('J'[tumor]*' pre')) + 
  scale_fill_brewer(palette="Paired") + #, breaks=c("CON", "MEN", "GLI")) 
  scale_shape_manual(values=c(21,24,22))


# 2) Box plot post-pre differences per group
p4b<-ggplot(results, aes(x=group, y=JiT_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = expression('J'[tumor]*' difference')) 


## JiNT (uncorrected)

# 1) Scatter plot pre vs. post with main diagonal
range(JiNT_t1)
range(JiNT_t2)

p5a<-ggplot(results, aes(x=JiNT_t1,y=JiNT_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(1.2,1.55) + ylim(1.2,1.55) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = expression('J'[non-tumor]*' post'), x = expression('J'[non-tumor]*' pre')) + 
  scale_fill_brewer(palette="Paired") + #, breaks=c("CON", "MEN", "GLI")) 
  scale_shape_manual(values=c(21,24,22))


# 2) Box plot post-pre differences per group
p5b<-ggplot(results, aes(x=group, y=JiNT_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = expression('J'[non-tumor]*' difference')) 

png('results_modelparams_JRes_all2.png', width=1000, height=1200)
grid.arrange(p2b,p2a,p3b,p3a,p1b,p1a,ncol=2)
dev.off()

png('results_modelparams_JRaw_all2.png', width=1000, height=800)
grid.arrange(p4b,p4a,p5b,p5a,ncol=2)
dev.off()

rm(list=c('p1a','p1b','p2a','p2b','p3a','p3b','p4a','p4b','p5a','p5b'))


##################################################################################################################
###                                         PART 3: COGNITION DESCRIPTIVES                                     ###
##################################################################################################################
  
setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_cantab")

## Correlations with other covariates

# Continuous covariates
cbind(colnames(results))
#png('Cognition_covariates_t2.png')
corrplot(cor(results[,c(43:46,4,11,70)], use="pairwise"), type="lower", method="color", tl.col="black",
         diag=F, addCoef.col="black")
#dev.off()
#--> not so important: lesion volume

# Comparison with t1
par(mfrow=c(1,1))
#png('Cognition_covariates_t1.png')
corrplot(cor(results[,c(38:41,4,10,69)], use="pairwise"), type="lower", method="color", tl.col="black",
         diag=F, addCoef.col="black")
#dev.off()
#--> also at t1 lesion volume not so important

#png('Cognition_by_sex_t2.png')
par(mfrow=c(2,2))
boxplot(RTI_fiveRT_t2~sex); title('Reaction time')
boxplot(RVP_A_t2~sex); title('Sustained attention')
boxplot(SOC_prob_minmoves_t2~sex); title('Planning')
boxplot(SSP_spanlength_t2~sex); title('Working memory')
#dev.off()
#--> sex doesn't seem so important

# Okay so did it make a difference at t1? 
#png('Cognition_by_sex_t1.png')
par(mfrow=c(2,2))
boxplot(RTI_fiveRT_t1~sex); title('Reaction time')
boxplot(RVP_A_t1~sex); title('Sustained attention')
boxplot(SOC_prob_minmoves_t1~sex); title('Planning')
boxplot(SSP_spanlength_t1~sex); title('Working memory')
#dev.off()
#--> maybe a little...


##################################################################################################################
###                                        PART 4: SC TOPOLOGY DESCRIPTIVES                                    ###
##################################################################################################################

setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_GTA")

# Pairwise correlations among graph theory measures
results_SC = results[,c(58:68)]
colnames(results_SC)=
  c("density", "clustering", "Eloc", "Q", "Eglob", "comm", "degree", "strength", "EBC", "BC", "PC")
par(mfrow=c(1,1))
#png('GTAmetrics_corplot_t2.png')
corrplot(cor(results_SC), type='lower', diag=FALSE, method='color', tl.col='black', addCoef.col = 'black')
#dev.off()

# Remove variables with correlations > 0.80:
#Density & degree highly related --> remove degree
#Clust & Eloc highly related --> remove clust
#EBC & BC highly related --> remove EBC
#Eglob & comm highly related --> remove comm
#Eglob & strength highly related --> remove strength
#Eloc & Eglob highly related --> remove Eloc

# PCA
palette(c("blue4", "firebrick"))
SC_pca = prcomp(results_SC, center=TRUE, scale=TRUE)
#png('SC_pca1.png')
biplot(SC_pca, xlabs=rep("*", nrow(results)), cex=1.3)
#dev.off()

results_SC = results[,c(61,62,68)]
colnames(results_SC)=c("Q", "Eglob", "PC")
SC_pca = prcomp(results_SC, center=TRUE, scale=TRUE)
#png('SC_pca2.png')
biplot(SC_pca, xlabs=rep("*", nrow(results)), cex=1.3)
#dev.off()

#--> FINAL: Q, Eglob, PC (most segregated in PCA space)
rm(results_SC)
rm(SC_pca)


## Correlations with other covariates

cbind(colnames(results))
#png('GTAmetrics_covariates_t2.png')
corrplot(cor(results[,c(61,62,68,4,7,11,13,15,70)], use="pairwise"), type="lower", method="color", tl.col="black",
         diag=F, addCoef.col="black")
#dev.off()
#--> not so important: intensity normalization factor, motion
#png('GTAmetrics_t2_BySex.png')
par(mfrow=c(1,3))
boxplot(SC_Q_t2~sex); title('Q')
boxplot(SC_Eglob_t2~sex); title('Eglob')
boxplot(SC_PC_t2~sex); title('PC')
#dev.off()
#--> include sex


###################################################################################################################
###                                             PART 5: CLEAN DATASET                                           ###
###################################################################################################################

# RTI

lm_RTI_t2 = lm(RTI_fiveRT_t2 ~ STAI_t2 + age + sex + lesion_vol_cm3_t2)
summary(lm_RTI_t2)
shapiro.test(residuals(lm_RTI_t2)) 
par(mfrow=c(1,1))
hist(residuals(lm_RTI_t2))
results$RTI_res_t2 = residuals(lm_RTI_t2)

lm_RTI_t1 = lm(RTI_fiveRT_t1 ~ STAI_t1 + age + sex + lesion_vol_cm3)
summary(lm_RTI_t1)
shapiro.test(residuals(lm_RTI_t1)) 
results$RTI_res_t1 = residuals(lm_RTI_t1)

results$RTI_res_diff = results$RTI_res_t2 - results$RTI_res_t1

results$RTI_resnorm_t1 = (results$RTI_res_t1 - mean(results$RTI_res_t1[c(1:11)])) / 
  sqrt(var(results$RTI_res_t1[c(1:11)]))

results$RTI_resnorm_t2 = (results$RTI_res_t2 - mean(results$RTI_res_t1[c(1:11)])) / 
  sqrt(var(results$RTI_res_t1[c(1:11)]))

results$RTI_resnorm_diff = results$RTI_resnorm_t2-results$RTI_resnorm_t1

rm(list=c('lm_RTI_t2', 'lm_RTI_t1'))


# RVP

lm_RVP_t2 = lm(RVP_A_t2 ~ STAI_t2 + age + sex + lesion_vol_cm3_t2)
summary(lm_RVP_t2)
shapiro.test(residuals(lm_RVP_t2)) 
results$RVP_res_t2 = residuals(lm_RVP_t2)

lm_RVP_t1 = lm(RVP_A_t1 ~ STAI_t1 + age + sex + lesion_vol_cm3)
summary(lm_RVP_t1)
shapiro.test(residuals(lm_RVP_t1)) #p.03...
hist(residuals(lm_RVP_t1))
results$RVP_res_t1 = c(residuals(lm_RVP_t1)[1:11], NA, residuals(lm_RVP_t1)[12:26])

results$RVP_res_diff = results$RVP_res_t2 - results$RVP_res_t1

results$RVP_resnorm_t1 = (results$RVP_res_t1 - mean(results$RVP_res_t1[c(1:11)])) / 
  sqrt(var(results$RVP_res_t1[c(1:11)]))

results$RVP_resnorm_t2 = (results$RVP_res_t2 - mean(results$RVP_res_t1[c(1:11)])) / 
  sqrt(var(results$RVP_res_t1[c(1:11)]))

results$RVP_resnorm_diff = results$RVP_resnorm_t2-results$RVP_resnorm_t1

rm(list=c('lm_RVP_t2', 'lm_RVP_t1'))


# SOC

lm_SOC_t2 = lm(SOC_prob_minmoves_t2 ~ STAI_t2 + age + sex + lesion_vol_cm3_t2)
summary(lm_SOC_t2)
shapiro.test(residuals(lm_SOC_t2)) 
hist(residuals(lm_SOC_t2))
results$SOC_res_t2 = residuals(lm_SOC_t2)

lm_SOC_t1 = lm(SOC_prob_minmoves_t1 ~ STAI_t1 + age + sex + lesion_vol_cm3)
summary(lm_SOC_t1)
shapiro.test(residuals(lm_SOC_t1)) 
results$SOC_res_t1 = residuals(lm_SOC_t1)

results$SOC_res_diff = results$SOC_res_t2 - results$SOC_res_t1

results$SOC_resnorm_t1 = (results$SOC_res_t1 - mean(results$SOC_res_t1[c(1:11)])) / 
  sqrt(var(results$SOC_res_t1[c(1:11)]))

results$SOC_resnorm_t2 = (results$SOC_res_t2 - mean(results$SOC_res_t1[c(1:11)])) / 
  sqrt(var(results$SOC_res_t1[c(1:11)]))

results$SOC_resnorm_diff = results$SOC_resnorm_t2 - results$SOC_resnorm_t1

rm(list=c('lm_SOC_t2', 'lm_SOC_t1'))


# SSP

lm_SSP_t2 = lm(SSP_spanlength_t2 ~ STAI_t2 + age + sex + lesion_vol_cm3_t2)
summary(lm_SSP_t2)
shapiro.test(residuals(lm_SSP_t2)) 
results$SSP_res_t2 = residuals(lm_SSP_t2)

lm_SSP_t1 = lm(SSP_spanlength_t1 ~ STAI_t1 + age + sex + lesion_vol_cm3)
summary(lm_SSP_t1)
shapiro.test(residuals(lm_SSP_t1)) 
results$SSP_res_t1 = residuals(lm_SSP_t1)

results$SSP_res_diff = results$SSP_res_t2 - results$SSP_res_t1

results$SSP_resnorm_t1 = (results$SSP_res_t1 - mean(results$SSP_res_t1[c(1:11)])) / 
  sqrt(var(results$SSP_res_t1[c(1:11)]))

results$SSP_resnorm_t2 = (results$SSP_res_t2 - mean(results$SSP_res_t1[c(1:11)])) / 
  sqrt(var(results$SSP_res_t1[c(1:11)]))

results$SSP_resnorm_diff = results$SSP_resnorm_t2 - results$SSP_resnorm_t1

rm(list=c('lm_SSP_t2', 'lm_SSP_t1'))


# Q

lm_Q_t2 = lm(SC_Q_t2 ~ STAI_t2 + age + sex + hand_t2 + lesion_vol_cm3_t2 + intnorm_scaling_t2 + 
               fmri_MD_t2)
summary(lm_Q_t2)
shapiro.test(residuals(lm_Q_t2)) 
results$SC_Q_res_t2 = residuals(lm_Q_t2)

lm_Q_t1 = lm(SC_Q_t1 ~ STAI_t1 + age + sex + hand_t1 + lesion_vol_cm3 + intnorm_scaling_t1 + 
               fmri_MD_t1)
summary(lm_Q_t1)
shapiro.test(residuals(lm_Q_t1)) 
results$SC_Q_res_t1 = residuals(lm_Q_t1)

results$SC_Q_res_diff = results$SC_Q_res_t2 - results$SC_Q_res_t1

results$SC_Q_resnorm_t1 = (results$SC_Q_res_t1 - mean(results$SC_Q_res_t1[c(1:11)])) / 
  sqrt(var(results$SC_Q_res_t1[c(1:11)]))

results$SC_Q_resnorm_t2 = (results$SC_Q_res_t2 - mean(results$SC_Q_res_t1[c(1:11)])) / 
  sqrt(var(results$SC_Q_res_t1[c(1:11)]))

results$SC_Q_resnorm_diff = results$SC_Q_resnorm_t2 - results$SC_Q_resnorm_t1

rm(list=c('lm_Q_t2', 'lm_Q_t1'))


# Eglob

lm_Eglob_t2 = lm(SC_Eglob_t2 ~ STAI_t2 + age + sex + hand_t2 + lesion_vol_cm3_t2 + 
                   intnorm_scaling_t2 + fmri_MD_t2)
summary(lm_Eglob_t2)
shapiro.test(residuals(lm_Eglob_t2)) 
results$SC_Eglob_res_t2 = residuals(lm_Eglob_t2)

lm_Eglob_t1 = lm(SC_Eglob_t1 ~ STAI_t1 + age + sex + hand_t1 + lesion_vol_cm3 + 
                   intnorm_scaling_t1 + fmri_MD_t1)
summary(lm_Eglob_t1)
shapiro.test(residuals(lm_Eglob_t1)) #not really okay p=0.003
hist(residuals(lm_Eglob_t1)) #that's fine
results$SC_Eglob_res_t1 = residuals(lm_Eglob_t1)

results$SC_Eglob_res_diff = results$SC_Eglob_res_t2 - results$SC_Eglob_res_t1

results$SC_Eglob_resnorm_t1 = (results$SC_Eglob_res_t1 - mean(results$SC_Eglob_res_t1[c(1:11)])) / 
  sqrt(var(results$SC_Eglob_res_t1[c(1:11)]))

results$SC_Eglob_resnorm_t2 = (results$SC_Eglob_res_t2 - mean(results$SC_Eglob_res_t1[c(1:11)])) / 
  sqrt(var(results$SC_Eglob_res_t1[c(1:11)]))

results$SC_Eglob_resnorm_diff = results$SC_Eglob_resnorm_t2 - results$SC_Eglob_resnorm_t1

rm(list=c('lm_Eglob_t2', 'lm_Eglob_t1'))


# PC

lm_PC_t2 = lm(SC_PC_t2 ~ STAI_t2 + age + sex + hand_t2 + lesion_vol_cm3_t2 + 
                   intnorm_scaling_t2 + fmri_MD_t2)
summary(lm_PC_t2)
shapiro.test(residuals(lm_PC_t2)) 
results$SC_PC_res_t2 = residuals(lm_PC_t2)

lm_PC_t1 = lm(SC_PC_t1 ~ STAI_t1 + age + sex + hand_t1 + lesion_vol_cm3 + 
                   intnorm_scaling_t1 + fmri_MD_t1)
summary(lm_PC_t1)
shapiro.test(residuals(lm_PC_t1)) 
results$SC_PC_res_t1 = residuals(lm_PC_t1)

results$SC_PC_res_diff = results$SC_PC_res_t2 - results$SC_PC_res_t1

results$SC_PC_resnorm_t1 = (results$SC_PC_res_t1 - mean(results$SC_PC_res_t1[c(1:11)])) / 
  sqrt(var(results$SC_PC_res_t1[c(1:11)]))

results$SC_PC_resnorm_t2 = (results$SC_PC_res_t2 - mean(results$SC_PC_res_t1[c(1:11)])) / 
  sqrt(var(results$SC_PC_res_t1[c(1:11)]))

results$SC_PC_resnorm_diff = results$SC_PC_resnorm_t2 - results$SC_PC_resnorm_t1

rm(list=c('lm_PC_t2', 'lm_PC_t1'))

detach(results); attach(results)

# Save dataset
setwd('/home/hannelore/Documents/ANALYSES/TVB_post')
write.table(x=results, file='RESULTS_ALL_afterprep_20190411.csv',
            quote=TRUE, sep=';', dec='.', row.names=FALSE)



###################################################################################################################
###                               PART 6: GROUP DIFFERENCES IN COGNITION & SC TOPOLOGY                          ###
###################################################################################################################

results=read.table(file='RESULTS_ALL_afterprep_20190411.csv', header=TRUE, 
                   sep=";")
attach(results)

setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_cantab")


## --- RTI

# Pre
Anova(aov(RTI_fiveRT_t1~group)) #F(2,24)=0.72, p=0.50
shapiro.test(residuals(aov(RTI_fiveRT_t1~group))) #p=0.25

# Post
Anova(aov(RTI_fiveRT_t2~group)) #F(2,24)=0.28, p=0.76
shapiro.test(residuals(aov(RTI_fiveRT_t2~group))) #p<0.0001
kruskal.test(RTI_fiveRT_t2~group) #X²(2)=0.38, p=0.83

# Difference
RTI_diff=RTI_fiveRT_t2-RTI_fiveRT_t1
Anova(aov(RTI_diff~group)) #F(2,24)=0.06, p=0.94
shapiro.test(residuals(aov(RTI_diff~group))) #p<0.0001
kruskal.test(RTI_diff~group) #X²(2)=0.28, p=0.87

t.test(RTI_diff) #t(26)=-2.47, p=0.0206


## --- RTI resnorm

# Pre
Anova(aov(RTI_resnorm_t1~group)) #F(2,24)=0.59, p=0.56
shapiro.test(residuals(aov(RTI_resnorm_t1~group))) #p=0.44

# Post
Anova(aov(RTI_resnorm_t2~group)) #F(2,24)=0.13, p=0.88
shapiro.test(residuals(aov(RTI_resnorm_t2~group))) #p=0.0015
kruskal.test(RTI_resnorm_t2~group) #X²(2)=0.23, p=0.89

# Difference
Anova(aov(RTI_resnorm_diff~group)) #F(2,24)=0.17, p=0.85
shapiro.test(residuals(aov(RTI_resnorm_diff~group))) #p=0.0047
kruskal.test(RTI_resnorm_diff~group) #X²(2)=0.49, p=0.78

t.test(RTI_resnorm_diff) #t(26)=0, p=1


## --- RVP

# Pre
Anova(aov(RVP_A_t1~group)) #F(2,23)=0.21, p=0.82
shapiro.test(residuals(aov(RVP_A_t1~group))) #p=0.51

# Post
Anova(aov(RVP_A_t2~group)) #F(2,24)=0.40, p=0.68
shapiro.test(residuals(aov(RVP_A_t2~group))) #p=0.0016
kruskal.test(RVP_A_t2~group) #KW(2)=0.10, p=0.95

# Difference
RVP_diff=RVP_A_t2-RVP_A_t1
Anova(aov(RVP_diff~group)) #F(2,23)=0.50, p=0.61
shapiro.test(residuals(aov(RVP_diff~group))) #p=0.83

t.test(RVP_diff) #t(25)=0, p=1


## --- RVP resnorm

# Pre
Anova(aov(RVP_resnorm_t1~group)) #F(2,23)=0.74, p=0.49
shapiro.test(residuals(aov(RVP_resnorm_t1~group))) #p=0.28

# Post
Anova(aov(RVP_resnorm_t2~group)) #F(2,24)=0.03, p=0.97
shapiro.test(residuals(aov(RVP_resnorm_t2~group))) #p=0.34

# Difference
Anova(aov(RVP_resnorm_diff~group)) #F(2,23)=1.74, p=0.20
shapiro.test(residuals(aov(RVP_resnorm_diff~group))) #p=0.25

t.test(RVP_resnorm_diff) #t(25)=0.72, p=0.48


#--- SOC

# Pre
Anova(aov(SOC_prob_minmoves_t1~group)) #F(2,24)=0.08, p=0.92
shapiro.test(residuals(aov(SOC_prob_minmoves_t1~group))) #p=0.0210
kruskal.test(SOC_prob_minmoves_t1~group) #KW(2)=0.09, p=0.95

# Post
Anova(aov(SOC_prob_minmoves_t2~group)) #F(2,24)=1.47, p=0.25
shapiro.test(residuals(aov(SOC_prob_minmoves_t2~group))) #p=0.08

# Difference
SOC_diff=SOC_prob_minmoves_t2-SOC_prob_minmoves_t1
Anova(aov(SOC_diff~group)) #F(2,24)=0.73, p=0.49
shapiro.test(residuals(aov(SOC_diff~group))) #p=0.45

t.test(SOC_diff) #t(26)=1.32, p=0.20


## --- SOC resnorm

# Pre
Anova(aov(SOC_resnorm_t1~group)) #F(2,24)=0.45, p=0.64
shapiro.test(residuals(aov(SOC_resnorm_t1~group))) #p=0.91

# Post
Anova(aov(SOC_resnorm_t2~group)) #F(2,24)=0.64, p=0.53
shapiro.test(residuals(aov(SOC_resnorm_t2~group))) #p=0.07

# Difference
Anova(aov(SOC_resnorm_diff~group)) #F(2,24)=0.03, p=0.97
shapiro.test(residuals(aov(SOC_resnorm_diff~group))) #p=0.11

t.test(SOC_resnorm_diff) #t(26)=0, p=1


## --- SSP

# Pre
Anova(aov(SSP_spanlength_t1~group)) #F(2,24)=0.94, p=0.41
shapiro.test(residuals(aov(SSP_spanlength_t1~group))) #p=0.58

# Post
Anova(aov(SSP_spanlength_t2~group)) #F(2,24)=0.35, p=0.71
shapiro.test(residuals(aov(SSP_spanlength_t2~group))) #p=0.65

# Difference
SSP_diff=SSP_spanlength_t2-SSP_spanlength_t1
Anova(aov(SSP_diff~group)) #F(2,24)=1.49, p=0.24
shapiro.test(residuals(aov(SSP_diff~group))) #p=0.23

t.test(SSP_diff) #t(26)=1.14, p=0.26


## --- SSP resnorm

# Pre
Anova(aov(SSP_resnorm_t1~group)) #F(2,24)=1.23, p=0.31
shapiro.test(residuals(aov(SSP_resnorm_t1~group))) #p=0.68

# Post
Anova(aov(SSP_resnorm_t2~group)) #F(2,24)=0.27, p=0.76
shapiro.test(residuals(aov(SSP_resnorm_t2~group))) #p=0.24

# Difference
Anova(aov(SSP_resnorm_diff~group)) #F(2,24)=2.89, p=0.0751
shapiro.test(residuals(aov(SSP_resnorm_diff~group))) #p=0.28
TukeyHSD(aov(SSP_resnorm_diff~group))

t.test(SSP_resnorm_diff) #t(26)=0, p=1


## Plots

# RTI resnorm

# 1) Scatter plot pre vs. post with main diagonal
results$group=factor(results$group, ordered=TRUE, levels=c('MEN', 'GLI', 'CON'))
detach(results); attach(results)

range(RTI_resnorm_t1)
range(RTI_resnorm_t2)

p1a<-ggplot(results, aes(x=RTI_resnorm_t1,y=RTI_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-2.7,6.9) + ylim(-2.7,6.9) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Reaction time post', x = 'Reaction time pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))


# 2) Box plot post-pre differences per group
p1b<-ggplot(results, aes(x=group, y=RTI_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Reaction time difference') 


# RVP resnorm

# 1) Scatter plot pre vs. post with main diagonal
range(RVP_resnorm_t1,na.rm=T)
range(RVP_resnorm_t2,na.rm=T)

p2a<-ggplot(results, aes(x=RVP_resnorm_t1,y=RVP_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-2.7,2) + ylim(-2.7,2) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Sustained attention post', x = 'Sustained attention pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p2b<-ggplot(results, aes(x=group, y=RVP_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Sustained attention difference') 


# SOC resnorm

# 1) Scatter plot pre vs. post with main diagonal
range(SOC_resnorm_t1,na.rm=T)
range(SOC_resnorm_t2,na.rm=T)

p3a<-ggplot(results, aes(x=SOC_resnorm_t1,y=SOC_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-3.9,2.1) + ylim(-3.9,2.1) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Planning accuracy post', x = 'Planning accuracy pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p3b<-ggplot(results, aes(x=group, y=SOC_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Planning accuracy difference') 


# SSP resnorm

# 1) Scatter plot pre vs. post with main diagonal
range(SSP_resnorm_t1,na.rm=T)
range(SSP_resnorm_t2,na.rm=T)

p4a<-ggplot(results, aes(x=SSP_resnorm_t1,y=SSP_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-2.3,1.6) + ylim(-2.3,1.6) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Working memory post', x = 'Working memory pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p4b<-ggplot(results, aes(x=group, y=SSP_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Working memory difference') 

png('results_cognition_all2.png', width=1000, height=1500)
grid.arrange(p1b,p1a,p2b,p2a,p3b,p3a,p4b,p4a,ncol=2)
dev.off()

rm(list=c('p1a','p1b','p2a','p2b','p3a','p3b','p4a','p4b'))


#-------------------------------------------------------------------------------------------------------------------

setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_GTA")

## --- Q

# Pre
Anova(aov(SC_Q_t1~group)) #F(2,24)=0.61, p=0.55
shapiro.test(residuals(aov(SC_Q_t1~group))) #p=0.67

# Post
Anova(aov(SC_Q_t2~group)) #F(2,24)=0.58, p=0.57
shapiro.test(residuals(aov(SC_Q_t2~group))) #p=0.37

# Difference
Q_diff=SC_Q_t2-SC_Q_t1
Anova(aov(Q_diff~group)) #F(2,24)=0.09, p=0.91
shapiro.test(residuals(aov(Q_diff~group))) #p=0.23

t.test(Q_diff) #t(26)=-21.80, p<.0001


## --- Q resnorm

# Pre
Anova(aov(SC_Q_resnorm_t1~group)) #F(2,24)=0.75, p=0.48
shapiro.test(residuals(aov(SC_Q_resnorm_t1~group))) #p=0.91

# Post
Anova(aov(SC_Q_resnorm_t2~group)) #F(2,24)=1.07, p=0.36
shapiro.test(residuals(aov(SC_Q_resnorm_t2~group))) #p=0.69

# Difference
Anova(aov(SC_Q_resnorm_diff~group)) #F(2,24)=1.26, p=0.30
shapiro.test(residuals(aov(SC_Q_resnorm_diff~group))) #p=0.56

t.test(SC_Q_resnorm_diff) #t(26)=0, p=1


## --- Eglob

# Pre
Anova(aov(SC_Eglob_t1~group)) #F(2,24)=0.61, p=0.55
shapiro.test(residuals(aov(SC_Eglob_t1~group))) #p=0.24

# Post
Anova(aov(SC_Eglob_t2~group)) #F(2,24)=0.33, p=0.72
shapiro.test(residuals(aov(SC_Eglob_t2~group))) #p=0.0403
kruskal.test(SC_Eglob_t2~group) #X²(2)=1.34, p=0.51

# Difference
Eglob_diff=SC_Eglob_t2-SC_Eglob_t1
Anova(aov(Eglob_diff~group)) #F(2,24)=0.96, p=0.40
shapiro.test(residuals(aov(Eglob_diff~group))) #p=0.0017
kruskal.test(Eglob_diff~group) #X²(2)=2.59, p=0.27

t.test(Eglob_diff) #t(26)=15.62, p<.0001


##--- Eglob resnorm

# Pre
Anova(aov(SC_Eglob_resnorm_t1~group)) #F(2,24)=0.65, p=0.53
shapiro.test(residuals(aov(SC_Eglob_resnorm_t1~group))) #p=0.0003
kruskal.test(SC_Eglob_resnorm_t1~group) #X²(2)=4.92, p=0.0855

# Post
Anova(aov(SC_Eglob_resnorm_t2~group)) #F(2,24)=0.46, p=0.64
shapiro.test(residuals(aov(SC_Eglob_resnorm_t2~group))) #p=0.07
kruskal.test(SC_Eglob_resnorm_t2~group) #X²(2)=2.41, p=0.30

# Difference
Anova(aov(SC_Eglob_resnorm_diff~group)) #F(2,24)=0.95, p=0.40
shapiro.test(residuals(aov(SC_Eglob_resnorm_diff~group))) #p=0.0186
kruskal.test(SC_Eglob_resnorm_diff~group) #X²(2)=2.62, p=0.27

t.test(SC_Eglob_resnorm_diff) #t(26)=0, p=1



##--- PC

# Pre
Anova(aov(SC_PC_t1~group)) #F(2,24)=1.22, p=0.31
shapiro.test(residuals(aov(SC_PC_t1~group))) #p=0.23

# Post
Anova(aov(SC_PC_t2~group)) #F(2,24)=1.17, p=0.33
shapiro.test(residuals(aov(SC_PC_t2~group))) #p=0.71

# Difference
PC_diff=SC_PC_t2-SC_PC_t1
Anova(aov(PC_diff~group)) #F(2,24)=3.88, p=0.0347
shapiro.test(residuals(aov(PC_diff~group))) #p=0.77
TukeyHSD(aov(PC_diff~group)) 
#               diff         lwr          upr     p adj
#MEN-CON  0.00445000 -0.02685685  0.035756850 0.9330947
#GLI-CON -0.03178286 -0.06628133  0.002715611 0.0749673
#GLI-MEN -0.03623286 -0.07073133 -0.001734389 0.0382436
t.test(PC_diff) #t(26)=13.61, p<.0001


##--- PC resnorm

## Tests

# Pre
Anova(aov(SC_PC_resnorm_t1~group)) #F(2,24)=2.60, p=0.0947
shapiro.test(residuals(aov(SC_PC_resnorm_t1~group))) #p=0.90

# Post
Anova(aov(SC_PC_resnorm_t2~group)) #F(2,24)=0.48, p=0.62
shapiro.test(residuals(aov(SC_PC_resnorm_t2~group))) #p=0.76

# Difference
Anova(aov(SC_PC_resnorm_diff~group)) #F(2,24)=1.26, p=0.30
shapiro.test(residuals(aov(SC_PC_resnorm_diff~group))) #p=0.81

t.test(SC_PC_resnorm_diff) #t(26)=0, p=1


## Plots

# Q resnorm

# 1) Scatter plot pre vs. post with main diagonal
results$group=factor(results$group, ordered=TRUE, levels=c('MEN', 'GLI', 'CON'))
detach(results); attach(results)

range(SC_Q_resnorm_t1,na.rm=T)
range(SC_Q_resnorm_t2,na.rm=T)

p1a<-ggplot(results, aes(x=SC_Q_resnorm_t1,y=SC_Q_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-2.7,3.5) + ylim(-2.7,3.5) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Modularity post', x = 'Modularity pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p1b<-ggplot(results, aes(x=group, y=SC_Q_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Modularity difference') 


# Eglob resnorm

# 1) Scatter plot pre vs. post with main diagonal
range(SC_Eglob_resnorm_t1,na.rm=T)
range(SC_Eglob_resnorm_t2,na.rm=T)

p2a<-ggplot(results, aes(x=SC_Eglob_resnorm_t1,y=SC_Eglob_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-5.7,5.6) + ylim(-5.7,5.6) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Global efficiency post', x = 'Global efficiency pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p2b<-ggplot(results, aes(x=group, y=SC_Eglob_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Global efficiency difference') 


# PC resnorm

# 1) Scatter plot pre vs. post with main diagonal
range(SC_PC_resnorm_t1,na.rm=T)
range(SC_PC_resnorm_t2,na.rm=T)

p3a<-ggplot(results, aes(x=SC_PC_resnorm_t1,y=SC_PC_resnorm_t2, fill=group, shape=group)) +
  geom_abline(slope = 1, linetype='longdash', color='darkgray', size=1) +
  geom_point(size=6) +
  xlim(-2.9,4.3) + ylim(-2.9,4.3) +
  theme(axis.ticks = element_blank(), text = element_text(size = 25), legend.position="none") +
  labs(y = 'Participation coefficient post', x = 'Participation coefficient pre') + 
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

# 2) Box plot post-pre differences per group
p3b<-ggplot(results, aes(x=group, y=SC_PC_resnorm_diff, fill=group)) + 
  geom_hline(yintercept = 0, linetype = 'longdash', color='darkgray', size=0.8) +
  geom_boxplot(aes(fill=group), outlier.shape=NA) +
  geom_jitter(shape=16, position=position_jitter(0.2), size=2) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25),
        legend.position="none",
        legend.key.size=unit(3,"line")) +
  labs(x = "", y = 'Participation coefficient difference') 

png('results_GTA_all2.png', width=1000, height=1200)
grid.arrange(p1b,p1a,p2b,p2a,p3b,p3a,ncol=2)
dev.off()
rm(list=c('p1a','p1b','p2a','p2b','p3a','p3b'))

