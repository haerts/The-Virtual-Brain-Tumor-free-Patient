###################################################################################################################
###                                                                                                             ###
###                   MODELING BRAIN DYNAMICS AFTER TUMOR RESECTION USING THE VIRTUAL BRAIN                     ###
###               ==============================================================================                ###
###             PART 2: Regression models to check effect of GTA/cognition on modeling parameters               ###
###                                                                                                             ###
### Created by Hannelore Aerts                                                                                  ###
### Date last update: 21/08/2019                                                                                ###
###################################################################################################################

# Analysis remarks:
# --> use Ji corrected for ROI size only, not ROI size + in-strength      
#     UPDATE: use uncorrected Ji!
#     UPDATE 2: repeat analyses while correcting for ROI size, just to check robustness of results 
# --> use J=JiT=JiNT for controls
# --> use second batch of TVBii analyses with normalization factor=pre-op normalization factor
#     = 75K (instead of 45K in first batch)
# --> don't correct for motivation in cognitive performance scores
# --> 11/04/2019: resnorm use !!pre-operative!! mean and sd 

library('car')
library('lsr')
library('corrplot')
library(Hmisc)
library(ggplot2)
library('gridExtra')


### Read in data and prepare for analyses ------------------------------------------------------------------------#

# Wide data format

setwd("/home/hannelore/Documents/ANALYSES/TVB_post")
results=read.table(file="RESULTS_ALL_afterprep_20190411.csv", header=TRUE, sep=";")
str(results, list.len=ncol(results))

results=within(results, {
  subID=as.character(subID)
  date_t1 = as.Date(date_t1, "%d/%m/%y")
  date_t2 = as.Date(date_t2, "%d/%m/%y")
  group = factor(group, ordered=TRUE, levels=c('CON', 'MEN', 'GLI'))
  fmri_TR_t1 = as.factor(fmri_TR_t1)
  fmri_TR_t2 = as.factor(fmri_TR_t2)
})

attach(results)



### Regression models for G --------------------------------------------------------------------------------------#

# Pre CANTAB
lm_G_cantab_t1 = lm(G_t1 ~ RTI_resnorm_t1 + RVP_resnorm_t1 + SOC_resnorm_t1 + SSP_resnorm_t1)
shapiro.test(residuals(lm_G_cantab_t1)) #p=0.47
summary(lm_G_cantab_t1)
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     1.660329   0.053013  31.320   <2e-16 ***
#RTI_resnorm_t1 -0.003114   0.068485  -0.045    0.964    
#RVP_resnorm_t1  0.023881   0.078728   0.303    0.765    
#SOC_resnorm_t1  0.025245   0.058615   0.431    0.671    
#SSP_resnorm_t1 -0.035520   0.074647  -0.476    0.639
rm(lm_G_cantab_t1)

# Post CANTAB
lm_G_cantab_t2 = lm(G_t2 ~ RTI_resnorm_t2 + RVP_resnorm_t2 + SOC_resnorm_t2 + SSP_resnorm_t2)
shapiro.test(residuals(lm_G_cantab_t2)) #p=0.09
summary(lm_G_cantab_t2)
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     1.71252    0.06390  26.799   <2e-16 ***
#RTI_resnorm_t2  0.04932    0.04090   1.206    0.241    
#RVP_resnorm_t2  0.09485    0.09850   0.963    0.346    
#SOC_resnorm_t2 -0.04361    0.06635  -0.657    0.518    
#SSP_resnorm_t2 -0.03203    0.09768  -0.328    0.746   
rm(lm_G_cantab_t2)

# Difference CANTAB
lm_G_cantab_diff = lm(G_diff ~ RTI_resnorm_diff + RVP_resnorm_diff + SOC_resnorm_diff + SSP_resnorm_diff)
shapiro.test(residuals(lm_G_cantab_diff)) #p=0.78
summary(lm_G_cantab_diff)
#C                  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       0.047953   0.045606   1.051   0.3050  
#RTI_resnorm_diff  0.008101   0.030379   0.267   0.7923  
#RVP_resnorm_diff  0.090689   0.061673   1.470   0.1563  
#SOC_resnorm_diff -0.118760   0.042776  -2.776   0.0113 *
#SSP_resnorm_diff  0.096766   0.054549   1.774   0.0906 .
etaSquared(lm_G_cantab_diff, anova=F, type=3)
#                         eta.sq eta.sq.part
#RTI_resnorm_diff 0.002138712 0.003374801
#RVP_resnorm_diff 0.065032995 0.093354499
#SOC_resnorm_diff 0.231827588 0.268499614
#SSP_resnorm_diff 0.094645158 0.130322819
rm(lm_G_cantab_diff)

# Repeat without outlier
tmp=results
tmp$SOC_resnorm_diff[11]<-NA
lm_G_cantab_diff = lm(G_diff ~ RTI_resnorm_diff + RVP_resnorm_diff + SOC_resnorm_diff + SSP_resnorm_diff,
                      data=tmp)
summary(lm_G_cantab_diff)
#C                  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       0.045182   0.047669   0.948   0.3545  
#RTI_resnorm_diff  0.007244   0.031217   0.232   0.8189  
#RVP_resnorm_diff  0.089276   0.063271   1.411   0.1736  
#SOC_resnorm_diff -0.110717   0.052267  -2.118   0.0469 *
#SSP_resnorm_diff  0.097418   0.055834   1.745   0.0964 .
etaSquared(lm_G_cantab_diff, anova=F, type=3)
#                      eta.sq eta.sq.part
#RTI_resnorm_diff 0.001855949 0.002685059
#RVP_resnorm_diff 0.068622495 0.090533394
#SOC_resnorm_diff 0.154663435 0.183245986
#SSP_resnorm_diff 0.104929685 0.132105486
rm(tmp)


# Pre GTA
lm_G_GTA_t1 = lm(G_t1 ~ SC_Q_resnorm_t1 + SC_Eglob_resnorm_t1 + SC_PC_resnorm_t1)
shapiro.test(residuals(lm_G_GTA_t1)) #p=0.76
summary(lm_G_GTA_t1)
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.695151   0.046379  36.550   <2e-16 ***
#SC_Q_resnorm_t1      0.046394   0.059857   0.775   0.4462    
#SC_Eglob_resnorm_t1 -0.148612   0.058513  -2.540   0.0183 *  
#SC_PC_resnorm_t1    -0.007956   0.045318  -0.176   0.8622   
etaSquared(lm_G_GTA_t1, anova=F, type=3)
#                          eta.sq eta.sq.part
#SC_Q_resnorm_t1     0.0179563474 0.025454873
#SC_Eglob_resnorm_t1 0.1928110513 0.219035374
#SC_PC_resnorm_t1    0.0009213385 0.001338408
rm(lm_G_GTA_t1)

# Post GTA
lm_G_GTA_t2 = lm(G_t2 ~ SC_Q_resnorm_t2 + SC_Eglob_resnorm_t2 + SC_PC_resnorm_t2)
shapiro.test(residuals(lm_G_GTA_t2)) #p=0.0082
hist(residuals(lm_G_GTA_t2))
summary(lm_G_GTA_t2)
#                     Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          1.735281   0.056627  30.644   <2e-16 ***
#SC_Q_resnorm_t2     -0.017297   0.048176  -0.359   0.7228    
#SC_Eglob_resnorm_t2 -0.092214   0.033249  -2.773   0.0108 *  
#SC_PC_resnorm_t2    -0.007215   0.037514  -0.192   0.8492  
etaSquared(lm_G_GTA_t2, anova=F, type=3)
#                           eta.sq eta.sq.part
#SC_Q_resnorm_t2     0.003862355 0.005573713
#SC_Eglob_resnorm_t2 0.230458838 0.250619800
#SC_PC_resnorm_t2    0.001108296 0.001605749
rm(lm_G_GTA_t2)

# Post GTA without Eglob outlier
tmp=results
tmp$SC_Eglob_resnorm_t2[18]<-NA
lm_G_GTA_t2_v2=lm(tmp$G_t2~tmp$SC_Q_resnorm_t2+tmp$SC_Eglob_resnorm_t2+tmp$SC_PC_resnorm_t2)
summary(lm_G_GTA_t2_v2)
#                        Estimate Std. Error t value Pr(>|t|)    
#(Intercept)              1.72422    0.05767  29.896   <2e-16 ***
#tmp$SC_Q_resnorm_t2     -0.02735    0.04919  -0.556    0.584    
#tmp$SC_Eglob_resnorm_t2 -0.07801    0.03612  -2.160    0.042 *  
#tmp$SC_PC_resnorm_t2    -0.01882    0.03924  -0.479    0.636 
etaSquared(lm_G_GTA_t2_v2,anova=F,type=3)
#                             eta.sq eta.sq.part
#tmp$SC_Q_resnorm_t2     0.011544664  0.01385884
#tmp$SC_Eglob_resnorm_t2 0.174149554  0.17491519
#tmp$SC_PC_resnorm_t2    0.008584631  0.01034221
rm(lm_G_GTA_t2_v2)

# Difference GTA
lm_G_GTA_diff = lm(G_diff ~ SC_Q_resnorm_diff + SC_Eglob_resnorm_diff + SC_PC_resnorm_diff)
shapiro.test(residuals(lm_G_GTA_diff)) #p=0.86
summary(lm_G_GTA_diff)
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)  
#(Intercept)            0.06222    0.04732   1.315   0.2015  
#SC_Q_resnorm_diff     -0.00579    0.04702  -0.123   0.9031  
#SC_Eglob_resnorm_diff -0.06903    0.02942  -2.347   0.0279 *
#SC_PC_resnorm_diff    -0.02100    0.02523  -0.832   0.4138 
etaSquared(lm_G_GTA_diff, anova=F, type=3)
#                           eta.sq eta.sq.part
#SC_Q_resnorm_diff     0.0005204035 0.0006587128
#SC_Eglob_resnorm_diff 0.1890536601 0.1931949439
#SC_PC_resnorm_diff    0.0237830260 0.0292428544
rm(lm_G_GTA_diff)

# Diff GTA without Eglob outlier
tmp$SC_Eglob_resnorm_diff[18]<-NA
lm_G_GTA_dif_v2=lm(tmp$G_diff~tmp$SC_Q_resnorm_diff+tmp$SC_Eglob_resnorm_diff+tmp$SC_PC_resnorm_diff)
summary(lm_G_GTA_dif_v2)
#                           Estimate Std. Error t value Pr(>|t|)
#(Intercept)                0.041280   0.046966   0.879    0.389
#tmp$SC_Q_resnorm_diff     -0.009642   0.045150  -0.214    0.833
#tmp$SC_Eglob_resnorm_diff -0.033172   0.034980  -0.948    0.353
#tmp$SC_PC_resnorm_diff    -0.037270   0.025949  -1.436    0.165
etaSquared(lm_G_GTA_dif_v2, anova=F, type=3)
#                              eta.sq eta.sq.part
#tmp$SC_Q_resnorm_diff     0.001800258 0.002068638
#tmp$SC_Eglob_resnorm_diff 0.035499273 0.039270791
#tmp$SC_PC_resnorm_diff    0.081434204 0.085729584
rm(lm_G_GTA_dif_v2); rm(tmp)



##--- Regression models for JiT (in patients only!) 

# Pre CANTAB
lm_JiT_cantab_t1 = lm(JiT_res_t1[group!='CON'] ~ RTI_resnorm_t1[group!='CON'] + RVP_resnorm_t1[group!='CON'] + 
                        SOC_resnorm_t1[group!='CON'] + SSP_resnorm_t1[group!='CON'])
summary(lm_JiT_cantab_t1)
#                             Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                    -0.17097    0.16792  -1.018   0.3305  
#RTI_resnorm_t1[group != "CON"] -0.60852    0.23008  -2.645   0.0228 *
#RVP_resnorm_t1[group != "CON"] -0.34851    0.24675  -1.412   0.1855  
#SOC_resnorm_t1[group != "CON"] -0.03656    0.19628  -0.186   0.8556  
#SSP_resnorm_t1[group != "CON"]  0.04190    0.26138   0.160   0.8756  
rm(lm_JiT_cantab_t1)


# Post CANTAB
lm_JiT_cantab_t2 = lm(JiT_res_t2[group!='CON'] ~ RTI_resnorm_t2[group!='CON'] + RVP_resnorm_t2[group!='CON'] + 
                        SOC_resnorm_t2[group!='CON'] + SSP_resnorm_t2[group!='CON'])
summary(lm_JiT_cantab_t2)
# res                            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                     -0.1931     0.1575  -1.226   0.2437  
#RTI_resnorm_t2[group != "CON"]  -0.3480     0.1217  -2.859   0.0144 *
#RVP_resnorm_t2[group != "CON"]  -0.3068     0.2092  -1.467   0.1682  
#SOC_resnorm_t2[group != "CON"]  -0.1133     0.1535  -0.738   0.4746  
#SSP_resnorm_t2[group != "CON"]  -0.1053     0.2343  -0.449   0.6613
etaSquared(lm_JiT_cantab_t2, anova=F, type=3)
#RTI_resnorm_t2[group != "CON"] 0.394228158  0.40509071
#RVP_resnorm_t2[group != "CON"] 0.103797662  0.15202782
#SOC_resnorm_t2[group != "CON"] 0.026289165  0.04343551
#SSP_resnorm_t2[group != "CON"] 0.009735887  0.01653815
rm(lm_JiT_cantab_t2)


# Difference CANTAB
JiT_res_diff=JiT_res_t2-JiT_res_t1
lm_JiT_cantab_diff = lm(JiT_res_diff[group!='CON'] ~ RTI_resnorm_diff[group!='CON'] + 
                          RVP_resnorm_diff[group!='CON'] + SOC_resnorm_diff[group!='CON'] + 
                          SSP_resnorm_diff[group!='CON'])
summary(lm_JiT_cantab_diff)
#(Intercept)                       0.02497    0.05595   0.446    0.664
#RTI_resnorm_diff[group != "CON"]  0.03058    0.04563   0.670    0.517
#RVP_resnorm_diff[group != "CON"]  0.05113    0.06713   0.762    0.462
#SOC_resnorm_diff[group != "CON"] -0.01769    0.04182  -0.423    0.680
#SSP_resnorm_diff[group != "CON"]  0.03293    0.06959   0.473    0.645
rm(lm_JiT_cantab_diff)


# Pre GTA
lm_JiT_GTA_t1 = lm(JiT_res_t1[group!='CON'] ~ SC_Q_resnorm_t1[group!='CON'] + SC_Eglob_resnorm_t1[group!='CON'] + 
                     SC_PC_resnorm_t1[group!='CON'])
summary(lm_JiT_GTA_t1)
#(Intercept)                         -0.05173    0.27014  -0.191    0.851
#SC_Q_resnorm_t1[group != "CON"]      0.11410    0.27773   0.411    0.688
#SC_Eglob_resnorm_t1[group != "CON"] -0.14187    0.40600  -0.349    0.732
#SC_PC_resnorm_t1[group != "CON"]    -0.07441    0.21274  -0.350    0.732
rm(lm_JiT_GTA_t1)


# Post GTA
lm_JiT_GTA_t2 = lm(JiT_res_t2[group!='CON'] ~ SC_Q_resnorm_t2[group!='CON'] + SC_Eglob_resnorm_t2[group!='CON'] + 
                     SC_PC_resnorm_t2[group!='CON'])
summary(lm_JiT_GTA_t2)
#(Intercept)                         -0.07004    0.19090  -0.367    0.720
#SC_Q_resnorm_t2[group != "CON"]     -0.04484    0.14517  -0.309    0.762
#SC_Eglob_resnorm_t2[group != "CON"]  0.08596    0.11313   0.760    0.461
#SC_PC_resnorm_t2[group != "CON"]    -0.11132    0.12328  -0.903    0.383
rm(lm_JiT_GTA_t2)


# Difference GTA
lm_JiT_GTA_diff = lm(JiT_res_diff[group!='CON'] ~ SC_Q_resnorm_diff[group!='CON'] + 
                       SC_Eglob_resnorm_diff[group!='CON'] + SC_PC_resnorm_diff[group!='CON'])
summary(lm_JiT_GTA_diff)
#                                      Estimate Std. Error t value Pr(>|t|)
#(Intercept)                           0.054232   0.040791   1.330    0.207
#SC_Q_resnorm_diff[group != "CON"]     0.021467   0.036634   0.586    0.568
#SC_Eglob_resnorm_diff[group != "CON"] 0.004445   0.026764   0.166    0.871
#SC_PC_resnorm_diff[group != "CON"]    0.036360   0.021192   1.716    0.110
rm(lm_JiT_GTA_diff)



###--- Regression models for JiNT (PAT + CON)

# Pre CANTAB
lm_JiNT_cantab_t1 = lm(JiNT_res_t1 ~ RTI_resnorm_t1 + RVP_resnorm_t1 + SOC_resnorm_t1 + SSP_resnorm_t1)
shapiro.test(residuals(lm_JiNT_cantab_t1)) #p=0.37
summary(lm_JiNT_cantab_t1)
# res            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.565610   0.008987  62.939   <2e-16 ***
#RTI_resnorm_t1 -0.007391   0.011610  -0.637   0.5313    
#RVP_resnorm_t1  0.001770   0.013346   0.133   0.8957    
#SOC_resnorm_t1  0.013494   0.009936   1.358   0.1889    
#SSP_resnorm_t1 -0.028263   0.012654  -2.234   0.0365 *
rm(lm_JiNT_cantab_t1)

# Post CANTAB
lm_JiNT_cantab_t2 = lm(JiNT_res_t2 ~ RTI_resnorm_t2 + RVP_resnorm_t2 + SOC_resnorm_t2 + SSP_resnorm_t2)
summary(lm_JiNT_cantab_t2)
#(Intercept)     0.592778   0.010572  56.072   <2e-16 ***
#RTI_resnorm_t2 -0.005910   0.006766  -0.873    0.392    
#RVP_resnorm_t2 -0.002361   0.016295  -0.145    0.886    
#SOC_resnorm_t2 -0.008917   0.010977  -0.812    0.425    
#SSP_resnorm_t2 -0.007968   0.016159  -0.493    0.627  
rm(lm_JiNT_cantab_t2)


# Difference CANTAB
lm_JiNT_cantab_diff = lm(JiNT_res_diff ~ RTI_resnorm_diff + RVP_resnorm_diff + SOC_resnorm_diff + SSP_resnorm_diff)
summary(lm_JiNT_cantab_diff)
#                  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)       0.027621   0.009261   2.983   0.0071 **
#RTI_resnorm_diff -0.009352   0.006169  -1.516   0.1444   
#RVP_resnorm_diff  0.001100   0.012523   0.088   0.9309   
#SOC_resnorm_diff  0.002832   0.008686   0.326   0.7477   
#SSP_resnorm_diff -0.030164   0.011076  -2.723   0.0127 * 
etaSquared(lm_JiNT_cantab_diff)
#                       eta.sq  eta.sq.part
#RTI_resnorm_diff 0.0738231482 0.0986468715
#RVP_resnorm_diff 0.0002476933 0.0003670715
#SOC_resnorm_diff 0.0034136571 0.0050352770
#SSP_resnorm_diff 0.2382188265 0.2609892506
rm(lm_JiNT_cantab_diff)


# Pre GTA
lm_JiNT_GTA_t1 = lm(JiNT_res_t1 ~ SC_Q_resnorm_t1 + SC_Eglob_resnorm_t1 + SC_PC_resnorm_t1)
summary(lm_JiNT_GTA_t1)
# (Intercept)         0.562693   0.010757  52.311   <2e-16 ***
#SC_Q_resnorm_t1     0.003532   0.013883   0.254    0.801    
#SC_Eglob_resnorm_t1 0.010054   0.013571   0.741    0.466    
#SC_PC_resnorm_t1    0.005700   0.010511   0.542    0.593 
rm(lm_JiNT_GTA_t1)

# Post GTA
lm_JiNT_GTA_t2 = lm(JiNT_res_t2 ~ SC_Q_resnorm_t2 + SC_Eglob_resnorm_t2 + SC_PC_resnorm_t2)
summary(lm_JiNT_GTA_t2)
# res                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          0.593163   0.010513  56.422   <2e-16 ***
#SC_Q_resnorm_t2     -0.017426   0.008944  -1.948   0.0637 .  
#SC_Eglob_resnorm_t2 -0.004347   0.006173  -0.704   0.4883    
#SC_PC_resnorm_t2    -0.004581   0.006965  -0.658   0.5172  
etaSquared(lm_JiNT_GTA_t2,anova=FALSE,type=3)
#                        eta.sq eta.sq.part
#SC_Q_resnorm_t2     0.13946952  0.14167279
#SC_Eglob_resnorm_t2 0.01822092  0.02110858
#SC_PC_resnorm_t2    0.01589470  0.01846346
rm(lm_JiNT_GTA_t2)


# Difference GTA
lm_JiNT_GTA_diff = lm(JiNT_res_diff ~ SC_Q_resnorm_diff + SC_Eglob_resnorm_diff + SC_PC_resnorm_diff)
summary(lm_JiNT_GTA_diff)
#(Intercept)            0.027791   0.009615   2.890  0.00826 **
#SC_Q_resnorm_diff     -0.010475   0.009554  -1.096  0.28424   
#SC_Eglob_resnorm_diff  0.007561   0.005977   1.265  0.21851   
#SC_PC_resnorm_diff    -0.001015   0.005125  -0.198  0.84473  
rm(lm_JiNT_GTA_diff)



### Plot significant associations --------------------------------------------------------------------------------#
# when using JiRes

setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_regression")

results$group=factor(results$group, ordered=TRUE, levels=c('MEN', 'GLI', 'CON'))
detach(results); attach(results)

#--> Plot all together!

# 1) G post ~ Eglob post (all subjects)
p1<-ggplot(results, aes(x=SC_Eglob_resnorm_t2, y=G_t2)) +
  geom_smooth(method=lm, col="darkgray") + 
  geom_point(aes(fill=group, shape=group),size=6) +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25), legend.position="none") + 
  labs(y = "G (post)", x = "Global efficiency (post)") +
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))
  

# 2) G diff ~ SOC diff (all subjects)
p2<-ggplot(results, aes(x=SOC_resnorm_diff, y=G_diff)) + 
  geom_smooth(method=lm, col="darkgray") + 
  geom_point(aes(fill=group, shape=group), size=6) +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25), legend.position="none") + 
  labs(y = "G (post-pre)", x = "Planning accuracy (post-pre)") +
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))


# 3) JiT ~ RTI post (patients only)
results_PAT=results[group!='CON',]
detach(results); attach(results_PAT)

p3<-ggplot(results_PAT, aes(x=RTI_resnorm_t2, y=JiT_t2)) +
  geom_smooth(method=lm, col="darkgray") + 
  geom_point(aes(fill=group, shape=group), size=6) +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25), legend.position="none") + 
  labs(y = expression('J'[tumor]*' (post)'), x = "Reaction time (post)") +
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))


# 4) JiNT ~ SSP diff (patients + controls)
detach(results_PAT); attach(results)
p4<-ggplot(results, aes(x=SSP_resnorm_diff, y=JiNT_res_diff)) + 
  geom_smooth(method=lm, col="darkgray") + 
  geom_point(aes(fill=group,shape=group), size=6) +
  theme(axis.ticks = element_blank(),
        text = element_text(size = 25), legend.position="none") + 
  labs(y = expression('J'[non-tumor]*' (post-pre)'), x = "Working memory (post-pre)") +
  scale_fill_brewer(palette="Paired") +
  scale_shape_manual(values=c(21,24,22))

png('results_regression_all.png', width=1200, height=800)
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()
rm(list=c('p1','p2','p3','p4','results_PAT'))


