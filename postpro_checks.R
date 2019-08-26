###################################################################################################################
###                                                                                                             ###
###                   MODELING BRAIN DYNAMICS AFTER TUMOR RESECTION USING THE VIRTUAL BRAIN                     ###
###               ==============================================================================                ###
###                                           PART 3: Additional checks                                         ###
###                                                                                                             ###
### Created by Hannelore Aerts                                                                                  ###
### Date last update: 18/10/2018                                                                                ###
###################################################################################################################

## Some meningioma patients have very low JiT values. These patients have in common that they have rather small,
## localized tumors. Given that Ji distribution is skewed (majority rather low, some very high values) these 
## very high values could maybe be some kind of outlier that could not be diminished because of small number
## of ROIs?
## --> check relation between number of tumornodes and JiT!

setwd("/home/hannelore/Documents/ANALYSES/TVB_post")
tumorrois=read.table(file="tumorrois_Min1voxelTumorOverlap_scale68.csv", header=F, sep=",")
results=read.table(file="RESULTS_ALL_afterprep_20190411.csv", header=T, sep=";")
attach(results)
results_subset=results[group!='CON',]
results_subset$group = factor(results_subset$group, ordered=TRUE, levels=c('MEN', 'GLI'))
tumorrois_subset=tumorrois[group!='CON']

detach(results); attach(results_subset)

n_tumorrois=colSums(tumorrois_subset>0)

palette(c("cyan", "magenta"))
setwd("/home/hannelore/Documents/ANALYSES/TVB_post/results_TVBii")
png('JiT_NumberOfTumorrois.png')
plot(n_tumorrois, JiT_res_t2, pch=21, bg=group, cex=1, xlab="Number of tumorregions", ylab=expression("J"[tumor]))#,
#     cex.axis=2, cex.lab=3)
legend(16, -1, levels(group), xpd=TRUE, pt.bg=1:length(group), box.lty=0, pch=21)#, cex=3)
dev.off()
