
##########################################################################################
###### Stat test comparing proportion of proteins with hits to the COG database between cassettes and metagenomes ###### 

###Check normality assumption for Two-sample T-test
with(COG_freq, shapiro.test(Percent[Type == "Cassettes"])) #W = 0.94096, p-value = 0.5106 - Normal distribution
with(COG_freq, shapiro.test(Percent[Type == "Metagenomes"])) #W = 0.88707, p-value = 0.108 - Normal distribution

###Check homogeniety of variance assumption for Two-sample T-test
res.ftest <- var.test(Percent ~ Type, data = COG_freq)
res.ftest #p-value = 0.03786 - Does not meet assumption

###Use non-parametric Two-sample Wilcoxon test
res <- wilcox.test(Percent ~ Type, data = COG_freq, exact = FALSE, alternative = "less")
res #W = 0, p-value = 1.829e-05


##########################################################################################
#Chi-squared tests comparing distribution of COG categories between cassettes and metagenomes ######  
chi_HE10<-chisq.test(HE10)
chi_HE29<-chisq.test(HE29)
chi_O150<-chisq.test(O150)
chi_X1634<-chisq.test(X1634)
chi_MQ06<-chisq.test(MQ06)
chi_MQ12<-chisq.test(MQ12)
chi_MCB.C<-chisq.test(MCB.C)
chi_MCB.D<-chisq.test(MCB.D)
chi_LCS2.3<-chisq.test(LCS2.3)
chi_LCS3.1<-chisq.test(LCS3.1)
chi_AS.1<-chisq.test(AS.1)
chi_BM.3<-chisq.test(BM.3)

### Extract Pearson's residuals calculated from Chi-squared contingency tables
write.table(data.frame(chi_HE10$residual), file = "HE10.residuals.txt")
write.table(data.frame(chi_HE29$residual), file = "HE29.residuals.txt")
write.table(data.frame(chi_O150$residual), file = "O150.residuals.txt")
write.table(data.frame(chi_X1634$residual), file = "X1634.residuals.txt")
write.table(data.frame(chi_MQ06$residual), file = "MQ06.residuals.txt")
write.table(data.frame(chi_MQ12$residual), file = "MQ12.residuals.txt")
write.table(data.frame(chi_MCB.C$residual), file = "MCB-C.residuals.txt")
write.table(data.frame(chi_MCB.D$residual), file = "MCB-D.residuals.txt")
write.table(data.frame(chi_LCS2.3$residual), file = "LCS2.3.residuals.txt")
write.table(data.frame(chi_LCS3.1$residual), file = "LCS3.1.residuals.txt")
write.table(data.frame(chi_AS.1$residual), file = "AS-1.residuals.txt")
write.table(data.frame(chi_BM.3$residual), file = "BM-3.residuals.txt")



#############################################################################################
#Chi-squared tests comparing distribution of the COG categories 
#('Defence'- V; & Secondary metabolite biosynthesis, transport and catabolism' - Q)
#between cassettes and metagenomes ######  

chi_HE10<-chisq.test(HE10)
chi_HE29<-chisq.test(HE29)
chi_O150<-chisq.test(O150)
chi_X1634<-chisq.test(X1634)
chi_MQ06<-chisq.test(MQ06)
chi_MQ12<-chisq.test(MQ12)
chi_MCB.C<-chisq.test(MCB.C)
chi_MCB.D<-chisq.test(MCB.D)
chi_LCS2.3<-chisq.test(LCS2.3)
chi_LCS3.1<-chisq.test(LCS3.1)
chi_AS.1<-chisq.test(AS.1)
chi_BM.3<-chisq.test(BM.3)

### Extract Pearson's residuals calculated from Chi-squared contingency tables
write.table(data.frame(chi_HE10$residual), file = "HE10.residuals.txt")
write.table(data.frame(chi_HE29$residual), file = "HE29.residuals.txt")
write.table(data.frame(chi_O150$residual), file = "O150.residuals.txt")
write.table(data.frame(chi_X1634$residual), file = "X1634.residuals.txt")
write.table(data.frame(chi_MQ06$residual), file = "MQ06.residuals.txt")
write.table(data.frame(chi_MQ12$residual), file = "MQ12.residuals.txt")
write.table(data.frame(chi_MCB.C$residual), file = "MCB-C.residuals.txt")
write.table(data.frame(chi_MCB.D$residual), file = "MCB-D.residuals.txt")
write.table(data.frame(chi_LCS2.3$residual), file = "LCS2.3.residuals.txt")
write.table(data.frame(chi_LCS3.1$residual), file = "LCS3.1.residuals.txt")
write.table(data.frame(chi_AS.1$residual), file = "AS-1.residuals.txt")
write.table(data.frame(chi_BM.3$residual), file = "BM-3.residuals.txt")


####################################################################################################################################################################
###### Stat test comparing SignalP freq between cassettes and metagenomes ###### 

###Check normality assumption for Two-sample T-test
with(sp_freq, shapiro.test(Percent[Type == "Cassettes"])) #W = 0.8793, p-value = 0.08585 - Normal distribution
with(sp_freq, shapiro.test(Percent[Type == "Metagenomes"])) #W = 0.97611, p-value = 0.9632 - Normal distribution

###Check homogeniety of variance assumption for Two-sample T-test
res.ftest <- var.test(Percent ~ Type, data = sp_freq)
res.ftest #p-value = 0.006063 - DOes not meet assumption

###Use non-parametric Two-sample Wilcoxon test
res <- wilcox.test(Percent ~ Type, data = sp_freq, exact = FALSE, alternative = "greater")
res #W = 144, p-value = 1.829e-05

####################################################################################################################################################################

#### Stats for proportion of AMR genes cassettes vs meatgenomes
prop.test(x = c(54,25), n = c(4618,21963), alternative = 'greater') #AS-1
prop.test(x = c(32,2), n = c(3035,7157), alternative = 'greater') #BM-3
prop.test(x = c(28,164), n = c(5817,217132), alternative = 'greater') #HE10
prop.test(x = c(41,103), n = c(6291,152853), alternative = 'greater') #HE29
prop.test(x = c(13,11), n = c(1604,8183), alternative = 'greater') #LCS2.3
prop.test(x = c(9,28), n = c(1193,24207), alternative = 'greater') #LCS3.1
prop.test(x = c(27,69), n = c(3457,78924), alternative = 'greater') #MCB-C
prop.test(x = c(33,10), n = c(2357,4867), alternative = 'greater') #MCB-D
prop.test(x = c(30,50), n = c(6989,70759), alternative = 'greater') #MQ06
prop.test(x = c(8,37), n = c(2157,56492), alternative = 'greater') #MQ12
prop.test(x = c(5,69), n = c(3240,93832), alternative = 'greater') #O150
prop.test(x = c(18,60), n = c(4611,81186), alternative = 'greater') #X1634

##########################################################################################
##### AMR gene cassette accumulation curve using Vegan ##### 

library(vegan)

sp1<- specaccum(AMR_abund)

## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, model="lomolino")
coef(mod1)
fitted(mod1)
plot(sp1)
plot(mod1, add = TRUE, col=2, lwd=2)

