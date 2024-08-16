#Load bsseq and methylSig packages

library(bsseq)
library(methylSig)

#Place all Male .bedGraph (pipeline A) or .cov (pipeline B) to be analyzed in a single working directory (MethylMWD)
setwd('MethylMWD')
files <- list.files(pattern = ".bedGraph")
#If pipeline B change to
#files <- list.files(pattern = ".cov")

#Create a phenotype file (PM3) for Males
Temp1<-c("H","H","H","H","H","H","H")
Temp2<-c("L","L","L","L","L","L","L")
Sex1<-c("M","M","M","M","M","M","M")
P1<-cbind.data.frame(Temp1,Sex1)
P2<-cbind.data.frame(Temp2,Sex1)
colnames(P2)<-c("Temp1", "Sex1")
PM3<-rbind(P1,P2)
rownames(PM3)<-c("H1M1", "H1M2","H1M3", "H1M4","H1M5","H1M6","H1M7",
"L1M1","L1M2","L1M3","L1M4","L1M5","L1M6","L1M7")
PM3


#Load in files to a bsseq object, MethylDackyl files should be strandCollapse=FALSE
bs1 = bsseq::read.bismark(
    files = files,
    colData = PM3,
    rmZeroCov = FALSE,
    strandCollapse = FALSE)

#Filter by min and max coverage
bs2 =filter_loci_by_coverage(bs1, min_count = 10, max_count = 200)

#Retain only sites with sufficient coverage in at least 4 individuals per group
bs3 = filter_loci_by_group_coverage(
    bs = bs2,
    group_column = 'Temp1',
    c('H' = 4, 'L' = 4))

#Tile methylation across 100 bp windows
bs4 = tile_by_windows(bs = bs3, win_size = 100)

#Retain only windows with sufficient coverage in at least 4 individuals per group
bsfin = filter_loci_by_group_coverage(
    bs = bs4,
    group_column = 'Temp1',
    c('H' = 4, 'L' = 4))

#Check number of windows with GenomicRanges
GenomicRanges::granges(bsfin)

#Create a simple DSS fit model in methylSig
diff_fit_simple = diff_dss_fit(
    bs = bsfin,
    design = bsseq::pData(bsfin),
    formula = as.formula('~ Temp1'))

#Set a simple contrast
simple_contrast = matrix(c(0,1), ncol = 1)

#Compare methylation among temperature groups, with H set to case
diff_simple_M = diff_dss_test(
    bs = bsfin,
    diff_fit = diff_fit_simple,
    contrast = simple_contrast,
    methylation_group_column = 'Temp1',
    methylation_groups = c('case' = 'H', 'control' = 'L'))
write.csv(diff_simple_M, file="MethylOutputMales.csv")

#Remove those windows with 0 difference in methylation AND 0 methylation in both treatments (case/control)
diff_simple_M1<-as.data.frame(diff_simple_M)
diff_simple_M1$absMD<-abs(diff_simple_M1$meth_diff)
diff_simple_M1$All<-diff_simple_M1$meth_case+diff_simple_M1$meth_control+diff_simple_M1$absMD
diff_simple_MRem<-subset(diff_simple_M1,diff_simple_M1$All!=0)

#Perform fdr adjustment after removal of windows with 0 methylation
fdr2<-p.adjust(diff_simple_MRem$pvalue, method = "fdr")
diff_simple_MRemFin<-as.data.frame(cbind(diff_simple_MRem,fdr2))
write.csv(diff_simple_MRemFin,file="Male4_7_zerorem_fdr.csv", row.names=FALSE)

diff_simple_MRemFin4<-subset(diff_simple_MRemFin,diff_simple_MRemFin$fdr2<0.05)
write.csv(diff_simple_MRemFin4,file="Male4_7_zerorem_fdr_05.csv", row.names=FALSE)


#Place all Female bedGraphs .bedGraph (pipeline A) or .cov (pipeline B) to be analyzed in a single working directory (MethylFWD)
setwd('MethylFWD')
files <- list.files(pattern = ".bedGraph")
#If pipeline B, change to: 
#files <- list.files(pattern = ".cov")

#Create a phenotype file (PF3) for Females
Temp1<-c("H","H","H","H","H","H","H")
Temp2<-c("L","L","L","L","L")
Sex3<-c("F","F","F","F","F","F","F")
Sex2<-c("F","F","F","F","F")
P4<-cbind.data.frame(Temp1,Sex3)
P5<-cbind.data.frame(Temp2,Sex2)
colnames(P5)<-c("Temp1", "Sex3")
PF3<-rbind(P4,P5)
rownames(PF3)<-c("H1F1", "H1F2","H1F3", "H1F4","H1F5","H1F6","H1F7",
"L1F1","L1F2","L1F3","L1F4","L1F5","L1F6","L1F7")
PF3

#Repeatsteps
bsf1 = bsseq::read.bismark(
    files = files,
    colData = PF3,
    rmZeroCov = FALSE,
    strandCollapse = FALSE)

#Filter by min and max coverage
bsf2 =filter_loci_by_coverage(bsf1, min_count = 10, max_count = 200)

#Retain only sites with sufficient coverage in at least 4 individuals per group
bsf3 = filter_loci_by_group_coverage(
    bs = bsf2,
    group_column = 'Temp1',
    c('H' = 4, 'L' = 4))

#Tile methylation across 100 bp windows
bsf4 = tile_by_windows(bs = bsf3, win_size = 100)

#Retain only windows with sufficient coverage in at least 4 individuals per group
bsfinfem = filter_loci_by_group_coverage(
    bs = bsf4,
    group_column = 'Temp1',
    c('H' = 4, 'L' = 4))

#Check number of windows with GenomicRanges
GenomicRanges::granges(bsfinfem)

#Create a simple DSS fit model in methylSig
diff_fit_simpleF = diff_dss_fit(
    bs = bsfinfem,
    design = bsseq::pData(bsfinfem),
    formula = as.formula('~ Temp1'))

#Set a simple contrast
simple_contrast = matrix(c(0,1), ncol = 1)

#Compare methylation among tempearture groups, with H set to case
diff_simple_F = diff_dss_test(
    bs = bsfinfem,
    diff_fit = diff_fit_simpleF,
    contrast = simple_contrast,
    methylation_group_column = 'Temp1',
    methylation_groups = c('case' = 'H', 'control' = 'L'))
write.csv(diff_simple_F, file="MethylOutputFemales.csv")

#Remove those windows with 0 difference in methylation AND 0 methylation in both treatments (case/control)
diff_simple_F1<-as.data.frame(diff_simple_F)
diff_simple_F1$absMD<-abs(diff_simple_F1$meth_diff)
diff_simple_F1$All<-diff_simple_F1$meth_case+diff_simple_F1$meth_control+diff_simple_F1$absMD
diff_simple_FRem<-subset(diff_simple_F1,diff_simple_F1$All!=0)

#Perform fdr adjustment after removal of windows with 0 methylation
fdr3<-p.adjust(diff_simple_FRem$pvalue, method = "fdr")
diff_simple_FRemFin<-as.data.frame(cbind(diff_simple_FRem,fdr3))
write.csv(diff_simple_FRemFin,file="Female4_7_zerorem_fdr.csv", row.names=FALSE)

diff_simple_FRemFin4<-subset(diff_simple_FRemFin,diff_simple_FRemFin$fdr2<0.05)
write.csv(diff_simple_FRemFin4,file="Female4_7_zerorem_fdr_05.csv", row.names=FALSE)
