library(data.table)
library(Rcpp)
library(limma)
FDR=0.2
openmp = 1
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp("ISLE/src/essentialityTestPair.cpp", rebuild=T)
args = commandArgs()
f = args[4]
of = paste0(f, '-dep-mutual')
load('isle-gene-sample-table.rda')
gt = alias2SymbolTable(prob$genes)
sl = fread(f, head=T)
sl$g1 = match(sl$g1, gt)
sl$g2 = match(sl$g2, gt)
sl = sl[!is.na(sl$g1) & !is.na(sl$g2)]
sl = data.matrix(sl)
load("ISLE/data/prob.ach.old1.RData") # Cheung et al. PNAS (2011)
load("ISLE/data/prob.ach.new1.RData") # Cowley et al. Sci. Data. (2014)
load("ISLE/data/prob.ach.new2.RData") # Auirre et al. Cancer Discov (2016)
load("ISLE/data/prob.mar.old1.RData") # Marcotte et al. Cancer Discov (2012)
load("ISLE/data/prob.mar.new1.RData") # Marcotte et al. Cell (2016)
pcl=NULL
for (k in seq(5)){
	if (k==1) probc=prob.ach.old1		
	if (k==2) probc=prob.ach.new1
	if (k==3) probc=prob.ach.new2
	if (k==4) probc=prob.mar.old1
	if (k==5) probc=prob.mar.new1
	mRNAq=probc$mRNAq2
	scnaq=probc$scnaq2
	mat=probc$mat
	p = essentialityTestPair(sl, t(scnaq), t(mRNAq), t(mat))
	pcl=cbind(pcl,p)
}
pcl[pcl< -999] = NA
pcl[pcl==0.5] = NA
hd = colnames(sl)
sl = data.table(sl, pcl)
colnames(sl) = c(hd, paste0('V', 1:ncol(pcl)))
sl$g1 = gt[sl$g1]
sl$g2 = gt[sl$g2]
write.table(sl, file=of, quote=F, row.names=F, sep='\t')
