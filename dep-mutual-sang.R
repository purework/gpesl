library(limma)
library(data.table)
library(Rcpp)
FDR=0.2
openmp = 1
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp("ISLE/src/essentialityTestPair.cpp", rebuild=T)
args = commandArgs()
f = args[4]
of = paste0(f, '-dep-mutual-new2')
load("dp.crispr.sanger.rda")
gt = alias2SymbolTable(dp.crispr.sanger$genes) #
sl = fread(f, head=T)
sl$g1 = match(alias2SymbolTable(sl$g1), gt)
sl$g2 = match(alias2SymbolTable(sl$g2), gt)
sl = sl[!is.na(sl$g1) & !is.na(sl$g2)]
sl = data.matrix(sl)
pcl=NULL
for (k in seq(1)){
	if (k==1) probc=dp.crispr.sanger
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
