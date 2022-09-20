library(SummarizedExperiment)
library(magrittr)
library(data.table)
library(dplyr)
library(limma)
library(igraph)
library(rhdf5)
library(bigstatsr)
library(BiocParallel)
cpu = 8
register(MulticoreParam(workers=cpu), default=T)
count2tpm = function(mat, bp){
  a = assay(mat) / bp * 1e3
  b = colSums(a) / 1e6
  assay(mat) = t(t(a) / b)
  rm(a, b)
  gc()
  mat
}
gv25 = fread('gencode.v25.annotation.gtf-enstype.txt', head=F)
colnames(gv25) = c('ensgene', 'symbol', 'biotype')
setkey(gv25, ensgene)
load('rc/gtex-rse_gene.Rdata')
gtex = rse_gene
load('rc/tcga-rse_gene.Rdata')
tcga = rse_gene
bp = as.data.frame(rowData(gtex))$bp_length
gtex = count2tpm(gtex, bp)
tcga = count2tpm(tcga, bp)
idx = (rownames(rse_gene) %in% gv25$ensgene[gv25$biotype == 'protein_coding'])
gtex = gtex[idx,]
tcga = tcga[idx,]
tcga.meta = fread('rc/tcga-phenotype.tsv', sep='\t')
tcga.meta = tcga.meta[, c('gdc_file_id', 'gdc_cases.samples.submitter_id')]
tmp = strsplit(tcga.meta$gdc_cases.samples.submitter_id, '-')
tcga.meta$code = sapply(tmp, function(s) as.integer(substr(s[4], 1, 2)))
tcga.meta$gdc_file_id = toupper(tcga.meta$gdc_file_id)
tcga.meta.tumor = tcga.meta[code < 10]
tcga.tumor = tcga[,tcga.meta.tumor$gdc_file_id]
# gtex
nr = nrow(gtex)
nc = ncol(gtex)
step = ceiling(nc / cpu)
mat = FBM(nr, nr, init=0)
reslist = bplapply(seq(1, nc, step), function(x){
  for(i in x:min(x+step-1, nc)){
    gid = which(assay(gtex[,i]) <= 0.5)
    mat[gid, gid] = 1
  }
  NULL
})
mat = mat[]
diag(mat) = 0
mat = 1 - mat
# tcga
nr = nrow(tcga.tumor)
nc = ncol(tcga.tumor)
step = ceiling(nc / cpu)
mat3 = FBM(nr, nr, init=0)
reslist = bplapply(seq(1, nc, step), function(x){
  for(i in x:min(x+step-1, nc)){
    gid = which(assay(tcga.tumor[,i]) <= 2)
    mat3[gid, gid] = 1
  }
  NULL
})
mat3 = mat3[]
diag(mat3) = 0
mat3 = 1 - mat3
mat4 = mat3 - mat
mat4[mat4<0] = 0
rownames(mat4) = rownames(gtex)
colnames(mat4) = rownames(gtex)
g = graph.adjacency(mat4, mode='undirected')
elist = get.edgelist(g)
elist = elist[elist[,1] != elist[,2],]
write.table(elist, file='tx-spec.tsv', sep='\t', quote=F, row.names = F, col.names = F)

