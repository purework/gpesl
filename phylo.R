library(data.table)
library(magrittr)
args = commandArgs()
f = args[4]
newgene = 'old'
if(length(args) > 4){ newgene = args[5] }
rk = 'norank'
if(length(args) > 5) { rk = args[6] }
gm = fread('all-gene-alias-map.tsv', head=T)
setkey(gm, 'old')
# function to identify phylogenetic distance of a pair of genes
phylo.profile = function(sl.gene.all){
  sl.gene1 = sl.gene.all
  phylo$genes = gm[as.character(phylo$genes),]$new
  sl.phylo =  cbind(match(sl.gene1[,1], phylo$genes), match(sl.gene1[,2], phylo$genes))
  featureMat = data.matrix((phylo[sl.phylo[,1],-(1:3)] - phylo[sl.phylo[,2],-(1:3)])^2)
  featureMat %*% t(feature.weight)
}
# The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
load("ISLE/data/yuval.phylogenetic.profile.RData")
# the feature weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)
load("ISLE/data/feature.weight.RData")
pair = fread(f)
pair$pp=phylo.profile(data.frame(pair[,1:2], stringsAsFactors = F)) %>% as.numeric
# select candidate SL pairs that pass have small phylogenetic distance mesured by NMF
if(rk == 'norank'){
    res = pair[pp<10.5,]
}
if(rk == 'rank'){
    pair$pp.rank = rank(pair$pp, na.last="keep") / length(pair$pp)
    res = pair
}
write.table(res, file=paste0(f, '-phylo'), sep='\t', quote=F, row.names = F)
