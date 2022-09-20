library(data.table)
library(survival)
library(parallel)
library(limma)
args = commandArgs()
f = args[4]
of = paste0(f, '-cox')
sl = fread(f, head=T)
qnorm.array <- function(mat)
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}
load('isle-gene-sample-table.rda')
gt = alias2SymbolTable(prob$genes)
sl$g1 = match(sl$g1, gt)
sl$g2 = match(sl$g2, gt)
sl = sl[!is.na(sl$g1) & !is.na(sl$g2)]
sl = data.matrix(sl)
load("ISLE/prob.TCGA.RData")
numGenes = nrow(prob$scna)
numSamples = ncol(prob$scna)
mRNA.norm=prob$mRNA.norm
scna.norm=prob$scna.norm
surv.dt=prob$surv.dt
age=qnorm.array(prob$age)
sex=as.character(prob$sex)
race=prob$race
types=prob$types
mRNAq2=prob$mRNAq2
scnaq2=prob$scnaq2
GII=qnorm.array(prob$GII)
cox.pair.sl = function(pair,use.mRNA=F)
{
        if(use.mRNA){
                g1 = mRNA.norm[pair[1],]
                g2 = mRNA.norm[pair[2],]
                f1 = mRNAq2[pair[1],]
                f2 = mRNAq2[pair[2],]
        }else{  
                g1 = scna.norm[pair[1],]
                g2 = scna.norm[pair[2],]
                f1 = scnaq2[pair[1],]
                f2 = scnaq2[pair[2],]
        }
        res=rep(c(0,1),2)
        if (sum(!is.na(g1))>100 & sum(!is.na(g2))>100 & sum(!is.na(f1))>100 & sum(!is.na(f2))>100){
            cov = ifelse(f1 ==0 & f2==0,1,0 )
            dt1 = data.frame(cbind(surv.dt, cov, age, sex, GII, race, types, cbind(g1 , g2)))
            cox.out = coxph(Surv(time,status) ~ cov + strata(types), data=dt1)
            aa  = summary(cox.out)        
            cox.out = coxph(Surv(time,status) ~ cov + g1 + g2 +strata(types), data=dt1)
            bb  = summary(cox.out)        
            cox.out = coxph(Surv(time,status) ~ cov + GII+age+ strata(types,sex,race), data=dt1)
            cc  = summary(cox.out)        
            cox.out = coxph(Surv(time,status) ~ cov + g1 + g2 +GII+age+strata(types,sex,race), data=dt1)
            dd  = summary(cox.out)        
            res=c(aa$coefficients["cov",c(1,5)], bb$coefficients["cov",c(1,5)],cc$coefficients["cov",c(1,5)], dd$coefficients["cov",c(1,5)])
	}
	return(res)
}
cox.mRNA = mclapply(1:nrow(sl), function(tt) {print(tt); cox.pair.sl(sl[tt,], use.mRNA=T)}, mc.cores=24)
cox.sl.mRNA = t(do.call(cbind, cox.mRNA))
cox.scna = mclapply(1:nrow(sl), function(tt) {print(tt); cox.pair.sl(sl[tt,], use.mRNA=F)}, mc.cores=24)
cox.sl.scna = t(do.call(cbind, cox.scna))
sl = data.table(sl, cox.sl.mRNA[,c(1,2,7,8)], cox.sl.scna[,c(1,2,7,8)])
colnames(sl) = c('g1', 'g2', 'm1', 'm2', 'm7', 'm8', 'c1', 'c2', 'c7', 'c8')
sl$g1 = prob$genes[sl$g1]
sl$g2 = prob$genes[sl$g2]
write.table(sl, file=of, quote=F, row.names=F, sep='\t')
