library(qvalue)

args <- commandArgs(trailingOnly = TRUE)
df <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)
M<-df$M_sum
m<-df$m_sum

pvals<-c()

for(i in 1:length(M)){
	

	tdf <-  matrix(c(M[i], 0.5*(M[i]+m[i]), m[i], 0.5*(M[i]+m[i])), nrow = 2)       
        
	res <- fisher.test(tdf)
	if(res$p.value > 1) {
		pvals <-  append(pvals,  1)
	}else {
        	pvals <-  append(pvals,  res$p.value )
	}
}

#print(pvals)
FDR = lfdr(pvals, pi0 = 0.3)
df$pvals<-pvals
df$FDR<-FDR

write.table( df , args[2] , quote = F,row.names =F ,sep="\t")
