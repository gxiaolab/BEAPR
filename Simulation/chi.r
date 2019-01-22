library(qvalue)

args <- commandArgs(trailingOnly = TRUE)
df <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)
M<-df$M_sum
m<-df$m_sum

pvals<-c()

for(i in 1:length(M)){
	       
        res<-chisq.test(c(M[i],m[i]))

        pvals <-  append(pvals,  res$p.value )
}

print(pvals)
#FDR = p.adjust(pvals, "BH")
FDR = lfdr(pvals, pi0 = 0.3)
df$pvals<-pvals
df$FDR<-FDR

write.table( df , args[2] , quote = F,row.names =F ,sep="\t")
