library(qvalue)

args <- commandArgs(trailingOnly = TRUE)
df <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)
rep_num <-as.numeric(args[3])

M_mean<-df$M_mean
m_mean<-df$m_mean

std<-df$std
pvals<-c()

for(i in 1:length(M_mean)){	
	#print(paste( i,"," ,M_mean[i],",",m_mean[i], ",", std[i]))
	#print(pnorm( m_mean[i] , M_mean[i], sqrt(std[i]*std[i]/rep_num) ))
	pvals <- append(pvals, 2*pnorm( m_mean[i] , M_mean[i], sqrt(std[i]*std[i]/rep_num)  ))
}

FDR = lfdr(pvals, pi0 = 0.3)

#out_df <- data.frame(df$chr , df$pos , df$strand , df$label , df$M_mean, df$m_mean, df$std, pvals, FDR)
df$pvals<-pvals
df$FDR<-FDR
#out_fname <- paste(args[2] ,".out" ,sep="")
write.table( df , args[2] , quote = F,row.names =F ,sep="\t")
#print(out_fname);


