library(trust)
library(aster)
args <- commandArgs(trailingOnly = TRUE)
df <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)
rep_num <- as.numeric(args[2])
lib_size <- as.numeric(args[3])
Mcv2<-df$Mcv2
mcv2<-df$mcv2
M_mtx <- matrix(0,nrow = length(df$M_mean),ncol=rep_num )
m_mtx <- matrix(0,nrow = length(df$M_mean),ncol=rep_num )
M_sum <- c()
M_mean <- c()
m_sum <- c()
m_mean <- c()

for(i in 1:length(df$M_mean)){

  M <- rktnb(rep_num, mu = df$M_mean[i], size =1/Mcv2[i] , k =0  )
  m <- rktnb(rep_num, mu = df$m_mean[i], size =1/mcv2[i] , k =0  )

  if(mean(M)< mean(m)){
	tmpvec <- M 
        M<-m
        m<-tmpvec
  }

 for(j in 1:rep_num){
  M_mtx[i,j] <-M[j] 
 }

 M_mean <- append(M_mean ,mean(M) )
 M_sum <- append(M_sum ,sum(M)  )

 for(j in 1:rep_num){
  m_mtx[i,j] <-m[j]
 }


 m_mean <- append(m_mean ,mean(m) )
 m_sum <- append(m_sum ,sum(m)  )


}

M_mean <- M_mean /(0.25*lib_size)
m_mean <- m_mean /(0.25*lib_size)

out_df <- data.frame(df$chr , df$pos , df$strand , df$label ,M_mean,m_mean ,M_sum,m_sum  )
colnames(out_df)<-c("chr" , "pos" , "strand" , "label" ,"M_mean", "m_mean" , "M_sum" ,"m_sum" )
out_fname <- paste(args[1] ,".out" ,sep="")
write.table( out_df , out_fname , quote = F,row.names =F ,sep="\t")

out_fname <- paste(args[1] ,".M_rc" ,sep="")
write.table( M_mtx , out_fname , quote = F,row.names =F , col.names =F ,sep="\t")
out_fname <- paste(args[1] ,".m_rc" ,sep="")
write.table( m_mtx , out_fname , quote = F,row.names =F , col.names =F ,sep="\t")
