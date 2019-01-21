library(ggplot2)
library(RColorBrewer)
args <- commandArgs(trailingOnly = TRUE)

ExportPlot <- function(gplot, filename, width=5, height=4.5) {
        # Export plot in PDF and EPS.
        # Notice that A4: width=11.69, height=8.27
        pdf(file = paste(filename, '.pdf', sep=""), width = width, height = height)
        print(gplot)
        dev.off()
        postscript(file = paste(filename, '.eps', sep=""), width = width, height = height)
        print(gplot)
        dev.off()
}

df <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)

if(nrow(df) > 10000 ){
	df <- df[sample(nrow(df),10000),]

}


print(max(df$xmean))
Max_xmean <- max(df$xmean)
print(min(df$xmean))
Min_xmean <- min(df$xmean)

myloess <- loess(df$cv2 ~ log2(df$xmean), data = df)

predY <- predict(myloess,log2(df$xmean), se = TRUE)

#print(predY$s)

# expVar <- df$xmean*df$xmean*(predY$fit)
# std <- sqrt(expVar)
# #Res <- data.frame(df$id,df$label,df$allele,df$xmean,df$cv2,predY$fit+2*predY$s,expVar,std)
# Res <- data.frame(df$id,df$xmean,df$cv2,predY$fit ,expVar,std)
#colnames(Res)<-c("ID","xmean","cv2","pred_cv2","expectedVar","ex_std")
#write.table(Res, args[2]  ,sep="\t",quote=F, col.names = T, row.names = F)

p1 <- ggplot(df, aes(x = log2(df$xmean), y = df$cv2 ))+geom_hex(bins = 15, alpha = 0.75) +geom_line(aes(y=predY$fit),colour="red")+
        scale_fill_gradient2(high  = "#2c7fb8",
                             mid  = "#2c7fb8",
                             low = "#ece7f2",  midpoint = 400) +
      labs(x = "log2(mean)",y="CV2")+
      theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=10),
        plot.title=element_text(size=8, face="bold"),
        strip.text.x = element_text(size=10), strip.background = element_rect(fill='white'),
        legend.text=element_text(size=10), legend.key=element_rect(fill="white"),
        legend.title=element_text(size=10), #,face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )

#outpname <- paste(args[2],".png",sep="");
#ggsave(outpname, dpi = 300)

 ExportPlot(p1,args[1])

 df_pred <- read.table(args[2], header=T, sep = "\t",as.is=TRUE)
 rep_num <- as.numeric(args[4])

 M_mean<-df_pred$M
 m_mean<-df_pred$m

 pvals<-c()
 M <-c()
 
 for(i in 1:length(M_mean)){
   M_val <- M_mean[i]
   if(M_mean[i]>Max_xmean){
      M_val <- Max_xmean
   }
   if(M_mean[i]<Min_xmean){
      M_val <- Min_xmean
   }

   M <- append(M, M_val ) 
 }

predM <- predict(myloess,log2(M), se = TRUE)
expVar <- M*M*(predM$fit)
std <- sqrt(expVar)
out_std <- c()

for(i in 1:length(M)){
   
   pvals <- append(pvals, 2*pnorm( m_mean[i] , M_mean[i], sqrt(std[i]*std[i]/rep_num) )) 
   out_std <- append(out_std, sqrt(std[i]*std[i]/rep_num) )
 }
FDR <- p.adjust(pvals, "BH")

#out_df <- data.frame( df_pred$id ,  df_pred$M, df_pred$m, out_std, pvals,FDR)
out_df <- df_pred
out_df$expVar<-expVar
out_df$exp_std<-out_std
out_df$pvals<-pvals
out_df$fdr<-FDR
#colnames(out_df) <- c("id" ,  "M_mean", "m_mean", "exp_std", "pvals","FDR")
write.table( out_df , args[3] , quote = F,row.names =F ,sep="\t")

