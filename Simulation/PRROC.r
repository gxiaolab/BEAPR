library(PRROC)
library(ggplot2)
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
        #jpeg(file = paste(filename, '.jpg', sep=""), width = width, height = height)
        #print(gplot)
        #dev.off()
}


#d_posdata <- as.matrix(read.table(args[1], header=F, sep = "\t",row.names = 1,as.is=TRUE))
#d_negdata <- as.matrix(read.table(args[2], header=F, sep = "\t",row.names = 1,as.is=TRUE))

#d_pr<-pr.curve(scores.class0 =d_posdata[,1], scores.class1 = d_negdata[,1], curve=TRUE,max.compute = TRUE)
#d_pr

#plot(d_pr,col=2)


#asbfname <- paste(args[1], "/asbclip.txt", sep="")

#df_asb <- read.table(asbfname, header=T, sep = "\t",as.is=TRUE)
df_asb <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)
#df_asb$FDR <-  p.adjust(df_asb$pvals, "BH")

neg_asb <- df_asb[df_asb$label ==0 ,]
pos_asb <- df_asb[df_asb$label ==1 ,]

#print(FDR)
#length(pos_asb$FDR)

s_pr<-pr.curve(scores.class0 =1-pos_asb$pvals, scores.class1 = 1-neg_asb$pvals, curve=TRUE)
s_pr
print(s_pr$auc.integral)
plot(s_pr)

#s_data <-data.frame(s_pr$curve)
#s_data["Methods"] <- rep("ACLIP",nrow(s_pr$curve) ) 
#print(s_data)

#chifname <- paste(args[1], "/chi.txt", sep="")
#df_chi <- read.table(chifname, header=T, sep = "\t",as.is=TRUE)
df_chi <- read.table(args[2], header=T, sep = "\t",as.is=TRUE)
#df_chi$FDR <-  p.adjust(df_chi$pvals, "BH")

neg_chi <- df_chi[df_chi$label ==0 ,]
pos_chi <- df_chi[df_chi$label ==1 ,]
#print(neg_chi$FDR)

t_pr<-pr.curve(scores.class0 = 1-pos_chi$pvals, scores.class1 =  1-neg_chi$pvals, curve=TRUE)
t_pr
print(t_pr$auc.integral)
#plot(t_pr)

#t_data <-data.frame(t_pr$curve)
#t_data["Methods"] <- rep("Chi-squared",nrow(t_pr$curve) )

df_fish <- read.table(args[3], header=T, sep = "\t",as.is=TRUE)
neg_fish <- df_fish[df_fish$label ==0 ,]
pos_fish <- df_fish[df_fish$label ==1 ,]
f_pr<-pr.curve(scores.class0 = 1-pos_fish$pvals, scores.class1 =  1-neg_fish$pvals, curve=TRUE)
f_pr
print(f_pr$auc.integral)

df <- rbind(s_pr$auc.integral, t_pr$auc.integral, f_pr$auc.integral)

out_fname <- paste(args[1] ,".pr" ,sep="")
write.table( df , out_fname , quote = F,row.names =F ,sep="\t")


#print(PR_data)

#df$Methods <- factor(df$Methods)
#p <- ggplot(df, aes(group=Methods, x=X1,y=X2 ,colour=Methods) )

#p <- p + theme(axis.text=element_text(size=18),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1) , axis.title=element_text(size=18),plot.title=element_text(size=14, face="bold"),legend.text=element_text(size=14), legend.title=element_text(size=14), #,face="bold"),
#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#panel.background = element_blank(),
#axis.line.x = element_line(color="black", size = .5),
#axis.line.y = element_line(color="black", size = .5)
#)

#p <- p  + #xlim( "[10,20)","[20,30)","[30,40)","[40,50)","[50,60)","[60,70)","[70,80)","[80,90)","[90,100)","[100,-)" )+
#    geom_line() +
    #scale_x_discrete()+
#    ylab("Precision") +xlab("Recall")
#ExportPlot(p,"PRROC")



