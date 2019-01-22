library(pROC)

args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1], header=T, sep = "\t",as.is=TRUE)

label <- df$label
FDR <- df$FDR

ROC1 <- roc(df$label,predictor= df$FDR, partial.auc = c(1, 0.9), partial.auc.correct = TRUE)
plot(ROC1, col = "grey")

AUC1 <- auc(ROC1)

print(AUC1)


TPR95<-coords(ROC1, 0.95, "specificity")
TPR95

SPE95<-coords(ROC1, 0.95, "sensitivity")
SPE95

out_df <- data.frame(AUC1,TPR95[3] , SPE95[2])

out_fname <- paste(args[1] ,".matrices" ,sep="")
write.table( out_df , out_fname , quote = F,row.names =F ,sep="\t")