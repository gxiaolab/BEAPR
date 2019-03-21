./genMean 10 0.8
Rscript Sim.r data/mean 2 10
# The following command lines are for reproducing results on simulated data without cross-linking bias 
#Rscript Reg.r reg.rc data/mean.out data/asbclip.in
#Rscript Pred.r data/asbclip.in data/asbclip.out 2
#Rscript chi.r data/mean.out data/chi.out
#Rscript fisher.r data/mean.out data/fisher.out
#Rscript pROC_new.r data/asbclip.out
#Rscript pROC_new.r data/chi.out
#Rscript pROC_new.r data/fisher.out
#Rscript PRROC.r data/asbclip.out  data/chi.out data/fisher.out