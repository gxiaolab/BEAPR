./genMean 10 0.8
Rscript Sim.r /home/yyang/data/eclip/K562/SRSF1_test/sim/mean 2 10
Rscript Reg.r /home/yyang/data/eclip/K562/SRSF1_test/tmp/reg.rc /home/yyang/data/eclip/K562/SRSF1_test/sim/mean.out /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.in
Rscript Pred.r /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.in /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out 2
Rscript chi.r /home/yyang/data/eclip/K562/SRSF1_test/sim/mean.out /home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out
Rscript pROC_new.r /home/yyang/data/eclip/K562/SRSF1_test/sim/asbclip.out
Rscript pROC_new.r /home/yyang/data/eclip/K562/SRSF1_test/sim/chi.out