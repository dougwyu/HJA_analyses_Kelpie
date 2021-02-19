

getwd()

getwd()
wd <- here::here()
setwd(wd)
dir()

load("Hmsc_CD/oregon_ada/code_sjSDM/r20210217/results/sp_results.rdata")


png("Hmsc_CD/local/plots/eval_metrics_pairs_test.png")
plot(sp.res.test[[1]])
mtext("test")
dev.off()

png("Hmsc_CD/local/plots/eval_metrics_pairs_train.png")
plot(sp.res.train[[1]])
mtext("test")
dev.off()


png("Hmsc_CD/local/plots/eval_metrics_auc.png")
plot(sp.mn.train$auc, sp.mn.test$auc, xlim = c(0,1), ylim = c(0,1))
abline(0,1)
dev.off()




