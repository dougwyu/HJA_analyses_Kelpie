#### Comparsion gpu cpu analysis


wd <- here::here()
setwd(wd)
getwd()
rm(wd)

# save(model_cpu_lm, model_gpu_lm, coef_cpu_lm, coef_gpu_lm,
#      model_cpu_DNN, model_gpu_DNN, coef_cpu_DNN, coef_gpu_DNN,
#      file ="results_sjSDM/cpu_gpu_cf.rdata")
load("Hmsc_CD/oregon_ada/results_sjSDM/cpu_gpu_cf.rdata")
load("Hmsc_CD/oregon_ada/results_sjSDM/cpu_gpu_cf2.rdata")
str(model_cpu_lm,max.level = 1)

# compare linear models
str(coef_cpu_lm, max.level = 2)

coef_cpu_lm$env[[1]][1:10, 1:10]
coef_gpu_lm$env[[1]][1:10, 1:10]

# coefficient names
preds <- c(model_cpu_lm$names, "spatial")

range(coef_gpu_lm$env[[1]])
range(coef_cpu_lm$env[[1]])

range(coef_gpu_lm$spatial[[1]])
range(coef_cpu_lm$spatial[[1]])

# number of coefficients per species
n <- ncol(coef_gpu_lm$env[[1]]) + 1

# make single matrix of env and spatial coefficients
cpu <- cbind(coef_cpu_lm$env[[1]], coef_cpu_lm$spatial[[1]])
gpu <- cbind(coef_gpu_lm$env[[1]], coef_gpu_lm$spatial[[1]])

cpu1 <- cbind(coef_cpu_lm1$env[[1]], coef_cpu_lm1$spatial[[1]])
cpu2 <- cbind(coef_cpu_lm2$env[[1]], coef_cpu_lm2$spatial[[1]])

# plot(1, 1, xlim = c(0, n+1), ylim = range(c(coef_cpu_lm$env[[1]][,1],coef_gpu_lm$env[[1]][,1])), type= "n")
# segments(x0 = 1:n-0.5, y0 = coef_cpu_lm$env[[1]][,1], x1 = 1:n+0.5[1], y1 = coef_gpu_lm$env[[1]][,1])
# segments(x0 = 0.5, y0 = coef_cpu_lm$env[[1]][,1], x1 = 1.5, y1 = coef_gpu_lm$env[[1]][,1])
# segments(x0 = 2, y0 = coef_cpu_lm$env[[1]][,1], x1 = 3.5, y1 = coef_cpu_lm$env[[1]][,1])

summary(cpu)

diff <- (cpu-gpu)/cpu
range(diff)
summary(diff)
brq <- quantile(diff, probs = seq(0,1,0.01))
brq

cols <- colorRampPalette(c("red", "white", "blue"))
#breaks <- c(-2, -1, -0.2, -0.1, 0, 0.1, 0.2, 1, 2)
breaks <- c(-2, -1, 0, 1, 2)
breaks <- quantile(c(cpu, gpu))

breaks <- quantile(c(cpu, cpu1, cpu2))

nCol <- length(breaks)-1

symbols(1:nCol, rep(1, nCol), squares = rep(4, nCol), bg = cols(nCol))

png("Hmsc_CD/local/plots/gpu_cpu_compare_brksQ.png", width = 250, height = 150, units = "mm", res = 100)
op <- par(mfrow = c(1,3), mar = c(12, 5, 2, 1))
image(t(cpu), col = cols(nCol), breaks = breaks, ylab = "species", yaxt = "n", xaxt = "n", main = "cpu")
axis(side = 1, at = seq(0,1, length.out = n), labels = preds, las = 2)
image(t(gpu), col = cols(nCol), breaks = breaks, xlab = "predictors", yaxt = "n", xaxt = "n", main = "gpu")
image(t(diff), col = cols(length(brq)-1), breaks = brq, xlab = "% difference", yaxt = "n", xaxt = "n")
dev.off()

par(op)

# compare three identical runs

cpu.arr <- array(c(cpu, cpu1, cpu2), dim = c(dim(cpu), 3))

cpu[1:5,1:5]
cpu.arr[1:5,1:5,1]

cpu1[1:5,1:5]
cpu.arr[1:5,1:5,2]

mn.cpu <- apply(cpu.arr, c(1,2), mean)
mean(c(cpu[1,1], cpu1[1,1], cpu2[1,1]))
mn.cpu[1,1]

sd.cpu <- apply(cpu.arr, c(1,2), sd)
cv.cpu <- sd.cpu/mn.cpu
summary(cv.cpu)
range(cv.cpu)
quantile(cv.cpu)

breaks.cv <- c(range(cv.cpu)[1], -0.25, 0.25, range(cv.cpu)[2])

png("Hmsc_CD/local/plots/gpu_cpu_compare_3runs_lm_cpu.png", width = 250, height = 150, units = "mm", res = 100)
op <- par(mfrow = c(1,4), mar = c(12, 5, 2, 1))
image(t(cpu), col = cols(nCol), breaks = breaks, ylab = "species", yaxt = "n", xaxt = "n", main = "cpu0")
axis(side = 1, at = seq(0,1, length.out = n), labels = preds, las = 2)
image(t(cpu1), col = cols(nCol), breaks = breaks, xlab = "predictors", yaxt = "n", xaxt = "n", main = "cpu1")
image(t(cpu2), col = cols(nCol), breaks = breaks, yaxt = "n", xaxt = "n", main = "cpu2")
image(t(cv.cpu), col = cols(length(breaks.cv)-1), breaks = breaks.cv, yaxt = "n", xaxt = "n", main = "cpu cv")
dev.off()



## Look at DNN coefficients
str(coef_cpu_DNN, max.level = 1)
str(coef_cpu_DNN$env, max.level = 1)

coef_cpu_DNN$env[[8]]
coef_gpu_DNN$env[[8]]

range(coef_cpu_DNN$env[[8]])
range(coef_gpu_DNN$env[[8]])
breaks <- quantile(c(coef_cpu_DNN$env[[8]], coef_gpu_DNN$env[[8]]))
breaks

range(coef_cpu_DNN$env[[7]])
range(coef_gpu_DNN$env[[7]])

summary(coef_cpu_DNN$env[[7]])

breaks <- quantile(c(coef_cpu_DNN$env[[7]], coef_gpu_DNN$env[[7]]))

nCol = length(breaks)-1

png("Hmsc_CD/local/plots/gpu_cpu_DNN_compare_brksQ.png", width = 250, height = 150, units = "mm", res = 100)
op <- par(mfrow = c(1,2), mar = c(12, 5, 2, 1))
image(t(coef_cpu_DNN$env[[7]]), col = cols(nCol), breaks = breaks, ylab = "species", yaxt = "n", xaxt = "n", main = "cpu")
image(t(coef_gpu_DNN$env[[7]]), col = cols(nCol), breaks = breaks, xlab = "DNN factors", yaxt = "n", xaxt = "n", main = "gpu")
dev.off()

png("Hmsc_CD/local/plots/gpu_cpu_DNN_compare_spp_Coef.png", width = 250, height = 150, units = "mm", res = 100)
op <- par(mfrow = c(1,2), mar = c(2, 4, 2, 1))
image(t(coef_cpu_DNN$env[[8]]), col = cols(nCol), breaks = breaks, ylab = "species", yaxt = "n", xaxt = "n", main = "cpu")
image(t(coef_gpu_DNN$env[[8]]), col = cols(nCol), breaks = breaks, yaxt = "n", xaxt = "n", main = "gpu")
dev.off()

# save(results_new_vars_cv, results_old_vars_cv, file ="results_sjSDM/oregon_trial_cf_vars_cv.rdata")
load("Hmsc_CD/oregon_ada/results_sjSDM/oregon_trial_cf_vars_cv.rdata")

str(results_new_vars_cv, max.level = 1)
str(results_new_vars_cv$tune_results, max.level = 1)
str(results_new_vars_cv$tune_results[[1]][[1]], max.level = 1)

head(results_new_vars_cv$summary, 30)

head(results_new_vars_cv$short_summary[order(results_new_vars_cv$short_summary$AUC_test, decreasing = T),])

head(results_old_vars_cv$short_summary[order(results_old_vars_cv$short_summary$AUC_test, decreasing = T),])




