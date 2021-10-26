## PA DATA
var <- apply(otu.pa.csv, 2, var)
mn <- colMeans(otu.pa.csv)

plot(mn, var, log = "xy")

## QP DATA
var <- apply(otu.qp.csv, 2, var)
mn <- colMeans(otu.qp.csv)

plot(mn, var, log = "xy")
