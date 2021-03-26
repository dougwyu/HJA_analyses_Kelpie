

# ms <- difftime(strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
#          strptime("25/02/2021 20:56:09", format = "%d/%m/%Y %H:%M:%S", tz = "GMT"),
#          units = "mins")
# 
# 10000/ (6873 / as.numeric(ms))
# 
# 
# 
# 
# ms <- difftime(strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
#                strptime("27/02/2021 17:05:20", format = "%d/%m/%Y %H:%M:%S", tz = "GMT"),
#                units = "mins")
# 
# 81/as.numeric(ms)
# 
# 2000/3/60

t1 <- "26/03/2021 13:00:31"
t2 <- "26/03/2021 13:40:00"
#n <- 50*5 # number of models, nSteps * k
n <- 103

ms <- difftime(strptime(t2, format = "%d/%m/%Y %H:%M:%S", tz = "GMT"),
               strptime(t1, format = "%d/%m/%Y %H:%M:%S", tz = "GMT"),
               units = "mins")

ms
(tpm <- as.numeric(ms)/n)

k <- 5
nSteps <- 1000

# Total run time
(k*nSteps * tpm)/60
