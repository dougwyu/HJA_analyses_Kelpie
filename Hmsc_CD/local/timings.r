

ms <- difftime(strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
         strptime("25/02/2021 20:56:09", format = "%d/%m/%Y %H:%M:%S", tz = "GMT"),
         units = "mins")

10000/ (6873 / as.numeric(ms))




ms <- difftime(strptime(Sys.time(), format = "%Y-%m-%d %H:%M:%S"),
               strptime("27/02/2021 17:05:20", format = "%d/%m/%Y %H:%M:%S", tz = "GMT"),
               units = "mins")

81/as.numeric(ms)

2000/3/60
