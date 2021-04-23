

## Create new data... 





## 3. Export data for prediction - get newData #####

## CHANGE RESOLUTION HERE/./// or above in resample


# allData <- values(allStack)
# dim(allData)
# indNA <- rowSums(is.na(allData)) == ncol(allData) # TRUE where all NAs
# 
# # remove NAs, convert to data frame and save index to replace after prediction
# allData <- data.frame(allData[!indNA, ])
# 
# ## change categorical to predictor values
# allData[, "insideHJA"] <- ifelse(is.na(allData[, "insideHJA"]), "no", "yes")
# table(allData[,"insideHJA"])
# 
# # cut.r == YrsSinceDisturb
# 
# summary(allData)
