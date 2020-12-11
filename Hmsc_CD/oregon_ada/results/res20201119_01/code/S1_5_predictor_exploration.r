### Exploring predictors

# setwd("J:/UEA/Oregon")


## predictor choices - see exel doc for description
# cbind(X_cols_use)
# [1,] "elevation"         
# [2,] "canopy.ht"         
# [3,] "min.T"             
# [4,] "max.T"             
# [5,] "precipitation"     
# [6,] "metre.road"        
# [7,] "metre.stream"      
# [8,] "yrs.disturb.min"   
# [9,] "hja"    # hja inside (no logging really, nearer to primary forest) or outside experimental forest (logging)
# make a domindant land  cover variable around each point. % of structure... 

# [10,] "lysis.ratio"   #  ignore    lysis buffer ratio. different samples are different sizes - control for this?
# [11,] "spike.mean"        # 
# [12,] "l_Cover_2m_4m"     
# [13,] "l_Cover_2m_4m_all" 
# [14,] "l_Cover_2m_max"    
# [15,] "l_Cover_2m_max_all"
# [16,] "l_Cover_4m_16m"    # understorey
# [17,] "l_p25"             
# [18,] "l_p25_all"         
# [19,] "l_p95"             
# [20,] "l_p95_all"         
# [21,] "l_rumple"          
# [22,] "mean.NDVI"         
# [23,] "mean.EVI"          
# [24,] "mean.bright"       
# [25,] "mean.green"        
# [26,] "mean.wet" 