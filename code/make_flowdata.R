# Jane Rogosch
# 1 June 2017
# Create flow metrics for known important fish-flow relationships
library(plyr)
library(dplyr)

 ?ddply
# ?match
# ?duplicated
# ?rle
# ?merge
# ?write.csv

# Make flowdata file with spring and summer floods and minimum baseflow duration

gagedata <- read.csv(file ="data/Paulden_USGS_gage1963-2017.csv", header = T)
head(gagedata)

maxcfs <- ddply(gagedata, .(year), summarize, cfs = max(cfs))

# Spring max flow
SPgagedata <- gagedata[gagedata$month >= 1 & gagedata$month <=4,]
SPmaxcfs <- ddply(SPgagedata, "year", summarize, cfs = max(cfs))
head(SPmaxcfs)
head(SPgagedata)

SPmaxcfs_merge<-merge(SPgagedata, SPmaxcfs, by = c("year","cfs") )
head(SPmaxcfs_merge)

SPmaxcfs_merge_unique <- SPmaxcfs_merge[!duplicated(SPmaxcfs_merge[,c('cfs','year')], fromLast=F),]

# Summer baseflow duration
# 25% flow is 22 cfs. So days less than or equal to that will be counted as baseflow.

SUgagedata <- gagedata[gagedata$month >=5 & gagedata$month <=9,]
length(SUgagedata$cfs[SUgagedata$cfs <= 22]) 

SUduration <- ddply(SUgagedata, "year", summarize, days=sum(cfs <= 22))
rle(SUgagedata$cfs<22)

SUduration2 <- ddply(SUgagedata, "year", summarize, length_rep = rle(cfs < 22)$lengths, rep_tf = rle(cfs < 22)$values)
SUduration2.5 <- SUduration2[SUduration2$rep_tf == TRUE,]
SUduration3 <- ddply(SUduration2.5, "year", summarize, max_baseflow_dur = max(length_rep)) 
SUduration3
SPmaxgage_SUdur_merge <- merge(SPmaxcfs_merge_unique, SUduration3,  by = "year", all.x = T)
head(SPmaxgage_SUdur_merge)

SPmaxgage_SUdur_merge[is.na(SPmaxgage_SUdur_merge)] <- 0
# ?boxplot.stats
quantile(SPmaxgage_SUdur_merge$max_baseflow_dur)

# Summer max flow
SUgagedata <- gagedata[gagedata$month >= 5 & gagedata$month <=9,]
SUmaxcfs <- ddply(SUgagedata, "year", summarize, cfs = max(cfs))

SUmaxcfs_merge<-merge(SUgagedata, SUmaxcfs, by =c("year","cfs") )

SUmaxcfs_merge_unique <- SUmaxcfs_merge[!duplicated(SUmaxcfs_merge[,c('cfs','year')], fromLast=F),]

# Altogether now
flowdata_Verde3 <- merge(SPmaxgage_SUdur_merge, SUmaxcfs_merge_unique, by = "year")
colnames(flowdata_Verde3)
flowdata_Verde2 <- flowdata_Verde3[,-c(6, 12, 13, 15)]
colnames(flowdata_Verde2) <- c("year", "sp_max_cfs", "agency", "site_no", "sp_datetime", "sp_calendar_day", "sp_water_day",
                               "sp_month", "basedur", "su_max_cfs", "su_datetime", "su_calendar_day", "su_water_day", "su_month")
flowdata_Verde2
flowdata_Verde <- flowdata_Verde2[,-c(3, 4, 5, 6, 8, 11, 12, 14)]
colnames(flowdata_Verde) <- c("Year", "SpFloodMag", "SpFloodDate", "BaseDur", "SuFloodMag", "SuFloodDate")
flowdata_Verde
### write.csv(flowdata_Verde, file = "data/flowdata_Verde.csv")                                  
