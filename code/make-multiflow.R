## This is code to make a list of all possible future flow regimes.
## This will then get read into the main model as a list. 
#rm(list = ls()) # clearing the workspace 
library(plyr)
library(tidyverse)

## Loading functions from functions.R file--------------------------------------
source('code/functions.R')

## Make flowdata file with spring and summer floods and minimum baseflow dur
## Read data in
gagedata <- read.csv(file ="data/Paulden_USGS_gage1963-2017.csv", header = T)

## Remove 1963
gagedata <- gagedata %>%
    filter(year != 1963)

## max flow of year ------------------------------------------------------------
maxcfs <- gagedata %>%
    group_by(year) %>%
    summarise(cfs = max(cfs))

## Spring max flow -------------------------------------------------------------
## Filtering to just get spring months
SPgagedata <- gagedata %>%
    filter(month >= 1 & month <=4)

## Calculating max flow in these months each year
SPmaxcfs <- SPgagedata %>%
    group_by(year) %>%
    summarise(cfs = max(cfs))
    
## merging with gage data (note US spelling leftover from Jane) ;) 
SPmaxcfs_merge <- inner_join(SPgagedata, SPmaxcfs, by = c('year', 'cfs'))

## Removing duplicates of year and cfs
## Some cfs are repeated across multiple dates but it's only for the low flows
##so not a problem here
## i.e. big max flows only have one peak
SPmaxcfs_merge_unique <- SPmaxcfs_merge %>%
    group_by(year, cfs) %>%
    slice(1) %>%
    ungroup

# Summer baseflow duration -----------------------------------------------------
# 25% flow is 22 cfs. So days less than that will be counted as baseflow.

## Filtering to just get spring months
SUgagedata <- gagedata %>%
    filter(month >= 5 & month <=9)

## Calculating baseflow duration
## Less than 22 cfs
## Use rle to calculate sequences. Select those that meet the requirements.
## This leads to multiple values for some years and none for others. So, I do
## this inside 'do'. For some reason 'do' can't keep zero-length groups so
## I have to right join the results of the summarize where I take the max
## length of each year with a vector of all years then fill the NAs with 0.
SUduration <- SUgagedata %>%
    group_by(year) %>%
    do({tmp <- rle(.$cfs < 22)
        data.frame(lengths_T = tmp$lengths[tmp$values == TRUE])
    } 
    ) %>%
    summarise(max_baseflow_dur = max(lengths_T)) %>%
    right_join(data.frame(year = unique(SUgagedata$year))) %>%
    mutate(max_baseflow_dur = ifelse(is.na(max_baseflow_dur),
                                     0,
                                     max_baseflow_dur))

## Joining these results with the earlier and renaming
SPmax_SUdur_merge <- left_join(SPmaxcfs_merge_unique, SUduration) %>%
    select(year, cfs, water_day, max_baseflow_dur) %>%
    rename(SpFloodMag = cfs,
           SpFloodDate = water_day,
           BaseDur = max_baseflow_dur)

# Summer max flow --------------------------------------------------------------
SUmaxcfs <- SUgagedata %>%
    group_by(year) %>%
    summarise(cfs = max(cfs))

## Joining with gage data
SUmaxcfs_merge <- inner_join(SUgagedata, SUmaxcfs, by = c('year', 'cfs'))

## Taking unique year-cfs combos. Note slice again. 
SUmaxcfs_merge_unique <- SUmaxcfs_merge %>%
    group_by(year, cfs) %>%
    slice(1) %>%
    ungroup %>%
    select(year, cfs, water_day) %>%
    rename(SuFloodMag = cfs,
           SuFloodDate = water_day)

## Joining it all together -----------------------------------------------------
flowdata_Verde <- left_join(SPmax_SUdur_merge, SUmaxcfs_merge_unique)

## Exporting
write.csv(flowdata_Verde, file = "data/flowdata_Verde.csv")    

## -----------------------------------------------------------------------------
## Making flow yeartypes for simulations ---------------------------------------
## -----------------------------------------------------------------------------

## Bringing in flow data
## Verde flow data at Paulden 7/17/1963-2017, 54 years continuous
raw.flow <- read.csv("data/flowdata_Verde.csv") 

## FloodMag - magnitude of flood in cfs
## BaseDur - summer baseflow duration in days
## FloodDate is peak dates of all floods (Oct 1 = 1) because using water year

## Visualising
raw.flow %>%
    gather(metric,
           val,
           SpFloodMag:SuFloodMag,
           -Year,
           -X,
           -SpFloodDate,
           -SuFloodDate) %>%
    ggplot(aes(Year, val)) +
    facet_wrap(~metric) +
    geom_line() +
    geom_hline(yintercept = c(700, 200, 220), colour = 'red')

## Flow thresholds -------------------------------------------------------------

## Magnitude of peak flow over which is considered a large flood event
SP_highfloodcutoff = 700
## SP: Jan 1 - Apr 30 (water_day 93-213).
## Floods (cfs) at 4-yr recurrence interval or greater (30 times median flow).
SU_highfloodcutoff = 200
## SU: May 1 - Sep 30 (Water_day 214-365).
## Floods (Cfs) at 4-year recurrence interval for summer
medfloodcutoff = 220
## SP Floods at bankfull flood, 2.5-yr recurrence interval (10 X median flow)
## Non-event is same for drought (threshold below rather than above)
mindroughtlength = 40
## Threshold length of low-flow in days
## (greater than 75th percentile of low flow duration) -- determined as number
## of consecutive days with discharge less than 22 cfs
## (25th percentile of flows) May 1 - Sep 30
## This threshold is compared to baseflow duration in any given year --
## a vector specified in csv (baseflow_dur).


## Defining year types here ----------------------------------------------------
SP_highflood <- SP_highflood_func(raw.flow$SpFloodMag) 

SU_highflood <- SU_highflood_func(raw.flow$SuFloodMag)

medflood <- medflood_func(raw.flow$SpFloodMag)

drought <- drought_func(Spfl = raw.flow$SpFloodMag,
                        BD = raw.flow$BaseDur,
                        Sufl = raw.flow$SuFloodMag)

nonevent <- nonevent_func(Spfl = raw.flow$SpFloodMag,
                          BD = raw.flow$BaseDur,
                          Sufl = raw.flow$SuFloodMag) 


## year type (yt) flow
natural.flow <- as.data.frame(cbind(year = raw.flow$Year,
                                    SP_highflood,
                                    SU_highflood,
                                    medflood,
                                    drought,
                                    nonevent)) 
natural.flow$rep <- as.numeric(as.character(rownames(natural.flow)))

## Reordering
natural.flow <- natural.flow %>%
    select(year, rep, SP_highflood:nonevent)

apply(natural.flow[,-c(1:2)], 1, sum)
apply(natural.flow[,-c(1:2)], 2, sum)

natural.flow %>% mutate(sum = rowSums(.[,3:7]))

## medflood is never found with summer highflood, so I won't allow that combo.
## We thus have nonevent, drought, medflood, and any combo of the 2 highfloods.
## i.e. each individually, and combined

## Making flow scenarios -------------------------------------------------------
## Sequentially converting each year-type to 0 or 1

head(natural.flow)
str(natural.flow)

## Increased drought
drought.list <- list()

for(i in 1:length(natural.flow[,1])) {
    
    nm <- paste0('drought_', i)
    
    dat <- natural.flow %>% 
        mutate(SP_highflood = ifelse(
                   rep <= i, 0, SP_highflood),
               SU_highflood = ifelse(
                   rep <= i, 0, SU_highflood),
               medflood = ifelse(
                   rep <= i, 0, medflood),
               drought = ifelse(
                   rep <= i, 1, drought),
               nonevent = ifelse(
                   rep <= i, 0, nonevent)
               )
    
    drought.list[[nm]] <- dat
    
}

## Checking to see it's working.
## Should finish w/ 54 drought years and 0 others
plyr::ldply(drought.list, function(x) apply(x, 2, sum))

## Removing duplicated dataframes. They can be as adding a year type can be the
## same as it already was. 
drought.list.u <- drought.list[!duplicated(drought.list)]
plyr::ldply(drought.list.u, function(x) apply(x, 2, sum))

## Increased nonevent
nonevent.list <- list()

for(i in 1:length(natural.flow[,1])) {
    
    nm <- paste0('nonevent_', i)
    
    dat <- natural.flow %>% 
        mutate(SP_highflood = ifelse(
                   rep <= i, 0, SP_highflood),
               SU_highflood = ifelse(
                   rep <= i, 0, SU_highflood),
               medflood = ifelse(
                   rep <= i, 0, medflood),
               drought = ifelse(
                   rep <= i, 0, drought),
               nonevent = ifelse(
                   rep <= i, 1, nonevent)
               )
    
    nonevent.list[[nm]] <- dat
    
}

## Checking to see it's working.
## Should finish w/ 54 nonevent years and 0 others
plyr::ldply(nonevent.list, function(x) apply(x, 2, sum))

## Removing duplicated dataframes. They can be as adding a year type can be the
## same as it already was. 
nonevent.list.u <- nonevent.list[!duplicated(nonevent.list)]
plyr::ldply(nonevent.list.u, function(x) apply(x, 2, sum))


## Increased medflood
medflood.list <- list()

for(i in 1:length(natural.flow[,1])) {
    
    nm <- paste0('medflood_', i)
    
    dat <- natural.flow %>% 
        mutate(SP_highflood = ifelse(
                   rep <= i, 0, SP_highflood),
               SU_highflood = ifelse(
                   rep <= i, 0, SU_highflood),
               medflood = ifelse(
                   rep <= i, 1, medflood),
               drought = ifelse(
                   rep <= i, 0, drought),
               nonevent = ifelse(
                   rep <= i, 0, nonevent)
               )
    
    medflood.list[[nm]] <- dat
    
}

## Checking to see it's working.
## Should finish w/ 54 medflood years and 0 others
plyr::ldply(medflood.list, function(x) apply(x, 2, sum))

## Removing duplicated dataframes. They can be as adding a year type can be the
## same as it already was. 
medflood.list.u <- medflood.list[!duplicated(medflood.list)]
plyr::ldply(medflood.list.u, function(x) apply(x, 2, sum))

## Increased SP_highflood
SP_highflood.list <- list()

for(i in 1:length(natural.flow[,1])) {
    
    nm <- paste0('SP_highflood_', i)
    
    dat <- natural.flow %>% 
        mutate(SP_highflood = ifelse(
                   rep <= i, 1, SP_highflood),
               SU_highflood = ifelse(
                   rep <= i, 0, SU_highflood),
               medflood = ifelse(
                   rep <= i, 0, medflood),
               drought = ifelse(
                   rep <= i, 0, drought),
               nonevent = ifelse(
                   rep <= i, 0, nonevent)
               )
    
    SP_highflood.list[[nm]] <- dat
    
}

## Checking to see it's working.
## Should finish w/ 54 SP_highflood years and 0 others
plyr::ldply(SP_highflood.list, function(x) apply(x, 2, sum))

## Removing duplicated dataframes. They can be as adding a year type can be the
## same as it already was. 
SP_highflood.list.u <- SP_highflood.list[!duplicated(SP_highflood.list)]
plyr::ldply(SP_highflood.list.u, function(x) apply(x, 2, sum))

## Increased SU_highflood
SU_highflood.list <- list()

for(i in 1:length(natural.flow[,1])) {
    
    nm <- paste0('SU_highflood_', i)
    
    dat <- natural.flow %>% 
        mutate(SP_highflood = ifelse(
                   rep <= i, 0, SP_highflood),
               SU_highflood = ifelse(
                   rep <= i, 1, SU_highflood),
               medflood = ifelse(
                   rep <= i, 0, medflood),
               drought = ifelse(
                   rep <= i, 0, drought),
               nonevent = ifelse(
                   rep <= i, 0, nonevent)
               )
    
    SU_highflood.list[[nm]] <- dat
    
}

## Checking to see it's working.
## Should finish w/ 54 SU_highflood years and 0 others
plyr::ldply(SU_highflood.list, function(x) apply(x, 2, sum))
 
## Removing duplicated dataframes. They can be as adding a year type can be the
## same as it already was. 
SU_highflood.list.u <- SU_highflood.list[!duplicated(SU_highflood.list)]
plyr::ldply(SU_highflood.list.u, function(x) apply(x, 2, sum))

## Increased SPSU_highflood
SPSU_highflood.list <- list()

for(i in 1:length(natural.flow[,1])) {
    
    nm <- paste0('SPSU_highflood_', i)
    
    dat <- natural.flow %>% 
        mutate(SP_highflood = ifelse(
                   rep <= i, 1, SP_highflood),
               SU_highflood = ifelse(
                   rep <= i, 1, SU_highflood),
               medflood = ifelse(
                   rep <= i, 0, medflood),
               drought = ifelse(
                   rep <= i, 0, drought),
               nonevent = ifelse(
                   rep <= i, 0, nonevent)
               )
    
    SPSU_highflood.list[[nm]] <- dat
    
}

## Checking to see it's working.
## Should finish w/ 54 SP and SU_highflood years and 0 others
plyr::ldply(SPSU_highflood.list, function(x) apply(x, 2, sum))

## Removing duplicated dataframes. They can be as adding a year type can be the
## same as it already was. 
SPSU_highflood.list.u <- SPSU_highflood.list[!duplicated(SPSU_highflood.list)]
plyr::ldply(SPSU_highflood.list.u, function(x) apply(x, 2, sum))

## Putting all in a single list
all.scenarios.list <- c(drought.list.u,
                        nonevent.list.u,
                        medflood.list.u,
                        SP_highflood.list.u,
                        SU_highflood.list.u,
                        SPSU_highflood.list.u
                        )

map(all.scenarios.list, head)
map(all.scenarios.list, class)
names(all.scenarios.list)

all.scenarios.list$natural.flow <- natural.flow

## Save natural flow through to 2017
write.csv(natural.flow, 'data/naturalflow.csv', row.names = FALSE)

## Save the objects as .rds files - then use loadRDS in other file. 
saveRDS(all.scenarios.list, 'data/all_scenarios_list.rds')

