## This is code to make a list of all possible future flow regimes.
## This will then get read into the main model as a list. 
#rm(list = ls()) # clearing the workspace 

library(tidyverse)

## Loading functions from functions.R file--------------------------------------
source('code/functions.R')


## Bringing in flow data
## Verde flow data at Paulden 7/17/1963-2017, 54 years continuous
raw.flow <- read.csv("data/flowdata_Verde.csv") 

## str(flowdata)
head(raw.flow)
## FloodMag - magnitude of flood in cfs
## BaseDur - baseflow duration in days
## flooddate is peak dates of all floods (Oct 1 = 1) because using water year

str(raw.flow)

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

natural.flow <- natural.flow[,c(1,7,2:6)]
natural.flow
head(natural.flow)
apply(natural.flow[,-c(1:2)], 1, sum)
apply(natural.flow[,-c(1:2)], 2, sum)

natural.flow %>% mutate(sum = rowSums(.[,3:7]))

## medflood is never found with summer highflood, so I won't allow that combo.
## We therefore have nonevent, drought, medflood, and any combo of the two highfloods.
## i.e. each individually, and combined

## Making flow scenarios
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

drought.list

## Checking to see it's working.
## Should finish w/ 54 drought years and 0 others
plyr::ldply(drought.list, function(x) apply(x, 2, sum))

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

nonevent.list

## Checking to see it's working.
## Should finish w/ 54 nonevent years and 0 others
plyr::ldply(nonevent.list, function(x) apply(x, 2, sum))


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

medflood.list

## Checking to see it's working.
## Should finish w/ 54 medflood years and 0 others
plyr::ldply(medflood.list, function(x) apply(x, 2, sum))


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

SP_highflood.list

## Checking to see it's working.
## Should finish w/ 54 SP_highflood years and 0 others
plyr::ldply(SP_highflood.list, function(x) apply(x, 2, sum))


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

SU_highflood.list

## Checking to see it's working.
## Should finish w/ 54 SU_highflood years and 0 others
plyr::ldply(SU_highflood.list, function(x) apply(x, 2, sum))
 

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

SPSU_highflood.list

## Checking to see it's working.
## Should finish w/ 54 SP and SU_highflood years and 0 others
plyr::ldply(SPSU_highflood.list, function(x) apply(x, 2, sum))


all.scenarios.list <- c(drought.list,
                        nonevent.list,
                        medflood.list,
                        SP_highflood.list,
                        SU_highflood.list,
                        SPSU_highflood.list)

map(all.scenarios.list, head)
map(all.scenarios.list, class)
names(all.scenarios.list)


## Save natural flow through to 2017
saveRDS(natflow_full, 'data/natflow_full.rds')


## Save the objects as .rds files - then use loadRDS in other file. 
saveRDS(flowlist, 'data/flowlist.rds')

