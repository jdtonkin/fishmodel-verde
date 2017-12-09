## Verde Fish Model
## Jane Rogosch, Jono Tonkin, et al.
## 01-Mar-17
## 26-Jul-17
## 22-Aug-17
## 17-Oct-17

## Required libraries
library(ggplot2)
library(plyr)
library(tidyverse)
library(popbio)

## SETUP ----------------------------------------------------------------------------------

rm(list = ls()) # clearing the workspace 




## bringing in flow data
## Verde flow data at Paulden 7/17/1963-2017, 54 years continuous
flowdata <- read.csv("data/flowdata_Verde.csv") 

## str(flowdata)
## head(flowdata)
## FloodMag - magnitude of flood in cfs
## BaseDur - baseflow duration in days
## flooddate is peak dates of all floods (Oct 1 = 1) because using water year

count <- 54 #54 years in flow record, if count = 45 goes to 2008 
iterations <- 1000 # number of replicate projections to run (mid loop)


## Modifiers
modifiers <- read.csv('data/modifiers-all-spp.csv')

## adding 'Modifier' value from csv to 'Code' in csv
for(j in 1:length(modifiers[,1])) {
    nam <- paste(modifiers[j,4])
    assign(nam, modifiers[j,5]) # SHOULD BE [j,5] FOR REAL VALUES, [j,6] for null values!!!!
}

## Vital rates
## Baseline maturation probabiliity, aCACL3 (adult senescence rate)
## Background mortality
## Initial volume in grams in 100-m reach
## Fecundity based on year type and GSI
## Stage specific densities (ind./g)

vitalrates <- read.csv('data/vital-rates.csv') 

## assigning vital rate values from column 3 to 'code' in column 2
for(k in 1:length(vitalrates[,1])) {
    nam <- paste(vitalrates[k,2])
    assign(nam, vitalrates[k,3]) 
}

## Key ------------------------------------------------------------
## CACL (Catostomus clarki) – desert sucker 
## GIRO (Gila robusta) – roundtail chub
## LECY (Lepomis cyanellus) – green sunfish
## CAIN (Catostomus insignis) – sonora sucker
## MIDO (Micropterus dolomieu) - smallmouth bass
## CYLU (Cyprinella lutrensis) – red shiner
## AMNA (Ameiurus natalis) – yellow bullhead

## Average total volume of water per 100 m reach in m3: 307
## Average total fish biomass per 100 m reach in g: 4766
## Average total biomass Bonar 2004 in g/100m2: 606
## Max for a 100 m rech in Gibson samples (excluding GAAF): 6996

## vector of species names
sppnames <- c('CACL', 'GIRO', 'LECY', 'CAIN', 'MIDO', 'CYLU', 'AMNA')

K = 47660 # mean for 1-km reach across 6 replicate reaches

## Bunch of functions -----------------------------------------
## =====================
## XXXXX - add all these to source file and read in using 'source(functions)'
## =====================

## checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
    ifelse(x < 0, 0, x)
}

## FLOOD THRESHOLD FUNCTIONS --------------------------------------------------------------
## flowdata$FloodMag - vector containing peak flood magnitude 

## Magnitude of peak flow over which is considered a large flood event
SP_highfloodcutoff = 700 # SP: Jan 1 - Apr 30 (water_day 93-213). Floods (cfs) at 4-yr recurrence interval or greater (30 times median flow).
SU_highfloodcutoff = 200 # SU: May 1 - Sep 30 (Water_day 214-365). Floods (Cfs) at 4-year recurrence interval for summer
medfloodcutoff = 220 # SP Floods (cfs) at bankfull flood, 2.5-yr recurrence interval (10 times median flow)
## Non-event is same for drought (threshold below rather than above)
mindroughtlength = 40 # threshold length of low-flow in days (greater than 75th percentile of low flow duration) -- 
## determined as number of consecutive days with discharge (cfs) less than 22 cfs (25th percentile of flows) May 1 - Sep 30
## This threshold is compared to baseflow duration in any given year - a vector specified in csv (baseflow_dur).

## Convert peak discharge values into a vector of floods/no floods, 
SP_highflood_func <- function(x) {
    ifelse(x >= SP_highfloodcutoff, 1, 0)
}

SU_highflood_func <- function(x) {
    ifelse(x >= SU_highfloodcutoff, 1, 0)
}

medflood_func <- function(x) {
    ifelse(x >= medfloodcutoff & x < SP_highfloodcutoff, 1, 0)
}

nonevent_func <- function(Spfl, BD, Sufl) {
    ifelse(Spfl < medfloodcutoff & BD < mindroughtlength & Sufl < SU_highfloodcutoff, 1, 0)
}

drought_func <- function(Spfl, BD, Sufl) {
    ifelse(Spfl < medfloodcutoff & BD >= mindroughtlength & Sufl < SU_highfloodcutoff, 1, 0)
}


## Defining year types here 
SP_highflood <- SP_highflood_func(flowdata$SpFloodMag) 

SU_highflood <- SU_highflood_func(flowdata$SuFloodMag)

medflood <- medflood_func(flowdata$SpFloodMag)

drought <- drought_func(Spfl = flowdata$SpFloodMag, BD = flowdata$BaseDur, Sufl = flowdata$SuFloodMag)

nonevent <- nonevent_func(Spfl = flowdata$SpFloodMag, BD = flowdata$BaseDur, Sufl = flowdata$SuFloodMag) 


## ITERATION PARAMETERS -------------------------------------------------------------------
## Setting up arrays/vectors to fill with data from loops

## Mid loop details -----------------------------------------------------------------------
## 'iterations' - number of replicate flow sequences to run for averaging, SE, etc.

years <- as.character(seq(1964, 2017, by = 1))
stages <- as.character(c("S1", "S2", "S3"))
## CACLrep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
## GIROrep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
## LECYrep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
## CAINrep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
## MIDOrep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
## CYLUrep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
## AMNArep <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))  
Total.N <- array(0, dim = c(54, iterations), dimnames = list(years, 1:iterations))

## === replace total.n with a df instead. redo stuff following too. 
str(Total.N)
#### HERE -- need to replace the stuff below - e.g. caclrep etc.
## ========================
## ========================
## the following replaces the commented out region above
## ========================

## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
reparray <- array(0, dim = c(54, 3, iterations), dimnames = list(years, stages, 1:iterations))

## Repeating the array 7 times 
replist <- rep(list(reparray), 7)

## Assigning names to each array from sppnames vector
names(replist) <- sppnames

## Inner loop details ---------------------------------------------------------------------
## 'count' - number of years to project simulations (inner loop)

## Output of biomass and no. ind. for each age class for each year projected  
## An array with 3 columns (each stage class) and however many rows there are years projected 
## =====================
## need to replace all teh stuff below in the remaining code
## =====================
## CACLoutput.N <- array(0, dim = c(count, 3))
## GIROoutput.N <- array(0, dim = c(count, 3))
## LECYoutput.N <- array(0, dim = c(count, 3))
## CAINoutput.N <- array(0, dim = c(count, 3))
## MIDOoutput.N <- array(0, dim = c(count, 3))
## CYLUoutput.N <- array(0, dim = c(count, 3))
## AMNAoutput.N <- array(0, dim = c(count, 3))

## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
output.N.array <- array(0, dim = c(count, 3))

## Repeating the array 7 times 
output.N.list <- rep(list(output.N.array), 7)

## Assigning names to each array from sppnames vector
names(output.N.list) <- sppnames



## CACL.lambda <- array(0, dim = c(count,1))
## GIRO.lambda <- array(0, dim = c(count,1))
## LECY.lambda <- array(0, dim = c(count,1))
## CAIN.lambda <- array(0, dim = c(count,1))
## MIDO.lambda <- array(0, dim = c(count,1))
## CYLU.lambda <- array(0, dim = c(count,1))
## AMNA.lambda <- array(0, dim = c(count,1))

## Create a df to fill in w/ lambda values
lambda.df <- data.frame(matrix(ncol = 7, nrow = count))
names(lambda.df) <- sppnames

## Biomass 
## CYLUoutput.biom <- array(0, dim = c(count, 3))
## MIDOoutput.biom <- array(0, dim = c(count, 3))
## CAINoutput.biom <- array(0, dim = c(count, 3))
## LECYoutput.biom <- array(0, dim = c(count, 3))
## GIROoutput.biom <- array(0, dim = c(count, 3))
## CACLoutput.biom <- array(0, dim = c(count, 3))
## AMNAoutput.biom <- array(0, dim = c(count, 3))

## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
output.biom.array <- array(0, dim = c(count, 3))

## Repeating the array 7 times 
output.biom.list <- rep(list(output.biom.array), 7)

## Assigning names to each array from sppnames vector
names(output.biom.list) <- sppnames


## Total biomass as % of K 
## CACLbiomoutput <- numeric(length = count)
## GIRObiomoutput <- numeric(length = count)
## LECYbiomoutput <- numeric(length = count)
## CAINbiomoutput <- numeric(length = count)
## MIDObiomoutput <- numeric(length = count)
## CYLUbiomoutput <- numeric(length = count)
## AMNAbiomoutput <- numeric(length = count)

## Creating a list of 7 vectors to fill in. One for each spp. 
## Create a vector to be repeated
biomoutput.vector <- numeric(length = count)

## Repeating the vector 7 times 
biomoutput.list <- rep(list(biomoutput.vector), 7)

## Assigning names to each vector from sppnames vector
names(biomoutput.list) <- sppnames




## Flood and drought settings for each year projected into the future (i.e. 0 or 1) 
## SPhighfloodoutput <- numeric(length = count) 
## SUhighfloodoutput <- numeric(length = count)
## medfloodoutput <- numeric(length = count) 
## droughtoutput <- numeric(length = count) 
## noneventoutput <- numeric(length = count) 

## Create a data frame with 5 columns and 'count' rows to fill in with flow results
flowresults <- data.frame(matrix(ncol = 5, nrow = count))
names(flowresults) <- c('SPhighflood', 'SUhighflood', 'medflood', 'drought', 'nonevent')
flowresults


## Fecundities
## FCACLoutput <- numeric(length = count)
## FGIROoutput <- numeric(length = count)
## FLECYoutput <- numeric(length = count)
## FCAINoutput <- numeric(length = count)
## FMIDOoutput <- numeric(length = count)
## FCYLUoutput <- numeric(length = count)
## FAMNAoutput <- numeric(length = count)

## Creating a list of 7 vectors to fill in. One for each spp. 
## Create a vector to be repeated
fec.vector <- numeric(length = count)

## Repeating the vector 7 times 
fec.list <- rep(list(fec.vector), 7)

## Assigning names to each vector from sppnames vector
names(fec.list) <- sppnames


## Mid loop ###############################################################################
## Middle loop uses iterator "iter" to get "iterations" for suming S2 and S3
for(iter in 1:iterations) {

    ## USE THIS to examine different flow year types ++++++++++++++++++++++++++++++++++++++++++++++++++++  
    ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
                                         # # All 2010 SPflood + SUflood years 
                                         # z <- rep(47, 84)
    
                                         # # All Spflood 1993
                                         # z <- rep(30, 84)
    
                                         # # All drought Y2K
                                         # z <- rep(37, 84)
    
                                         # # Nonevent 1985
                                         # z <- rep(22, 84)
    
                                         # # SUflood
                                         # z <- rep(21, 84)
    
                                         # # Medflood
                                         # z <- rep(25, 84)
    
    ## Need to read in initial biom every time so starting biomass is reset each iteration
    ##  N gives the total number of individuals for each age class.
    ##   Initially here, this is found by multiplying the number of g occupied by a given class
    ##   by the density per g
    ##   biom = g/m3
    ##   den = indiv/g
    
    ## To have different initial starting population sizes for each iteration, taking biom of stage 3 from negative binomial
    ##   distribution, where the parameter (lambda = mean) and K (dispersion) is calculated 
    ##   from mean and variance in abundance across seven sites in Verde River from 94-08, and scaled to biomass from Gibson 2012 survey 
    ##   in file "Rinne Verde River Data 1994-2008-.xlsx"
    
    biomCACL <- c(biomCACL1, biomCACL2, rnbinom(1, size = 1.52, mu = 5284)) 
    
    biomGIRO <- c(biomGIRO1, biomGIRO2, rnbinom(1, 0.44, mu = 2376)) 
    
    biomLECY <- c(biomLECY1, biomLECY2, rnbinom(1, 0.34, mu = 164)) 
    
    biomCAIN <- c(biomCAIN1, biomCAIN2, rnbinom(1, 1.33, mu = 34068)) 
    
    biomMIDO <- c(biomMIDO1, biomMIDO2, rnbinom(1, 0.66, mu = 4202)) 
    
    biomCYLU <- c(biomCYLU1, biomCYLU2, rnbinom(1, 1.78, mu = 238)) 
    
    biomAMNA <- c(biomAMNA1, biomAMNA2, rnbinom(1, 0.36, mu = 1306))
    
    ## Inner loop #############################################################################
    for(i in 1:count) {

        ## CHANGE WHAT 'y' IS TO SIMULATE DIFFERENT FLOW REGIMES ACROSS THE 54 YEARS 
        y = i # follow flow record sequence
        ##y = z[i] # to examine specific flow year type defined by z above

        ## y is directly taken from flow vector. 

        ## VITAL RATE DEFINITIONS: Desert Sucker  -----------------------------------------------------
        ## G is prob. of transition to next stage
        ## P is prob. of remaining in that stage
        ## Baseline mortality vital rate (from file object: 'vitalrate') is multiplied by modifier (from file object: 'modifiers') based on yeartype as specified above

        ## Natives
        ## Desert Sucker - 
        ## YOY (GCACL1) survival and recruitment depends on spring flows (J is modifier for stage 1, A is modifier for stage 2&3)
        ## GCACL1 <- aCACL1 * denCACLJ * (1/denCACL2) *
        ##     (1 - (SP_highflood[y] * STMortCACL * CACL_J_HF)) *
        ##     (1 - (SU_highflood[y] * STMortCACL * CACL_J_NE)) *
        ##     (1 - (medflood[y] * STMortCACL * CACL_J_MF)) *
        ##     (1 - (nonevent[y] * STMortCACL* CACL_J_NE))  *
        ##     (1 - (drought[y] * STMortCACL* CACL_J_DR)) 
        
        ## GCACL2 <- aCACL2 * denCACL2 * (1/denCACL3) *
        ##     (1 - (SP_highflood[y] * STMortCACL* CACL_A_HF)) *
        ##     (1 - (SU_highflood[y] * STMortCACL * CACL_A_NE)) *
        ##     (1 - (medflood[y] * STMortCACL * CACL_A_MF))*
        ##     (1 - (nonevent[y] * STMortCACL * CACL_A_NE))  *
        ##     (1 - (drought[y] * STMortCACL * CACL_A_DR))  

        ## PCACL3 <- (1 - aCACL3) *
        ##     (1 - (SP_highflood[y] * STMortCACL * CACL_A_HF)) *
        ##     (1 - (SU_highflood[y] * STMortCACL * CACL_A_NE)) *
        ##     (1 - (medflood[y] * STMortCACL * CACL_A_MF))  *
        ##     (1 - (nonevent[y] * STMortCACL * CACL_A_NE))  *
        ##     (1 - (drought[y] * STMortCACL * CACL_A_DR)) 
        

        ## ## Chub
        ## GGIRO1 <- aGIRO1 * denGIROJ * (1/denGIRO2) *
        ##     (1 - (SP_highflood[y] * STMortGIRO * GIRO_J_HF)) *
        ##     (1 - (SU_highflood[y] * STMortGIRO * GIRO_J_NE)) *
        ##     (1 - (medflood[y] * STMortGIRO * GIRO_J_MF)) *
        ##     (1 - (nonevent[y] * STMortGIRO * GIRO_J_NE))  *
        ##     (1 - (drought[y] * STMortGIRO * GIRO_J_DR)) 
        
        ## GGIRO2 <- aGIRO2 * denGIRO2 * (1/denGIRO3) *
        ##     (1 - (SP_highflood[y] * STMortGIRO)) * GIRO_A_HF *
        ##     (1 - (SU_highflood[y] * STMortGIRO * GIRO_A_NE)) *
        ##     (1 - (medflood[y] * STMortGIRO * GIRO_A_MF)) *
        ##     (1 - (nonevent[y] * STMortGIRO * GIRO_A_NE)) *
        ##     (1 - (drought[y] * STMortGIRO * GIRO_A_DR)) 

        ## PGIRO3 <- (1 - aGIRO3) *
        ##     (1 - (SP_highflood[y] * STMortGIRO * GIRO_A_HF)) *
        ##     (1 - (SU_highflood[y] * STMortGIRO * GIRO_A_NE)) *
        ##     (1 - (medflood[y] * STMortGIRO * GIRO_A_MF)) *
        ##     (1 - (nonevent[y] * STMortGIRO * GIRO_A_NE)) *
        ##     (1 - (drought[y] * STMortGIRO * GIRO_A_DR))
        
        ## ## Sonora sucker
        ## GCAIN1 <- aCAIN1 * denCAINJ * (1/denCAIN2) *
        ##     (1 - (SP_highflood[y] * STMortCAIN * CAIN_J_HF)) *
        ##     (1 - (SU_highflood[y] * STMortCAIN * CAIN_J_NE)) *
        ##     (1 - (medflood[y] * STMortCAIN * CAIN_J_MF)) *
        ##     (1 - (nonevent[y] * STMortCAIN * CAIN_J_NE)) *
        ##     (1 - (drought[y] * STMortCAIN * CAIN_J_DR))
        
        ## GCAIN2 <- aCAIN2 * denCAIN2 * (1/denCAIN3) *
        ##     (1 - (SP_highflood[y] * STMortCAIN * CAIN_A_HF)) *
        ##     (1 - (SU_highflood[y] * STMortCAIN * CAIN_A_NE)) *
        ##     (1 - (medflood[y] * STMortCAIN * CAIN_A_MF)) *
        ##     (1 - (nonevent[y] * STMortCAIN * CAIN_A_NE)) *
        ##     (1 - (drought[y] * STMortCAIN * CAIN_A_DR))
        
        ## PCAIN3 <- (1 - aCAIN3) *
        ##     (1 - (SP_highflood[y] * STMortCAIN * CAIN_A_HF)) *
        ##     (1 - (SU_highflood[y] * STMortCAIN * CAIN_A_NE)) *
        ##     (1 - (medflood[y] * STMortCAIN * CAIN_A_MF)) *
        ##     (1 - (nonevent[y] * STMortCAIN * CAIN_A_NE)) *
        ##     (1 - (drought[y] * STMortCAIN * CAIN_A_DR))
        
        ## ## Non-natives
        ## ## Green sunfish - note the "NN"
        ## GLECY1 <- aLECY1 * denLECYJ * (1/denLECY2) *
        ##     (1 - (SU_highflood[y] * STMortLECY * LECY_J_HF)) *
        ##     (1 - (SP_highflood[y] * STMortLECY * LECY_J_MF)) *
        ##     (1 - (medflood[y] * STMortLECY * LECY_J_MF)) *
        ##     (1 - (nonevent[y] * STMortLECY * LECY_J_NE)) *
        ##     (1 - (drought[y] * STMortLECY * LECY_J_DR))
        
        ## GLECY2 <- aLECY2 * denLECY2 * (1/denLECY3) *
        ##     (1 - (SU_highflood[y] * STMortLECY * LECY_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortLECY * LECY_A_MF)) *
        ##     (1 - (medflood[y] * STMortLECY * LECY_A_MF)) *
        ##     (1 - (nonevent[y] * STMortLECY * LECY_A_NE)) *
        ##     (1 - (drought[y] * STMortLECY * LECY_A_DR))

        ## PLECY3 <- (1 - aLECY3) *
        ##     (1 - (SU_highflood[y] * STMortLECY * LECY_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortLECY * LECY_A_MF)) *
        ##     (1 - (medflood[y] * STMortLECY * LECY_A_MF)) *
        ##     (1 - (nonevent[y] * STMortLECY * LECY_A_NE)) *
        ##     (1 - (drought[y] * STMortLECY * LECY_A_DR))
        
        ## ## Smallmouth bass
        ## GMIDO1 <- aMIDO1 * denMIDOJ * (1/denMIDO2) *
        ##     (1 - (SU_highflood[y] * STMortMIDO * MIDO_J_HF)) *
        ##     (1 - (SP_highflood[y] * STMortMIDO * MIDO_J_MF)) *
        ##     (1 - (medflood[y] * STMortMIDO * MIDO_J_MF)) *
        ##     (1 - (nonevent[y] * STMortMIDO * MIDO_J_NE)) *
        ##     (1 - (drought[y] * STMortMIDO * MIDO_J_DR))
        
        ## GMIDO2 <- aMIDO2 * denMIDO2 * (1/denMIDO3) *
        ##     (1 - (SU_highflood[y] * STMortMIDO * MIDO_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortMIDO * MIDO_A_MF)) *
        ##     (1 - (medflood[y] * STMortMIDO * MIDO_A_MF)) *
        ##     (1 - (nonevent[y] * STMortMIDO * MIDO_A_NE)) *
        ##     (1 - (drought[y] * STMortMIDO * MIDO_A_DR))
        
        ## PMIDO3 <- (1 - aMIDO3) *
        ##     (1 - (SU_highflood[y] * STMortMIDO * MIDO_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortMIDO * MIDO_J_MF)) *
        ##     (1 - (medflood[y] * STMortMIDO * MIDO_A_MF)) *
        ##     (1 - (nonevent[y] * STMortMIDO * MIDO_A_NE)) *
        ##     (1 - (drought[y] * STMortMIDO * MIDO_A_DR))
        
        ## ## Red shiner
        ## GCYLU1 <- aCYLU1 * denCYLUJ * (1/denCYLU2) *
        ##     (1 - (SU_highflood[y] * STMortCYLU * CYLU_J_HF)) *
        ##     (1 - (SP_highflood[y] * STMortCYLU * CYLU_J_MF)) *
        ##     (1 - (medflood[y] * STMortCYLU * CYLU_J_MF)) *
        ##     (1 - (nonevent[y] * STMortCYLU * CYLU_J_NE)) *
        ##     (1 - (drought[y] * STMortCYLU * CYLU_J_DR))
        
        ## GCYLU2 <- aCYLU2 * denCYLU2 * (1/denCYLU3) *
        ##     (1 - (SU_highflood[y] * STMortCYLU * CYLU_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortCYLU * CYLU_A_MF)) *
        ##     (1 - (medflood[y] * STMortCYLU * CYLU_A_MF)) *
        ##     (1 - (nonevent[y] * STMortCYLU * CYLU_A_NE)) *
        ##     (1 - (drought[y] * STMortCYLU * CYLU_A_DR))
        
        ## PCYLU3 <- (1 - aCYLU3) *
        ##     (1 - (SU_highflood[y] * STMortCYLU * CYLU_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortCYLU * CYLU_A_MF)) *
        ##     (1 - (medflood[y] * STMortCYLU * CYLU_A_MF)) *
        ##     (1 - (nonevent[y] * STMortCYLU * CYLU_A_NE)) *
        ##     (1 - (drought[y] * STMortCYLU * CYLU_A_DR))
        
        ## ## Yellow bullhead
        ## GAMNA1 <- aAMNA1 * denAMNAJ * (1/denAMNA2) *
        ##     (1 - (SU_highflood[y] * STMortAMNA * AMNA_J_HF)) *
        ##     (1 - (SP_highflood[y] * STMortAMNA * AMNA_J_MF)) *
        ##     (1 - (medflood[y] * STMortAMNA * AMNA_J_MF)) *
        ##     (1 - (nonevent[y] * STMortAMNA * AMNA_J_NE)) *
        ##     (1 - (drought[y] * STMortAMNA * AMNA_J_DR))
        
        ## GAMNA2 <- aAMNA2 * denAMNA2 * (1/denAMNA3) *
        ##     (1 - (SU_highflood[y] * STMortAMNA * AMNA_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortAMNA * AMNA_A_MF)) *
        ##     (1 - (medflood[y] * STMortAMNA * AMNA_A_MF)) *
        ##     (1 - (nonevent[y] * STMortAMNA * AMNA_A_NE)) *
        ##     (1 - (drought[y] * STMortAMNA * AMNA_A_DR))
        
        ## PAMNA3 <- (1 - aAMNA3) *
        ##     (1 - (SU_highflood[y] * STMortAMNA * AMNA_A_HF)) *
        ##     (1 - (SP_highflood[y] * STMortAMNA * AMNA_A_MF)) *
        ##     (1 - (medflood[y] * STMortAMNA * AMNA_A_MF)) *
        ##     (1 - (nonevent[y] * STMortAMNA * AMNA_A_NE)) *
        ##     (1 - (drought[y] * STMortAMNA * AMNA_A_DR))




###### PROBLEMS ########
        #### Janes highflood years have both MF and NE modifiers at the end - I've made the below 'HF'

        ## Stage 1 - G
        for(nm in sppnames) {
            assign(paste0('G', nm, '1'),

                   get(paste0('a', nm, '1')) * get(paste0('den', nm, 'J')) * (1 / get(paste0('den', nm, '2')))  *
                   (1 - (SU_highflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_J_HF')))) *
                   (1 - (SP_highflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_J_HF')))) *
                   (1 - (medflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_J_MF')))) *
                   (1 - (nonevent[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_J_NE')))) *
                   (1 - (drought[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_J_DR'))))
                   )          
        }        

        ## Stage 2 - G
        for(nm in sppnames) {
            assign(paste0('G', nm, '2'),

                   get(paste0('a', nm, '2')) * get(paste0('den', nm, '2')) * (1 / get(paste0('den', nm, '3')))  *
                   (1 - (SU_highflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_HF')))) *
                   (1 - (SP_highflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_HF')))) *
                   (1 - (medflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_MF')))) *
                   (1 - (nonevent[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_NE')))) *
                   (1 - (drought[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_DR'))))
                   )          
        }        

        ## Stage 3 - P
        for(nm in sppnames) {
            assign(paste0('P', nm, '3'),

                   (1 - get(paste0('a', nm, '3'))) *
                   (1 - (SU_highflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_HF')))) *
                   (1 - (SP_highflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_HF')))) *
                   (1 - (medflood[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_MF')))) *
                   (1 - (nonevent[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_NE')))) *
                   (1 - (drought[y] * get(paste0('STMort', nm)) * get(paste0(nm, '_A_DR'))))
                   )          
        }          
        
        ## Total grams occupied after year -----------------------------------------------
       
        totbiom <- ldply(sppnames, function(x)
            get(paste0('biom', x))[1] +
            get(paste0('biom', x))[2] +
            get(paste0('biom', x))[3]) %>%
            sum
        
        ## totbiom.CACL <- 
        ##     biomCACL[1] + #
        ##     biomCACL[2] + 
        ##     biomCACL[3] 

        ## totbiom.GIRO <- 
        ##     biomGIRO[1] + #
        ##     biomGIRO[2] + 
        ##     biomGIRO[3] 

        ## totbiom.LECY <- 
        ##     biomLECY[1] + #
        ##     biomLECY[2] + 
        ##     biomLECY[3] 

        ## totbiom.CAIN <- 
        ##     biomCAIN[1] + #
        ##     biomCAIN[2] + 
        ##     biomCAIN[3] 

        ## totbiom.MIDO <- 
        ##     biomMIDO[1] + #
        ##     biomMIDO[2] + 
        ##     biomMIDO[3] 

        ## totbiom.CYLU <- 
        ##     biomCYLU[1] + #
        ##     biomCYLU[2] + 
        ##     biomCYLU[3] 

        ## totbiom.AMNA <- 
        ##     biomAMNA[1] + #
        ##     biomAMNA[2] + 
        ##     biomAMNA[3] 

        
### Carrying capacity (K) is limiting spawning of all species based on the total biomass occupied at the end of the previous year. 
        ## i.e. if above K, no spp spawn in that year. If spawning occurs, they all spawn.

        ## Some slight differences in fecund so can't loop/lapply
        
        ## POTENTIAL CACL FECUNDITY ---------------------------------------------------------
        FCACL2 <- ((0.5 * GSI.CACL * (1 - S0MortCACL)) * checkpos((K - totbiom)/K)) *
            denCACL1 * (1/denCACLJ)
        
        FCACL3 <- ((0.5 * GSI.CACL * (1 - S0MortCACL)) * checkpos((K - totbiom)/K)) *
            denCACL1 * (1/denCACLJ)

        ## POTENTIAL GIRO FECUNDITY ---------------------------------------------------------
        FGIRO2 <- ((0.5 * GSI.GIRO * (1 - S0MortGIRO)) * checkpos((K - totbiom)/K)) * 
            denGIRO1 * (1/denGIROJ)
        
        FGIRO3 <- ((0.5 * GSI.GIRO * (1 - S0MortGIRO)) * checkpos((K - totbiom)/K)) *
            denGIRO1 * (1/denGIROJ)

        ## POTENTIAL CAIN FECUNDITY ---------------------------------------------------------
        FCAIN3 <- ((0.5 * GSI.CAIN * (1 - S0MortCAIN)) * checkpos((K - totbiom)/K)) *
            denCAIN1 * (1/denCAINJ)

        ## POTENTIAL LECY FECUNDITY ---------------------------------------------------------
        FLECY2 <- ((0.5 * GSI.LECY * (1 - S0MortLECY)) * checkpos((K - totbiom)/K)) *
            denLECY1 * (1/denLECYJ)
        
        FLECY3 <- ((0.5 * GSI.LECY * (1 - S0MortLECY)) * checkpos((K - totbiom)/K)) * 
            denLECY1 * (1/denLECYJ)
        
        ## POTENTIAL MIDO FECUNDITY ---------------------------------------------------------
        FMIDO3 <- ((0.5 * GSI.MIDO * (1 - S0MortMIDO)) * checkpos((K - totbiom)/K)) *
            denMIDO1 * (1/denMIDOJ)

        ## POTENTIAL CYLU FECUNDITY ---------------------------------------------------------
        ## because they are serial spawners, they are allowed to spawn twice a season in stage 2 and 3
        FCYLUJ <- ((0.5 * GSI.CYLU * (1 - S0MortCYLU)) * checkpos((K - totbiom)/K)) * 
            denCYLU1 * (1/denCYLUJ)
        
        FCYLU2 <- ((0.5 * 2 * GSI.CYLU * (1 - S0MortCYLU)) * checkpos((K - totbiom)/K)) *
            denCYLU1 * (1/denCYLUJ)
        
        FCYLU3 <- ((0.5 * 2 * GSI.CYLU * (1 - S0MortCYLU)) * checkpos((K - totbiom)/K)) *
            denCYLU1 * (1/denCYLUJ)

        ## POTENTIAL AMNA FECUNDITY ---------------------------------------------------------
        FAMNA3 <- ((0.5 * GSI.AMNA * (1 - S0MortAMNA)) * checkpos((K - totbiom)/K)) *
            denAMNA1 * (1/denAMNAJ)


        ## K --------------------------------------------------------------------------------------
        ## CACL
        ## gives total CACL biomass (g) as a percentage of K; 
        ## KCACL <- 100 * (biomCACL[1] + 
        ##                 biomCACL[2] + 
        ##                 biomCACL[3])/K

        ## ## GIRO
        ## ## gives total GIRO biomass as a percentage of K; 
        ## KGIRO <- 100 * (biomGIRO[1] + 
        ##                 biomGIRO[2] + 
        ##                 biomGIRO[3])/K

        ## ## LECY
        ## ## gives total LECY biomass as a percentage of K; 
        ## KLECY <- 100 * (biomLECY[1] + 
        ##                 biomLECY[2] + 
        ##                 biomLECY[3])/K

        ## ## CAIN
        ## ## gives total CAIN biomass as a percentage of K; 
        ## KCAIN <- 100 * (biomCAIN[1] + 
        ##                 biomCAIN[2] + 
        ##                 biomCAIN[3])/K

        ## ## MIDO
        ## ## gives total MIDO biomass as a percentage of K; 
        ## KMIDO <- 100 * (biomMIDO[1] + 
        ##                 biomMIDO[2] + 
        ##                 biomMIDO[3])/K

        ## ## CYLU
        ## ## gives total CYLU biomass as a percentage of K; 
        ## KCYLU <- 100 * (biomCYLU[1] + 
        ##                 biomCYLU[2] + 
        ##                 biomCYLU[3])/K

        ## ## AMNA
        ## ## gives total CYLU biomass as a percentage of K; 
        ## KAMNA <- 100 * (biomAMNA[1] + 
        ##                 biomAMNA[2] + 
        ##                 biomAMNA[3])/K

        ## Calculating the percentage of K occupied and addting to KCACL etc
        for(nm in sppnames) {
            assign(paste0('K', nm), 100 * (
                get(paste0('biom', nm))[1] +   
                get(paste0('biom', nm))[2] +
                get(paste0('biom', nm))[3])/K)                                  
        }
        
        
        ## TRANSITION MATRICES --------------------------------------------------------------------

        ## TRANSITION MATRIX FOR CACL -------------------------------------------------------
        ACACL1 <- c(0, FCACL2, FCACL3)
        ACACL2 <- c(GCACL1, 0, 0)
        ACACL3 <- c(0, GCACL2, PCACL3)
        ## Matrix
        ## ACACL <- rbind(ACACL1, ACACL2, ACACL3)

        ## TRANSITION MATRIX FOR GIRO -------------------------------------------------------
        AGIRO1 <- c(0, FGIRO2, FGIRO3)
        AGIRO2 <- c(GGIRO1, 0, 0)
        AGIRO3 <- c(0, GGIRO2, PGIRO3)
        ## Matrix
        ## AGIRO <- rbind(AGIRO1, AGIRO2, AGIRO3)

        ## TRANSITION MATRIX FOR LECY -------------------------------------------------------
        ALECY1 <- c(0, FLECY2, FLECY3)
        ALECY2 <- c(GLECY1, 0, 0)
        ALECY3 <- c(0, GLECY2, PLECY3)
        ## Matrix
        ## ALECY <- rbind(ALECY1, ALECY2, ALECY3)

        ## TRANSITION MATRIX FOR CAIN -------------------------------------------------------
        ACAIN1 <- c(0, 0, FCAIN3)
        ACAIN2 <- c(GCAIN1, 0, 0)
        ACAIN3 <- c(0, GCAIN2, PCAIN3)
        ## Matrix
        ## ACAIN <- rbind(ACAIN1, ACAIN2, ACAIN3)

        ## TRANSITION MATRIX FOR MIDO -------------------------------------------------------
        AMIDO1 <- c(0, 0, FMIDO3)
        AMIDO2 <- c(GMIDO1, 0, 0)
        AMIDO3 <- c(0, GMIDO2, PMIDO3)
        ## Matrix
        ## AMIDO <- rbind(AMIDO1, AMIDO2, AMIDO3)

        ## TRANSITION MATRIX FOR CYLU -------------------------------------------------------
        ACYLU1 <- c(FCYLUJ, FCYLU2, FCYLU3)
        ACYLU2 <- c(GCYLU1, 0, 0)
        ACYLU3 <- c(0, GCYLU2, PCYLU3)
        ## Matrix
        ## ACYLU <- rbind(ACYLU1, ACYLU2, ACYLU3)

        ## TRANSITION MATRIX FOR AMNA -------------------------------------------------------
        AAMNA1 <- c(0, 0, FAMNA3)
        AAMNA2 <- c(GAMNA1, 0, 0)
        AAMNA3 <- c(0, GAMNA2, PAMNA3)
        ## Matrix
        ## AAMNA <- rbind(AAMNA1, AAMNA2, AAMNA3)


        ## AAMNA1 <- c('a', 'b', 'c')
        ## AAMNA2 <- c('d', 'e', 'f')
        ## AAMNA3 <- c('g', 'h', 'i')

        ## ACACL3 <- c('a', 'b', 'c')
        ## ACACL2 <- c('d', 'e', 'f')
        ## ACACL1 <- c('g', 'h', 'i')

        ## spnmtest <- c('AMNA', 'CACL')

        ## rbinding the vectors from above into transition matrices
        for(nm in sppnames) {
            assign(paste0('A', nm), rbind(
                                        get(paste0('A', nm, '1')),
                                        get(paste0('A', nm, '2')),
                                        get(paste0('A', nm, '3'))
                                    ))
        }
        
        ## COMPILING OUTPUTS ----------------------------------------------------------------------


        ## Lambda values
        ## Filling in the df with lambda values for each species and each year
        ## Species as columns, years as rows
        ## This applies 'lambda(ACACL)' etc and adds to correct column each 'i' value (year)
        lambda.df[i,] <- sapply(mget(paste0('A', names(lambda.df))), lambda)

        ## Fecundity values
        ## Cant loop or anything as different for diff spp
        fec.list$CACL[i] <- FCACL3 + FCACL2
        fec.list$GIRO[i] <- FGIRO3 + FGIRO2
        fec.list$LECY[i] <- FLECY3 + FLECY2
        fec.list$CAIN[i] <- FCAIN3
        fec.list$MIDO[i] <- FMIDO3
        fec.list$CYLU[i] <- FCYLU3 + FCYLU2 + FCYLUJ
        fec.list$AMNA[i] <- FAMNA3


        ## biomass values into each df/array in the list
        for(nm in sppnames) {
            output.biom.list[[nm]][i,1:3] <- get(paste0('biom', nm))
        }
        

        ## N values into each df/array in the list
        
        for(nm in sppnames) {
            output.N.list[[nm]][i,1:3] <- c(get(paste0('biom', nm))[1] * get(paste0('den', nm, 'J')),
                                            get(paste0('biom', nm))[2] * get(paste0('den', nm, '2')),
                                            get(paste0('biom', nm))[3] * get(paste0('den', nm, '3')))
        }
        
        

        #map2(.x = output.biom.list, .y = names(output.biom.list), ~ get(paste0('biom', .y)))

        
        ## CACL
        ## CACLoutput.biom[i,1:3] <- biomCACL # array of biomass of each age class for each yr projected. biomCACL = total biomass for each age class
        ## ##CACLbiomoutput[i] <- KCACL # total biomass as % of K; this is in g
        ## output.N.list$CACL[i, 1:3] <- 
        ##     c(biomCACL[1] * denCACLJ,
        ##       biomCACL[2] * denCACL2,
        ##       biomCACL[3] * denCACL3)
        ## ##CACL.lambda[i] <- lambda(ACACL)
        
        ## GIRO
        ## GIROoutput.biom[i,1:3] <- biomGIRO # array of biomass of each age class for each yr projected. biomGIRO = total biomass for each age class
        ## ##GIRObiomoutput[i] <- KGIRO # total biomass as % of K; this is in g
        ## output.N.list$GIRO[i, 1:3] <- 
        ##     c(biomGIRO[1] * denGIROJ,
        ##       biomGIRO[2] * denGIRO2,
        ##       biomGIRO[3] * denGIRO3)
        ## ##GIRO.lambda[i] <- lambda(AGIRO)   

        ## ## LECY
        ## LECYoutput.biom[i,1:3] <- biomLECY # array of biomass of each age class for each yr projected. biomLECY = total biomass for each age class
        ## ##LECYbiomoutput[i] <- KLECY # total biomass as % of K; this is in g 
        ## output.N.list$LECY[i, 1:3] <- 
        ##     c(biomLECY[1] * denLECYJ,
        ##       biomLECY[2] * denLECY2,
        ##       biomLECY[3] * denLECY3)
        ## ##LECY.lambda[i] <- lambda(ALECY)

        ## ## CAIN
        ## CAINoutput.biom[i,1:3] <- biomCAIN # array of biomass of each age class for each yr projected. biomCAIN = total biomass for each age class
        ## ##CAINbiomoutput[i] <- KCAIN # total biomass as % of K; this is in g 
        ## output.N.list$CAIN[i, 1:3] <- 
        ##     c(biomCAIN[1] * denCAINJ,
        ##       biomCAIN[2] * denCAIN2,
        ##       biomCAIN[3] * denCAIN3)
        ## ##CAIN.lambda[i] <- lambda(ACAIN)

        ## ## MIDO
        ## MIDOoutput.biom[i,1:3] <- biomMIDO # array of biomass of each age class for each yr projected. biomMIDO = total biomass for each age class
        ## ##MIDObiomoutput[i] <- KMIDO # total biomass as % of K; this is in g 
        ## output.N.list$MIDO[i, 1:3] <- 
        ##     c(biomMIDO[1] * denMIDOJ,
        ##       biomMIDO[2] * denMIDO2,
        ##       biomMIDO[3] * denMIDO3)
        ## ##MIDO.lambda[i] <- lambda(AMIDO)

        ## ## CYLU
        ## CYLUoutput.biom[i,1:3] <- biomCYLU # array of biomass of each age class for each yr projected. biomCYLU = total biomass for each age class
        ## ##CYLUbiomoutput[i] <- KCYLU # total biomass as % of K; this is in g 
        ## output.N.list$CYLU[i, 1:3] <- 
        ##     c(biomCYLU[1] * denCYLUJ,
        ##       biomCYLU[2] * denCYLU2,
        ##       biomCYLU[3] * denCYLU3)
        ## ##CYLU.lambda[i] <- lambda(ACYLU)

        ## ## AMNA
        ## AMNAoutput.biom[i,1:3] <- biomAMNA # array of biomass of each age class for each yr projected. biomAMNA = total biomass for each age class
        ## ##AMNAbiomoutput[i] <- KAMNA # total biomass as % of K; this is in g 
        ## output.N.list$AMNA[i, 1:3] <- 
        ##     c(biomAMNA[1] * denAMNAJ,
        ##       biomAMNA[2] * denAMNA2,
        ##       biomAMNA[3] * denAMNA3)
        ## ##AMNA.lambda[i] <- lambda(AAMNA)

        ## Records flood settings of each particular projected year (0 for nonflood, 1 for flood)
        flowresults$SPhighflood[i] <- SP_highflood[y]
        flowresults$SUhighflood[i] <- SU_highflood[y]
        flowresults$medflood[i] <- medflood[y]
        flowresults$nonevent[i] <- nonevent[y]
        flowresults$drought[i] <- drought[y]

        
        ## MATRIX MULTIPLICATION ------------------------------------------------------------------
        ## can include rescue function for each with 0.5 chance of reach being colonized by 2 individuals
        ## ## CACL
        ## biomCACL <- ACACL %*% biomCACL # ACACL is transition matrix, biomCACL = total biomass for each age class,

        ## ## GIRO
        ## biomGIRO <- AGIRO %*% biomGIRO # AGIRO is transition matrix, biomGIRO = total biomass for each age class

        ## ## LECY
        ## biomLECY <- ALECY %*% biomLECY  # ALECY is transition matrix, biomLECY = total biomass for each age class

        ## ## CAIN
        ## biomCAIN <- ACAIN %*% biomCAIN # ACAIN is transition matrix, biomCAIN = total biomass for each age class

        ## ## MIDO
        ## biomMIDO <- AMIDO %*% biomMIDO # AMIDO is transition matrix, biomMIDO = total biomass for each age class

        ## ## CYLU
        ## biomCYLU <- ACYLU %*% biomCYLU # ACYLU is transition matrix, biomCYLU = total biomass for each age class

        ## ## AMNA
        ## biomAMNA <- AAMNA %*% biomAMNA # AAMNA is transition matrix, biomAMNA = total biomass for each age class


       
        ## biomCACL <- c(1,2,3)
        ## biomGIRO <- c(1,2,3)
        ## biomLECY <- c(1,2,3)
        ## biomCAIN <- c(1,2,3)
        ## biomMIDO <- c(1,2,3)
        ## biomCYLU <- c(1,2,3)
        ## biomAMNA <- c(1,2,3)
        
        ## ## Getting matrices just to test
        ## ACACL <- matrix(1:9,3)
        ## AGIRO <- matrix(1:3,3,3)
        ## ALECY <- matrix(1:4,3,3)
        ## ACAIN <- matrix(1:5,3,3)
        ## AMIDO <- matrix(1:6,3,3)
        ## ACYLU <- matrix(1:7,3,3)
        ## AAMNA <- matrix(1:8,3,3)

        ## Matrix multiplication loop - essentially ==
        ## biomAMNA <- AAMNA %*% biomAMNA # AAMNA is transition matrix, biomAMNA = total biomass for each age class
        for(nm in sppnames) {
             assign(paste0('biom', nm), get(paste0('A', nm)) %*% get(paste0('biom', nm)))
        }
        
        
    } # End of inner loop ####################################################################

    ## Mean values for each iteration run over each sequence of years
    ## GESPoutput.N[31:45, 2:3] to compare 94-08 with stage 2 and stage 3 added together
    ## replist$CACL[,,1] <- output.N.list$CACL
    ## replist$GIRO[,,1] <- output.N.list$GIRO
    ## replist$LECY[,,iter] <- output.N.list$LECY
    ## replist$CAIN[,,iter] <- output.N.list$CAIN
    ## replist$MIDO[,,iter] <- output.N.list$MIDO
    ## replist$CYLU[,,iter] <- output.N.list$CYLU
    ## replist$AMNA[,,iter] <- output.N.list$AMNA

    ## In a loop
    for(nm in sppnames) {
        replist[[nm]][,,iter] <- output.N.list[[nm]]    
    }
    
    ## Total.N[,iter] <- apply(
    ##     cbind(output.N.list$CACL[,2:3],
    ##           output.N.list$GIRO[,2:3],
    ##           output.N.list$LECY[,2:3],
    ##           output.N.list$CAIN[,2:3],
    ##           output.N.list$MIDO[,2:3],
    ##           output.N.list$CYLU[,2:3],
    ##           output.N.list$AMNA[,2:3]), 1, sum)

        
    ## Caculating Total.N for each year, and adding it to total.N data frame with however many iterations run.
    ## Total does not incl. juveniles.
    ## map is purrr version of lapply. Can pass fn using ~ and .x instead of function(x) x
    ## Gets list output of stages 2:3 for ea spp, then cbinds them all together, then caclcs sum. 
    Total.N[,iter] <- map(output.N.list, ~ .x[,2:3]) %>%
        do.call('cbind', .) %>%
        apply(1, sum)


} # End of mid loop ######################################################################
## ########################################################################################



## OUTPUTS --------------------------------------------------------------------------------
## ########################################################################################

## FINAL iteration data to examine plots ---------------------------------------
## Compiling abundance and biomass outputs into single dfs ----------------------------------------------------

### OLD
## Function to make longform dataframes of biomass and abundance  
## make.dat <- function(spec, biomN){

##     datname <- get(paste0(spec, 'output.', biomN))
    
##     dat <- as.data.frame(datname) %>%
##         rename(S1 = V1, S2 = V2, S3 = V3) %>%
##         mutate(rep = row.names(.)) %>%
##         gather(stage, val, -rep) %>%
##         mutate(spp = spec)            

##     return(dat)
## }

## Using ldply to call sppnames and create a new df element in list based on make.dat fn.
## Then convert back to one full df. 

## Biom
#ALLoutput.biom.DF <- ldply(sppnames, function(x) make.dat(x, biomN = 'biom')) %>%
 #   rename(g = val) 



### NEW

ALLoutput.biom.DF <- ldply(output.biom.list, function(x) {
    as.data.frame(x) %>%
        rename(S1 = V1, S2 = V2, S3 = V3) %>%
        mutate(rep = row.names(.)) %>%
        gather(stage, val, -rep) 
}) %>%
    rename(spp = `.id`, g = val)

ALLoutput.N.DF <- ldply(output.N.list, function(x) {
    as.data.frame(x) %>%
        rename(S1 = V1, S2 = V2, S3 = V3) %>%
        mutate(rep = row.names(.)) %>%
        gather(stage, val, -rep) 
}) %>%
    rename(spp = `.id`, N = val)

## N
#ALLoutput.N.DF <- ldply(sppnames, function(x) make.dat(x, biomN = 'N')) %>%
#    rename(N = val)

## Graph biomass
ggplot(ALLoutput.biom.DF, aes(as.numeric(rep), g, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(stage~spp, scales = "free")

## Graph abundance
ggplot(ALLoutput.N.DF, aes(as.numeric(rep), N, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(stage~spp, scales = "free")

## Graph flows
## ggplot(flowdata, aes(Year, SpFloodMag)) +
##     geom_point() + geom_path()
## ggplot(flowdata, aes(Year, BaseDur)) +
##     geom_point() + geom_path()

=====
redo
=====
flowresults.l <- flowresults %>%
    mutate(rep = as.numeric(row.names(.))) %>%  ## may want to change 'rep' here to years and to the df creation
    gather(metric, value, -rep)

ggplot(flowresults.l, aes(rep, value)) +
    geom_point() + geom_path() +
    facet_wrap(~metric)

## Graph all species together
ggplot(ALLoutput.biom.DF, aes(as.numeric(rep), g, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(~spp)

ggplot(ALLoutput.N.DF, aes(as.numeric(rep), N, colour = stage)) + # [ALLoutput.N.DF$spp != "CYLU",]
    geom_point() +
    geom_path() +
    facet_grid(~spp)


##---------------------------------------------------------------------------------------------------------
## Plot summary from all iterations of model run and compare to relative abundance from observed surveys
##---------------------------------------------------------------------------------------------------------
## Reading in observed field data
Verde <- read.csv("data/Rel_Abu_Verde_94-08.csv", header = T) 
head(Verde)

observed <- Verde[, c(1,2,4,5)]
head(observed)
names(observed) <- c('year', 'species', 'obs.mean.rel.abund', 'obs.se.rel.abund')
str(observed)
observed$year <- as.numeric(as.character(observed$year))

## Now have a 'replist' so no need to make the list
## replist <- list(AMNA = replist$AMNA,
##                 CACL = replist$CACL,
##                 CAIN = replist$CAIN,
##                 CYLU = replist$CYLU,
##                 GIRO = replist$GIRO,
##                 LECY = replist$LECY,
##                 MIDO = replist$MIDO)

repdf <- ldply(replist, function(x) {
    adply(x, c(1,2,3))
})
head(repdf)
str(repdf)

names(repdf) <- c('species', 'year', 'stage', 'rep', 'abund')
repdf <- filter(repdf, stage != 'S1')
repdf$year <- as.numeric(as.character(repdf$year))
str(repdf)
head(repdf)
tail(repdf)

str(Total.N)
totn <- adply(Total.N, c(1,2))
names(totn) <- c('year', 'rep', 'tot.abund')
totn$year <- as.numeric(as.character(totn$year))
str(totn)
head(totn)

repdf <- left_join(totn, repdf)
str(repdf)
head(repdf)

repdf <- mutate(repdf, rel.abund = abund/tot.abund)


means <- repdf %>%
    select(-tot.abund) %>%
    group_by(year, rep, species) %>%
    summarise(abund = sum(abund),
              rel.abund = sum(rel.abund)) %>%
    group_by(species, year) %>%
    summarise(mean.abund = mean(abund),
              sd.abund = sd(abund),
              se.abund = sd(abund)/sqrt(iterations),
              mean.rel.abund = mean(rel.abund),
              sd.rel.abund = sd(rel.abund),
              se.rel.abund = sd(rel.abund)/sqrt(iterations))
means

mean_end <- filter(means, year >= 1994)

mean_end <- left_join(mean_end, observed)

theme_classic_facet <- function() {
    theme_classic() +
        theme(strip.background = element_rect(colour = NA, fill = NA))
}

rel.abund.trends <- ggplot(mean_end, aes(year, mean.rel.abund, colour = species, fill = species)) +
    geom_ribbon(aes(ymin = mean.rel.abund - 2*se.rel.abund, ymax = mean.rel.abund + 2*se.rel.abund), colour = 'transparent', alpha = .5, show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    facet_wrap(~species, ncol = 2) +
    theme_classic_facet() +
    coord_cartesian(ylim = c(0,1)) +
    ylab('Relative abundance') +
    xlab('Year')

rel.abund.trends +
    geom_pointrange(aes(y = obs.mean.rel.abund, ymin = obs.mean.rel.abund - 2*obs.se.rel.abund, ymax = obs.mean.rel.abund + 2*obs.se.rel.abund), size = .1, show.legend = FALSE)

ggsave('export/multi-spp2.pdf', width = 4, height = 6)



                                         #par(mfrow=c(4,2), mar = c(4,3,2,1)+ 0.1, oma = c(0,0,0,0))
par(mfrow=c(1,1), mar = c(5,4,4,2), oma=c(0,0,0,0))
                                         # AMNA
AMNA_mean_run <- apply(replist$AMNA[,2:3,], 1, mean)
AMNA_sd_run <- apply(replist$AMNA[,2:3,], 1, sd)
AMNA_SE_run <- AMNA_sd_run/sqrt(iterations)

AMNA_RA <- apply(replist$AMNA[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
AMNA_RelAbu <- apply(AMNA_RA, 1, mean)
AMNA_RA_sd <- apply(AMNA_RA, 1, sd)
AMNA_RA_SE <- AMNA_RA_sd/sqrt(iterations)
AMNA_RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "AMNA"] - AMNA_RelAbu[31:45])^2))

plot(years[31:54], AMNA_RelAbu[31:54], type = "l", main = "AMNA", ylim = c(0, 0.18), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((AMNA_RelAbu[31:54] - 2*AMNA_RA_SE[31:54]), rev(AMNA_RelAbu[31:54] + 2*AMNA_RA_SE[31:54])), col = "#fee0b6", border = F)
lines(years, AMNA_RelAbu, lwd = 2, col = "#f1a340")
points(Verde$Year[Verde$SppCode == "AMNA"], Verde$MeanRelAbu[Verde$SppCode == "AMNA"], pch = 19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "AMNA"] - 2*Verde$SERelAbu[Verde$SppCode == "AMNA"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "AMNA"] + 2*Verde$SERelAbu[Verde$SppCode == "AMNA"], code = 3, length = 0.1, angle = 90, lwd = 2)

                                         #CACL
CACL_mean_run <- apply(replist$CACL[,2:3,], 1, mean)
CACL_sd_run <- apply(replist$CACL[,2:3,], 1, sd)
CACL_SE_run <- CACL_sd_run/sqrt(iterations)


CACL_RA <- apply(replist$CACL[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
CACL_RelAbu <- apply(CACL_RA, 1, mean)
CACL_RA_sd <- apply(CACL_RA, 1, sd)
CACL_RA_SE <- CACL_RA_sd/sqrt(iterations)
CACL_RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CACL"] - CACL_RelAbu[31:45])^2)) # RMSE: sqrt(mean((y-y_pred)^2))

plot(years[31:54], CACL_RelAbu[31:54], type = "l", main = "CACL", ylim = c(0, 0.5), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((CACL_RelAbu[31:54] - 2*CACL_RA_SE[31:54]), rev(CACL_RelAbu[31:54] + 2*CACL_RA_SE[31:54])), col = "#d8daeb", border = F)
lines(years, CACL_RelAbu, lwd = 2, col = "#542788")
points(Verde$Year[Verde$SppCode == "CACL"], Verde$MeanRelAbu[Verde$SppCode == "CACL"], pch = 19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CACL"] - 2*Verde$SERelAbu[Verde$SppCode == "CACL"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CACL"] + 2*Verde$SERelAbu[Verde$SppCode == "CACL"], code = 3, length = 0.1, angle = 90, lwd = 2)

                                         # CAIN
CAIN_mean_run <- apply(replist$CAIN[,2:3,], 1, mean)
CAIN_sd_run <- apply(replist$CAIN[,2:3,], 1, sd)
CAIN_SE_run <- CAIN_sd_run/sqrt(iterations)

CAIN_RA <- apply(replist$CAIN[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
CAIN_RelAbu <- apply(CAIN_RA, 1, mean)
CAIN_RA_sd <- apply(CAIN_RA, 1, sd)
CAIN_RA_SE <- CAIN_RA_sd/sqrt(iterations)

plot(years[31:54], CAIN_RelAbu[31:54], type = "l", main = "CAIN", ylim = c(0, 0.5), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((CAIN_RelAbu[31:54] - 2*CAIN_RA_SE[31:54]), rev(CAIN_RelAbu[31:54] + 2*CAIN_RA_SE[31:54])), col = "#d8daeb", border = F)
lines(years, CAIN_RelAbu, lwd = 2, col = "#542788")
points(Verde$Year[Verde$SppCode == "CAIN"], Verde$MeanRelAbu[Verde$SppCode == "CAIN"], pch=19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CAIN"] - 2*Verde$SERelAbu[Verde$SppCode == "CAIN"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CAIN"] + 2*Verde$SERelAbu[Verde$SppCode == "CAIN"], code = 3, length = 0.1, angle = 90, lwd = 2)


                                         # CYLU
CYLU_mean_run <- apply(replist$CYLU[,2:3,], 1, mean)
CYLU_sd_run <- apply(replist$CYLU[,2:3,], 1, sd)
CYLU_SE_run <- CYLU_sd_run/sqrt(iterations)

CYLU_RA <- apply(replist$CYLU[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
CYLU_RelAbu <- apply(CYLU_RA, 1, mean)
CYLU_RA_sd <- apply(CYLU_RA, 1, sd)
CYLU_RA_SE <- CYLU_RA_sd/sqrt(iterations)

plot(years[31:54], CYLU_RelAbu[31:54], type = "l", main = "CYLU", ylim = c(0, 0.8), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((CYLU_RelAbu[31:54] - 2*CYLU_RA_SE[31:54]), rev(CYLU_RelAbu[31:54] + 2*CYLU_RA_SE[31:54])), col = "#fee0b6", border = F)
lines(years, CYLU_RelAbu, lwd = 2, col = "#f1a340")
points(Verde$Year[Verde$SppCode == "CYLU"], Verde$MeanRelAbu[Verde$SppCode == "CYLU"], pch=19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CYLU"] - 2*Verde$SERelAbu[Verde$SppCode == "CYLU"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CYLU"] + 2*Verde$SERelAbu[Verde$SppCode == "CYLU"], code = 3, length = 0.1, angle = 90, lwd = 2)


                                         # GIRO
GIRO_mean_run <- apply(replist$GIRO[,2:3,], 1, mean)
GIRO_sd_run <- apply(replist$GIRO[,2:3,], 1, sd)
GIRO_SE_run <- GIRO_sd_run/sqrt(iterations)

GIRO_RA <- apply(replist$GIRO[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
GIRO_RelAbu <- apply(GIRO_RA, 1, mean)
GIRO_RA_sd <- apply(GIRO_RA, 1, sd)
GIRO_RA_SE <- GIRO_RA_sd/sqrt(iterations)

plot(years[31:54], GIRO_RelAbu[31:54], type = "l", main = "GIRO", ylim = c(0, 0.3), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((GIRO_RelAbu[31:54] - 2*GIRO_RA_SE[31:54]), rev(GIRO_RelAbu[31:54] + 2*GIRO_RA_SE[31:54])), col = "#d8daeb", border = F)
lines(years[31:54], GIRO_RelAbu[31:54], lwd = 2, col = "#542788")
points(Verde$Year[Verde$SppCode == "GIRO"], Verde$MeanRelAbu[Verde$SppCode == "GIRO"], pch=19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "GIRO"] - 2*Verde$SERelAbu[Verde$SppCode == "GIRO"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "GIRO"] + 2*Verde$SERelAbu[Verde$SppCode == "GIRO"], code = 3, length = 0.1, angle = 90, lwd = 2)


                                         # LECY
LECY_stage1 <- apply(replist$LECY[,1,], 1, mean)
LECY_stage2 <- apply(replist$LECY[,2,], 1, mean)
LECY_stage3 <- apply(replist$LECY[,3,], 1, mean)

plot(years, LECY_stage3, ylim = c(0, 70), xlab = NA, ylab = NA,
     type = "o", pch = 19, lwd = 2, cex = 1.5, col = "blue")
lines(years, LECY_stage1, type = "o", pch = 19, lwd = 2, cex = 1.5, col = "red")
lines(years, LECY_stage2, type = "o", pch = 19, lwd =2, cex = 1.5, col = "green3")
                                         #legend(locator(1), legend = c("S1", "S2", "S3"), pch = 19, col = c("red", "green3", "blue"), xpd = NA)

LECY_mean_run <- apply(replist$LECY[,2:3,], 1, mean)
LECY_sd_run <- apply(replist$LECY[,2:3,], 1, sd)
LECY_SE_run <- LECY_sd_run/sqrt(iterations)

LECY_RA <- apply(replist$LECY[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
LECY_RelAbu <- apply(LECY_RA, 1, mean)
LECY_RA_sd <- apply(LECY_RA, 1, sd)
LECY_RA_SE <- LECY_RA_sd/sqrt(iterations)

plot(years[31:54], LECY_RelAbu[31:54], type = "l", main = "LECY", ylim = c(0, 0.30), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((LECY_RelAbu[31:54] - 2*LECY_RA_SE[31:54]), rev(LECY_RelAbu[31:54] + 2*LECY_RA_SE[31:54])), col = "#fee0b6", border = F)
lines(years[31:54], LECY_RelAbu[31:54], lwd = 2, col = "#f1a340")
points(Verde$Year[Verde$SppCode == "LECY"], Verde$MeanRelAbu[Verde$SppCode == "LECY"], pch=19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "LECY"] - 2*Verde$SERelAbu[Verde$SppCode == "LECY"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "LECY"] + 2*Verde$SERelAbu[Verde$SppCode == "LECY"], code = 3, length = 0.1, angle = 90, lwd = 2)


                                         # MIDO
MIDO_mean_run <- apply(replist$MIDO[,2:3,], 1, mean)
MIDO_sd_run <- apply(replist$MIDO[,2:3,], 1, sd)
MIDO_SE_run <- MIDO_sd_run/sqrt(iterations)

MIDO_RA <- apply(replist$MIDO[,2:3,], c(1,3), sum)/Total.N # check: Total.N[1,,10] = 594.689
MIDO_RelAbu <- apply(MIDO_RA, 1, mean)
MIDO_RA_sd <- apply(MIDO_RA, 1, sd)
MIDO_RA_SE <- MIDO_RA_sd/sqrt(iterations)

plot(years[31:54], MIDO_RelAbu[31:54], type = "l", main = "MIDO", ylim = c(0, 0.5), cex.axis = 2, xlab = NA, ylab = NA)
polygon(c(years[31:54], rev(years[31:54])), c((MIDO_RelAbu[31:54] - 2*MIDO_RA_SE[31:54]), rev(MIDO_RelAbu[31:54] + 2*MIDO_RA_SE[31:54])), col = "#fee0b6", border = F)
lines(years[31:54], MIDO_RelAbu[31:54], lwd=2, col = "#f1a340")
points(Verde$Year[Verde$SppCode == "MIDO"], Verde$MeanRelAbu[Verde$SppCode == "MIDO"], pch=19, cex = 1.5)
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "MIDO"] - 2*Verde$SERelAbu[Verde$SppCode == "MIDO"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "MIDO"] + 2*Verde$SERelAbu[Verde$SppCode == "MIDO"], code = 3, length = 0.1, angle = 90, lwd = 2)

NN_Abu <- apply(cbind(AMNA_RelAbu[31:54], CYLU_RelAbu[31:54], LECY_RelAbu[31:54], MIDO_RelAbu[31:54]),1,sum)
N_Abu <- apply(cbind(CACL_RelAbu[31:54], CAIN_RelAbu[31:54], GIRO_RelAbu[31:54]),1,sum)
SE_NN <- apply(cbind(AMNA_RelAbu[31:54], CYLU_RelAbu[31:54], LECY_RelAbu[31:54], MIDO_RelAbu[31:54]),1, sd)/sqrt(4)
SE_N <- apply(cbind(CACL_RelAbu[31:54], CAIN_RelAbu[31:54], GIRO_RelAbu[31:54]),1, sd)/sqrt(4)
SE_dif <- sd(N_Abu - NN_Abu)/length(N_Abu)

N_Abu/NN_Abu

plot(years[31:54], N_Abu/NN_Abu, type = "l", ylim = c(0, 3), cex.axis = 2, lwd = 3, xlab = n, ylab = n)
abline(h = 1, lty = 2)

plot(years[31:54], N_Abu-NN_Abu, type = "l", ylim = c(-1, 1), cex.axis = 2, lwd = 3, xlab = n, ylab = n)
polygon(c(years[31:54], rev(years[31:54])), 
        c((N_Abu-NN_Abu - 2*SE_dif), rev(N_Abu-NN_Abu + 2*SE_dif)), col = "grey75", border = F)
lines(years[31:54], N_Abu-NN_Abu, lwd=2)
abline(h = 0, lty = 2)


str(replist$CACL)
str(CACL_RA)
write.csv

                                         # Correlation tests
Model.RelAbu <- rbind(AMNA_RelAbu, CACL_RelAbu, CAIN_RelAbu, CYLU_RelAbu, GIRO_RelAbu, LECY_RelAbu, MIDO_RelAbu)


spearman.cor <- c(
    cor(filter(Verde, Year == 1994)$MeanRelAbu, Model.RelAbu[,"1994"], method = "spearman"),
    cor(filter(Verde, Year == 1995)$MeanRelAbu, Model.RelAbu[,"1995"], method = "spearman"),
    cor(filter(Verde, Year == 1996)$MeanRelAbu, Model.RelAbu[,"1996"], method = "spearman"),
    cor(filter(Verde, Year == 1997)$MeanRelAbu, Model.RelAbu[,"1997"], method = "spearman"),
    cor(filter(Verde, Year == 1998)$MeanRelAbu, Model.RelAbu[,"1998"], method = "spearman"),
    cor(filter(Verde, Year == 1999)$MeanRelAbu, Model.RelAbu[,"1999"], method = "spearman"),
    cor(filter(Verde, Year == 2000)$MeanRelAbu, Model.RelAbu[,"2000"], method = "spearman"),
    cor(filter(Verde, Year == 2001)$MeanRelAbu, Model.RelAbu[,"2001"], method = "spearman"),
    cor(filter(Verde, Year == 2002)$MeanRelAbu, Model.RelAbu[,"2002"], method = "spearman"),
    cor(filter(Verde, Year == 2003)$MeanRelAbu, Model.RelAbu[,"2003"], method = "spearman"),
    cor(filter(Verde, Year == 2004)$MeanRelAbu, Model.RelAbu[,"2004"], method = "spearman"),
    cor(filter(Verde, Year == 2005)$MeanRelAbu, Model.RelAbu[,"2005"], method = "spearman"),
    cor(filter(Verde, Year == 2006)$MeanRelAbu, Model.RelAbu[,"2006"], method = "spearman"),
    cor(filter(Verde, Year == 2007)$MeanRelAbu, Model.RelAbu[,"2007"], method = "spearman"),
    cor(filter(Verde, Year == 2008)$MeanRelAbu, Model.RelAbu[,"2008"], method = "spearman"))
spear <- cbind(Verde$Year[1:15], spearman.cor)
round(spear, digits = 2)

spearman.cor.test <- c(
    cor.test(filter(Verde, Year == 1994)$MeanRelAbu, Model.RelAbu[,"1994"], method = "spearman")$p.value, #p-value = 0.0482, rho = 0.786
    cor.test(filter(Verde, Year == 1995)$MeanRelAbu, Model.RelAbu[,"1995"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 1996)$MeanRelAbu, Model.RelAbu[,"1996"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 1997)$MeanRelAbu, Model.RelAbu[,"1997"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 1998)$MeanRelAbu, Model.RelAbu[,"1998"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 1999)$MeanRelAbu, Model.RelAbu[,"1999"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2000)$MeanRelAbu, Model.RelAbu[,"2000"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2001)$MeanRelAbu, Model.RelAbu[,"2001"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2002)$MeanRelAbu, Model.RelAbu[,"2002"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2003)$MeanRelAbu, Model.RelAbu[,"2003"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2004)$MeanRelAbu, Model.RelAbu[,"2004"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2005)$MeanRelAbu, Model.RelAbu[,"2005"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2006)$MeanRelAbu, Model.RelAbu[,"2006"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2007)$MeanRelAbu, Model.RelAbu[,"2007"], method = "spearman")$p.value,
    cor.test(filter(Verde, Year == 2008)$MeanRelAbu, Model.RelAbu[,"2008"], method = "spearman")$p.value)

spear.test <-  cbind(Verde$Year[1:15],spearman.cor.test)
round(spear.test, digits = 3)

round(cbind(Verde$Year[1:15], spearman.cor, spearman.cor.test), 3)

                                         # Species level correlations
head(Verde)
Verde.spp.meanRelAbu <- aggregate(Verde$TotRelAbu, by = list(Verde$SppCode), FUN = mean)
model.spp.meanRelAbu <- rbind(mean(AMNA_RA), mean(CACL_RA), mean(CAIN_RA), mean(CYLU_RA), mean(GIRO_RA), mean(LECY_RA), mean(MIDO_RA))
cor(Verde.spp.meanRelAbu$x, model.spp.meanRelAbu, method = "spearman")
cor.test(Verde.spp.meanRelAbu$x, model.spp.meanRelAbu, method = "spearman")

AMNA_RelAbu[31:45]
Verde.tot <- Verde[,-c(4,5)]
Verde.t <- spread(Verde.tot, SppCode, TotRelAbu) 
spp.cor <- rbind(
    cor(AMNA_RelAbu[31:45], Verde.t$AMNA,  method = "pearson"),
    cor(CACL_RelAbu[31:45], Verde.t$CACL,  method = "pearson"),
    cor(CAIN_RelAbu[31:45], Verde.t$CAIN,  method = "pearson"),
    cor(CYLU_RelAbu[31:45], Verde.t$CYLU,  method = "pearson"),
    cor(GIRO_RelAbu[31:45], Verde.t$GIRO,  method = "pearson"),
    cor(LECY_RelAbu[31:45], Verde.t$LECY,  method = "pearson"),
    cor(MIDO_RelAbu[31:45], Verde.t$MIDO,  method = "pearson"))
rownames(spp.cor) <- c("AMNA", "CACL", "CAIN", "CYLU", "GIRO", "LECY", "MIDO")
round(spp.cor, digits = 2)

spp.cor.test <- rbind(
    cor.test(AMNA_RelAbu[31:45], Verde.t$AMNA,  method = "pearson")$p.value,
    cor.test(CACL_RelAbu[31:45], Verde.t$CACL,  method = "pearson")$p.value,
    cor.test(CAIN_RelAbu[31:45], Verde.t$CAIN,  method = "pearson")$p.value,
    cor.test(CYLU_RelAbu[31:45], Verde.t$CYLU,  method = "pearson")$p.value,
    cor.test(GIRO_RelAbu[31:45], Verde.t$GIRO,  method = "pearson")$p.value,
    cor.test(LECY_RelAbu[31:45], Verde.t$LECY,  method = "pearson")$p.value,
    cor.test(MIDO_RelAbu[31:45], Verde.t$MIDO,  method = "pearson")$p.value)
rownames(spp.cor.test) <- c("AMNA", "CACL", "CAIN", "CYLU", "GIRO", "LECY", "MIDO")
round(spp.cor.test, digits = 2)
round(cbind(spp.cor, spp.cor.test[,1]), digits = 2)

                                         # K Sensitivity analysis--------------------------------------------------------------------------------------------
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K28596.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K33362.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K38128.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K42894.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K47660.csv") #mean(Model.RelAbu["AMNA_RelAbu",])
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K52426.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K57192.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K61958.csv")
                                         # write.csv(Model.RelAbu, file = "output/K_sensitivity/RelAbu_K66724.csv")

RelAbu_47660 <- read.csv("output/K_sensitivity/RelAbu_K47660.csv", row.names = 1)
RelAbu_52426 <- read.csv("output/K_sensitivity/RelAbu_K52426.csv", row.names = 1)
RelAbu_57192 <- read.csv("output/K_sensitivity/RelAbu_K57192.csv", row.names = 1)
RelAbu_61958 <- read.csv("output/K_sensitivity/RelAbu_K61958.csv", row.names = 1)
RelAbu_66724 <- read.csv("output/K_sensitivity/RelAbu_K66724.csv", row.names = 1)
RelAbu_42894 <- read.csv("output/K_sensitivity/RelAbu_K42894.csv", row.names = 1)
RelAbu_38128 <- read.csv("output/K_sensitivity/RelAbu_K38128.csv", row.names = 1)
RelAbu_33362 <- read.csv("output/K_sensitivity/RelAbu_K33362.csv", row.names = 1)
RelAbu_28596 <- read.csv("output/K_sensitivity/RelAbu_K28596.csv", row.names = 1)

perdif_plus10 <- ((RelAbu_52426 - RelAbu_47660)/RelAbu_47660)*100
perdif_plus20 <- ((RelAbu_57192 - RelAbu_47660)/RelAbu_47660)*100
perdif_plus30 <- ((RelAbu_61958 - RelAbu_47660)/RelAbu_47660)*100
perdif_plus40 <- ((RelAbu_66724 - RelAbu_47660)/RelAbu_47660)*100
perdif_null <- ((RelAbu_47660 - RelAbu_47660)/RelAbu_47660)*100
perdif_minus10 <- ((RelAbu_42894 - RelAbu_47660)/RelAbu_47660)*100
perdif_minus20 <- ((RelAbu_38128 - RelAbu_47660)/RelAbu_47660)*100
perdif_minus30 <- ((RelAbu_33362 - RelAbu_47660)/RelAbu_47660)*100
perdif_minus40 <- ((RelAbu_28596 - RelAbu_47660)/RelAbu_47660)*100

                                         #AMNA.perdif_plus10 <- mean(as.matrix(perdif_plus10["AMNA_RelAbu",]), na.rm = T)

spp_mean_minus40 <- apply(as.matrix(perdif_minus40), 1, mean, na.rm = T)
spp_mean_minus30 <- apply(as.matrix(perdif_minus30), 1, mean, na.rm = T)
spp_mean_minus20 <- apply(as.matrix(perdif_minus20), 1, mean, na.rm = T)
spp_mean_minus10 <- apply(as.matrix(perdif_minus10), 1, mean, na.rm = T)
spp_mean_null <- apply(as.matrix(perdif_null), 1, mean, na.rm =T)
spp_mean_plus10 <- apply(as.matrix(perdif_plus10), 1, mean, na.rm = T)
spp_mean_plus20 <- apply(as.matrix(perdif_plus20), 1, mean, na.rm = T)
spp_mean_plus30 <- apply(as.matrix(perdif_plus30), 1, mean, na.rm = T)
spp_mean_plus40 <- apply(as.matrix(perdif_plus40), 1, mean, na.rm = T)

perdiff <- rbind(spp_mean_minus40, spp_mean_minus30, spp_mean_minus20, spp_mean_minus10, spp_mean_null,
                 spp_mean_plus10, spp_mean_plus20, spp_mean_plus30, spp_mean_plus40)

                                         #plot(1:2, perdiff[,"AMNA_RelAbu"])

par(mfrow = c(1,1), mar=c(5,4,4,2))


plot(seq(-40, 40, 10),
     perdiff[,"AMNA_RelAbu"],
     ylim = c(-60, 70), bty = "l",
     xlab = "% change in K from 47660", ylab = "% change in RelAbu",
     pch = 21, col = "#b35806", bg = "#b35806", cex = 1.5)
points(seq(-40, 40, 10), perdiff[,"CACL_RelAbu"],
       pch = 21, col = "#f1a340", bg = "#f1a340", cex = 1.5)
points(seq(-40, 40, 10), perdiff[,"CAIN_RelAbu"],
       pch = 21, col = "#fee0b6", bg = "#fee0b6", cex = 1.5)
points(seq(-40, 40, 10), perdiff[,"CYLU_RelAbu"],
       pch = 21, col = "#252525", bg = "#252525", cex = 1.5)
points(seq(-40, 40, 10), perdiff[,"GIRO_RelAbu"],
       pch = 21, col = "#d8daeb", bg = "#d8daeb", cex = 1.5)
points(seq(-40, 40, 10), perdiff[,"LECY_RelAbu"],
       pch = 21, col = "#998ec3", bg = "#998ec3", cex = 1.5)
points(seq(-40, 40, 10), perdiff[,"MIDO_RelAbu"],
       pch = 21, col = "#542788", bg = "#542788", cex = 1.5)
abline(h=5, lty=2)
abline(h=-5, lty=2)
legend(locator(1), legend = c("AMNA","CACL","CAIN", "CYLU", "GIRO", "LECY", "MIDO"), pch = c(rep(21,7)), xpd = NA,
       col = c("#b35806", "#f1a340", "#fee0b6", "#252525", "#d8daeb", "#998ec3", "#542788"),
       pt.bg = c("#b35806", "#f1a340", "#fee0b6", "#252525", "#d8daeb", "#998ec3", "#542788"))
