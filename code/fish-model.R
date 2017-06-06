
# Required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(popbio)

# SETUP ----------------------------------------------------------------------------------

rm(list = ls()) # clearing the workspace 

# bringing in flow data
# Verde flow data at Paulden 7/17/1963-2017, 54 years continuous
flowdata <- read.csv("data/flowdata_Verde.csv") 
str(flowdata)
head(flowdata)
# FloodMag - magnitude of flood in cfs
# BaseDur - baseflow duration in days
# flooddate is peak dates of all floods (Oct 1 = 1)

count <- length(flowdata$SpFloodMag) # inner loop - number of years to simulate; starting with natural sequence of flow record
# burnin = 100 # number of years to discard as burn in during long term mean estimation 
# outerreps <- 1 # number of iterations for outer loop that alters drought/flood frequency 
# replicates <- 100 # number of replicate projections to run (mid loop)

# Modifiers
modifiers <- read.csv('data/modifiers.csv')

# adding 'Modifier' value from csv to 'Code' in csv
for(j in 1:length(modifiers[,1])) {
    nam <- paste(modifiers[j,4])
    assign(nam, modifiers[j,6]) # SHOULD BE [j,5] FOR REAL VALUES, [j,6] for null values!!!!
}

# Vital rates
# Baseline maturation probabiliity, aCACL3 (adult senescence rate)
# Background mortality
# Initial volume in grams in 100-m reach
# Fecundity based on year type
# these are raw numbers taken from survey data: max, mean and min fecundities
# Stage specific densities (ind./g)

vitalrates <- read.csv('data/vital-rates.csv')

# assigning vital rate values from column 3 to 'code' in column 2
for(k in 1:length(vitalrates[,1])) {
    nam <- paste(vitalrates[k,2])
    assign(nam, vitalrates[k,3]) 
}

# Key ------------------------------------------------------------
# CACL: desert sucker
# GIRO: chub
# LECY: green sunfish

# Average total volume of water per reach in m3: 307
# Average total fish biomass per reach in g: 4766


K = 4766 # avg. for 100-m reach across 6 replicate reaches... in g/m3 this is 15.5 g/m3

# Bunch of functions -----------------------------------------

### :CLARIFY: This rescue function needs to be sorted out. Is it realistic?
# Random rescue function - calculates whether rescue occurs if less than 2 adults are present (or pop extinction occurs)
rescue <- rbinom(1, 1, .5)

# Checks to see if at least two adults are present to allow spawning
# If there are, it carries on, if not, there's a potential for population rescue, via the above random no. generator. i.e. if rescued, 2 individuals are present, if not, 0. 
adult_func <- function(x) { 
    ifelse(x > 1.99999, x, 2)
}


# checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
    ifelse(x < 0, 0, x)
}

# Keeps F..3 from dividing by zero by substituting an arbitrary nonzero number that will get
# multiplied by zero later anyway during matrix multiplication
nonind <- function(x) {
    ifelse(x == 0, 666, x)
}

# timing function for native fish spawning. floods have to occur within a set window for spawning to occur

spawnwindow.CACL_func <- function(x) {
    ifelse(flowdata$SpFloodDate >= 107 & flowdata$SpFloodDate <= 213, 1, 0) # flood can occur within 15 days of starting spawn (water day 124)
}

spawnwindow.GIRO_func <- function(x) {
    ifelse(flowdata$SpFloodDate >= 138 & flowdata$SpFloodDate <= 244, 1, 0) # flood can occur within 15 days of starting spawn (water day 153)
}

wipeoutwindow.LECY_func <- function(x) {
    ifelse(flowdata$SuFloodDate >= 214 & flowdata$SuFloodDate <= 305, 1, 0) # late floods stop spawning
}


spawnwindow.CACL <-
    spawnwindow.CACL_func(flowdata$SpFloodDate) # Spring flood date

spawnwindow.GIRO <-
    spawnwindow.GIRO_func(flowdata$SpFloodDate) # Spring flood date

wipeoutwindow.LECY <-
    wipeoutwindow.LECY_func(flowdata$SuFloodDate) # Summer flood date


# FLOOD THRESHOLD FUNCTIONS --------------------------------------------------------------
# flowdata$FloodMag - vector containing peak flood magnitude 

# Magnitude of peak flow over which is considered a large flood event
SP_highfloodcutoff = 700 # SP: Jan 1 - Apr 30 (water_day 93-213). Floods (cfs) at 4-yr recurrence interval or greater (30 times median flow).
SU_highfloodcutoff = 200 # SU: May 1 - Sep 30 (Water_day 214-365). Floods (Cfs) at 4-year recurrence interval for summer
medfloodcutoff = 220 # SP Floods (cfs) at 2.5-yr recurrence interval (10 times median flow)
# Non-event is same for drought (threshold below rather than above)
mindroughtlength = 40 # threshold length of low-flow in days (greater than 75th percentile of low flow duration) -- 
# determined as number of consecutive days with discharge (cfs) less than 22 cfs (25th percentile of flows) May 1 - Sep 30
# This threshold is compared to baseflow duration in any given year - a vector specified in csv (baseflow_dur).

# Convert peak discharge values into a vector of floods/no floods, 
SP_highflood_func <- function(x) {
    ifelse(x >= SP_highfloodcutoff, 1, 0)
}

SU_highflood_func <- function(x) {
  ifelse(x >= SU_highfloodcutoff, 1, 0)
}

medflood_func <- function(x) {
    ifelse(x >= medfloodcutoff & x <= SP_highfloodcutoff, 1, 0)
}

nonevent_func <- function(Spfl, BD, Sufl) {
    ifelse(Spfl < medfloodcutoff & BD < mindroughtlength & Sufl < SU_highfloodcutoff, 1, 0)
}

drought_func <- function(Spfl, BD, Sufl) {
    ifelse(Spfl < medfloodcutoff & BD >= mindroughtlength & Sufl < SU_highfloodcutoff, 1, 0)
}


# Defining year types here for now. May need to put into outer loop if simulating
SP_highflood <- SP_highflood_func(flowdata$SpFloodMag) 
#sum(SP_highflood == 1)
SU_highflood <- SU_highflood_func(flowdata$SuFloodMag)
#sum(SU_highflood ==1)
medflood <- medflood_func(flowdata$SpFloodMag)
#sum(medflood ==1)
drought <- drought_func(Spfl = flowdata$SpFloodMag, BD = flowdata$BaseDur, Sufl = flowdata$SuFloodMag)
#sum(drought ==1)
nonevent <- nonevent_func(Spfl = flowdata$SpFloodMag, BD = flowdata$BaseDur, Sufl = flowdata$SuFloodMag) 
#sum(nonevent ==1)
flood <- SP_highflood + medflood # simply whether it's a flood year or not
  ## Four years with SP_highflood and SU_highflood

# ITERATION PARAMETERS -------------------------------------------------------------------
# Setting up arrays/vectors to fill with data from loops

# Inner loop details ---------------------------------------------------------------------
# 'count' - number of years to project simulations (inner loop)

# Output of biomass and no. ind. for each age class for each year projected  
# An array with 3 columns (each stage class) and however many rows there are years projected 
CACLoutput.N <- array(0, dim = c(count, 3))
CACLoutput.biom <- array(0, dim = c(count, 3))
GIROoutput.N <- array(0, dim = c(count, 3))
GIROoutput.biom <- array(0, dim = c(count, 3))
LECYoutput.N <- array(0, dim = c(count, 3))
LECYoutput.biom <- array(0, dim = c(count, 3))

# Total biomass as % of K 
CACLbiomoutput <- numeric(length = count)
GIRObiomoutput <- numeric(length = count)
LECYbiomoutput <- numeric(length = count)

# Flood and drought settings for each year projected into the future (i.e. 0 or 1) 
SPhighfloodoutput <- numeric(length = count) 
SUhighfloodoutput <- numeric(length = count)
medfloodoutput <- numeric(length = count) 
droughtoutput <- numeric(length = count) 
noneventoutput <- numeric(length = count) 

# Fecundities
FCACLoutput <- numeric(length = count)
FGIROoutput <- numeric(length = count)
FLECYoutput <- numeric(length = count)


# N gives the total number of individuals for each age class.
# Initially here, this is found by multiplying the number of g occupied by a given class
# by the density per g
# biom = g/m3
# den = indiv/g

#NCACL <- c(biomCACL1 * denCACL1,
#           biomCACL2 * denCACL2,
#           biomCACL3 * denCACL3) 

biomCACL <- c(biomCACL1,
              biomCACL2,
              biomCACL3) 

#NGIRO <- c(biomGIRO1 * denGIRO1,
#           biomGIRO2 * denGIRO2,
#           biomGIRO3 * denGIRO3)

biomGIRO <- c(biomGIRO1,
              biomGIRO2,
              biomGIRO3)

#NLECY <- c(biomLECY1 * denLECY1,
#           biomLECY2 * denLECY2,
#           biomLECY3 * denLECY3) 

biomLECY <- c(biomLECY1,
              biomLECY2,
              biomLECY3) 


# Inner loop #############################################################################
for(i in 1:count) {

    y = sample(nrow(flowdata), 1) 
#y = i 
# y was a random number within the length of the flow data to randomly select a year from 
# the 'flood' and 'drought' vector. now it's directly taken from flow vector. 

# VITAL RATE DEFINITIONS: Desert Sucker  -----------------------------------------------------
# G is prob. of transition to next stage
# P is prob. of remaining in that stage

# flood is a vector of either medium or high magnitude
# Summer mortality vital rate is multiplied by multiplier based on yeartype as specified above
# Different from plant model as we have built in multiplier of severity of droughts or non-events on transition probabilities 

# Desert Sucker - 
# YOY (GCACL1) survival and recruitment depends on spring flows
    GCACL1 <- aCACL1 *
        ((1 - (SP_highflood[y] * S1MortCACL)) * CACL_Sp_HF) *
        ((1 - (medflood[y] * S1MortCACL)) * CACL_Sp_MF) *
        ((1 - (nonevent[y] * S1MortCACL)) * CACL_Sp_NE) *
        ((1 - (drought[y] * S1MortCACL)) * CACL_Sp_DR)
        
    GCACL2 <- aCACL2 *
        ((1 - (SP_highflood[y] * S2MortCACL)) * CACL_Su_HF) *
        ((1 - (medflood[y] * S2MortCACL)) * CACL_Su_MF) *
        ((1 - (nonevent[y] * S2MortCACL)) * CACL_Su_NE) *
        ((1 - (drought[y] * S2MortCACL))  * CACL_Su_DR)

    PCACL3 <- (1 - aCACL3) *
        ((1 - (SP_highflood[y] * S3MortCACL)) * CACL_Su_HF) *
        ((1 - (medflood[y] * S3MortCACL)) * CACL_Su_MF) *
        ((1 - (nonevent[y] * S3MortCACL)) * CACL_Su_NE) *
        ((1 - (drought[y] * S3MortCACL)) *  CACL_Su_DR)
    

# Chub

     GGIRO1 <- aGIRO1 *
        ((1 - (SP_highflood[y] * S1MortGIRO)) * GIRO_Sp_HF) *
        ((1 - (medflood[y] * S1MortGIRO)) * GIRO_Sp_MF) *
        ((1 - (nonevent[y] * S1MortGIRO)) * GIRO_Sp_NE) *
        ((1 - (drought[y] * S1MortGIRO)) * GIRO_Sp_DR)
        
    GGIRO2 <- aGIRO2 *
        ((1 - (SP_highflood[y] * S2MortGIRO)) * GIRO_Su_HF) *
        ((1 - (medflood[y] * S2MortGIRO)) * GIRO_Su_MF) *
        ((1 - (nonevent[y] * S2MortGIRO)) * GIRO_Su_NE) *
        ((1 - (drought[y] * S2MortGIRO)) * GIRO_Su_DR)

    PGIRO3 <- (1 - aGIRO3) *
        ((1 - (SP_highflood[y] * S3MortGIRO)) * GIRO_Su_HF) *
        ((1 - (medflood[y] * S3MortGIRO)) * GIRO_Su_MF) *
        ((1 - (nonevent[y] * S3MortGIRO)) * GIRO_Su_NE) *
        ((1 - (drought[y] * S3MortGIRO)) * GIRO_Su_DR)
    
 


# Green sunfish - note the "NN"

 GLECY1 <- aLECY1 *
        ((1 - (SU_highflood[y] * S1MortLECY)) * LECY_Sp_HF) *
        ((1 - (medflood[y] * S1MortLECY)) * LECY_Sp_MF) *
        ((1 - (nonevent[y] * S1MortLECY)) * LECY_Sp_NE) *
        ((1 - (drought[y] * S1MortLECY)) * LECY_Sp_DR)
        
    GLECY2 <- aLECY2 *
        ((1 - (SU_highflood[y] * S2MortLECY))  * LECY_Su_HF) *
        ((1 - (medflood[y] * S2MortLECY)) * LECY_Su_MF) *
        ((1 - (nonevent[y] * S2MortLECY)) * LECY_Su_NE) *
        ((1 - (drought[y] * S2MortLECY)) * LECY_Su_DR)

    PLECY3 <- (1 - aLECY3) *
        ((1 - (SU_highflood[y] * S3MortLECY)) * LECY_Su_HF) *
        ((1 - (medflood[y] * S3MortLECY)) * LECY_Su_MF) *
        ((1 - (nonevent[y] * S3MortLECY)) * LECY_Su_NE) *
        ((1 - (drought[y] * S3MortLECY)) * LECY_Su_DR)
    
 
# Total grams occupied after year -----------------------------------------------

totbiom.CACL <- 
    biomCACL[1] + #
    biomCACL[2] + 
    biomCACL[3] 

totbiom.GIRO <- 
    biomGIRO[1] + #
    biomGIRO[2] + 
    biomGIRO[3] 

totbiom.LECY <- 
    biomLECY[1] + #
    biomLECY[2] + 
    biomLECY[3] 



### :CLARIFY: at the moment carrying capacity is limiting spawning of all species based on the total biomass occupied at the end of the previous year. i.e. if above K, no spp spawn in that year. If spawning occurs, they can all spawn and there is no sequence, so overseeding COULD be massive.
#### Watch the brackets, especially related to floor(), too early rounding end in zeros.
# POTENTIAL CACL FECUNDITY ---------------------------------------------------------
    FCACL3 <- checkpos(####(adult_func(biomCACL[3] * denCACL3)) * # checks to see if at least 2 adults are present.
                       ####(1/nonind(biomCACL[3])) *
    ifelse(SP_highflood[y] == 1 & spawnwindow.CACL[y] == 1,
    (floor((biomCACL[3] * denCACL3) / 2 * meanfec.CACL / denCACL1)), # converting maxfec into g 
    ifelse(medflood[y] == 1 & spawnwindow.CACL[y] == 1,
           (floor((biomCACL[3] * denCACL3) / 2 * meanfec.CACL / denCACL1)),
           (floor((biomCACL[3] * denCACL3) / 2 * meanfec.CACL / denCACL1)))) * # fecundity function of yeartype and whether flood occurred during spawning window
                   ifelse((totbiom.CACL + totbiom.GIRO + totbiom.LECY) < K, 1, 0)) # if K is already occupied, then no fecundity
                     
# '(1/nonind(NCACL[3]))' keeps FCACL3 from dividing by zero by substituting an arbitrary non-0 
# number that will be multiplied by 0 later anyway during matrix multiplication
# This means it's fecund per individual, rather than having overall fecundity double multipled by overall biomass (i.e. here and in matrix mult)     

?floor
# POTENTIAL GIRO FECUNDITY ---------------------------------------------------------
    FGIRO3 <- checkpos(####(adult_func(biomGIRO[3] * denGIRO3)) *  # checks to see if at least 2 adults are present. 
                            ####(1/nonind(biomGIRO[3])) *
                       ifelse(SP_highflood[y] == 1 & spawnwindow.GIRO[y] == 1,
                       (floor((biomGIRO[3] * denGIRO3) / 2 * meanfec.GIRO / denGIRO1)), # converting maxfec into g 
                       ifelse(medflood[y] == 1 & spawnwindow.GIRO[y] == 1,
                       (floor((biomGIRO[3] * denGIRO3) / 2 * meanfec.GIRO / denGIRO1)),
                       (floor((biomGIRO[3] * denGIRO3) / 2 * meanfec.GIRO / denGIRO1)))) * # fecundity function of yeartype
                       ifelse((totbiom.CACL + totbiom.GIRO + totbiom.LECY) < K, 1, 0)) # if K is already occupied, then no fecundity

# POTENTIAL LECY FECUNDITY ---------------------------------------------------------
    FLECY3 <- checkpos(####(adult_func(biomLECY[3] * denLECY3)) * # checks to see if at least 2 adults are present.
                                  ####(1/nonind(biomLECY[3])) *
                       ifelse(SU_highflood[y] == 1 & wipeoutwindow.LECY[y] == 1, 0, # converting maxfec into g
                       ifelse(SU_highflood[y] == 1 & wipeoutwindow.LECY[y] == 0, (floor((biomLECY[3] * denLECY3) / 2 * meanfec.LECY / denLECY1)), # minfec
                       ifelse(medflood[y] == 1 & wipeoutwindow.LECY[y] == 1, (floor((biomLECY[3] * denLECY3) / 2 * meanfec.LECY / denLECY1)), (floor((biomLECY[3] * denLECY3) / 2 * meanfec.LECY / denLECY1))))) * # fecundity function of yeartype
                       ifelse((totbiom.CACL + totbiom.GIRO + totbiom.LECY) < K, 1, 0)) # if K is already occupied, then no fecundity


# K --------------------------------------------------------------------------------------
# CACL
# gives total CACL biomass (g) as a percentage of K; 
KCACL <- 100 * (biomCACL[1] + 
             biomCACL[2] + 
             biomCACL[3])/K

# GIRO
# gives total GIRO biomass as a percentage of K; 
KGIRO <- 100 * (biomGIRO[1] + 
             biomGIRO[2] + 
             biomGIRO[3])/K

# LECY
# gives total LECY biomass as a percentage of K; 
KLECY <- 100 * (biomLECY[1] + 
             biomLECY[2] + 
             biomLECY[3])/K
    
# TRANSITION MATRICES --------------------------------------------------------------------

# TRANSITION MATRIX FOR CACL -------------------------------------------------------
ACACL1 <- c(0, 0, FCACL3)
ACACL2 <- c(GCACL1, 0, 0)
ACACL3 <- c(0, GCACL2, PCACL3)
# Matrix
ACACL <- rbind(ACACL1, ACACL2, ACACL3)
lambda(ACACL) # Checking population growth rate
# TRANSITION MATRIX FOR GIRO -------------------------------------------------------
AGIRO1 <- c(0, 0, FGIRO3)
AGIRO2 <- c(GGIRO1, 0, 0)
AGIRO3 <- c(0, GGIRO2, PGIRO3)
# Matrix
AGIRO <- rbind(AGIRO1, AGIRO2, AGIRO3)
lambda(AGIRO) # Checking population growth rate

# TRANSITION MATRIX FOR LECY -------------------------------------------------------
ALECY1 <- c(0, 0, FLECY3)
ALECY2 <- c(GLECY1, 0, 0)
ALECY3 <- c(0, GLECY2, PLECY3)
# Matrix
ALECY <- rbind(ALECY1, ALECY2, ALECY3)
lambda(ALECY) # Checking population growth rate

# Manually adding fecundity

## ACACL1 <- c(0, 0, 0)
## ACACL2 <- c(GCACL1, 0, 0)
## ACACL3 <- c(0, GCACL2, PCACL3)
## # Matrix
## ACACL <- rbind(ACACL1, ACACL2, ACACL3)

## # TRANSITION MATRIX FOR GIRO -------------------------------------------------------
## AGIRO1 <- c(0, 0, 0)
## AGIRO2 <- c(GGIRO1, 0, 0)
## AGIRO3 <- c(0, GGIRO2, PGIRO3)
## # Matrix
## AGIRO <- rbind(AGIRO1, AGIRO2, AGIRO3)

## # TRANSITION MATRIX FOR LECY -------------------------------------------------------
## ALECY1 <- c(0, 0, 0)
## ALECY2 <- c(GLECY1, 0, 0)
## ALECY3 <- c(0, GLECY2, PLECY3)
## # Matrix
## ALECY <- rbind(ALECY1, ALECY2, ALECY3)


# COMPILING OUTPUTS ----------------------------------------------------------------------
# CACL
CACLoutput.biom[i,1:3] <- biomCACL # array of biomass of each age class for each yr projected. biomCACL = total biomass for each age class
CACLbiomoutput[i] <- KCACL # total biomass as % of K; this is in g
FCACLoutput[i] <- FCACL3
    
# GIRO
GIROoutput.biom[i,1:3] <- biomGIRO # array of biomass of each age class for each yr projected. biomGIRO = total biomass for each age class
GIRObiomoutput[i] <- KGIRO # total biomass as % of K; this is in g
FGIROoutput[i] <- FGIRO3
    
# LECY
LECYoutput.biom[i,1:3] <- biomLECY # array of biomass of each age class for each yr projected. biomLECY = total biomass for each age class
LECYbiomoutput[i] <- KLECY # total biomass as % of K; this is in g 
FLECYoutput[i] <- FLECY3


# Records flood settings of each particular projected year (0 for nonflood, 1 for flood)
SPhighfloodoutput[i] <- SP_highflood[y]
medfloodoutput[i] <- medflood[y]
noneventoutput[i] <- nonevent[y]
droughtoutput[i] <- drought[y]

# if we want to add total fecundity in afterwards
CACLplaceholder <- FCACL3
GIROplaceholder <- FGIRO3
LECYplaceholder <- FLECY3
    
# MATRIX MULTIPLICATION ------------------------------------------------------------------
# CACL
biomCACL <- ACACL %*% biomCACL # ACACL is transition matrix, biomCACL = total biomass for each age class
#    biomCACL[1] <- CACLplaceholder
    
#nCACLpost <- c(biomCACL[1] * denCACL1, biomCACL[2] * denCACL2, biomCACL[3] * denCACL3)
#nCACLpost <- ifelse(nCACLpost < 1, 0, 1)
#biomCACL <- biomCACL * nCACLpost

# GIRO
biomGIRO <- AGIRO %*% biomGIRO # AGIRO is transition matrix, biomGIRO = total biomass for each age class
#    biomGIRO[1] <- GIROplaceholder
    
#nGIROpost <- c(biomGIRO[1] * denGIRO1, biomGIRO[2] * denGIRO2, biomGIRO[3] * denGIRO3)
#nGIROpost <- ifelse(nGIROpost < 1, 0, 1)
#biomGIRO <- biomGIRO * nGIROpost

# LECY
biomLECY <- ALECY %*% biomLECY # ALECY is transition matrix, biomLECY = total biomass for each age class
#    biomLECY[1] <- LECYplaceholder
    
#nLECYpost <- c(biomLECY[1] * denLECY1, biomLECY[2] * denLECY2, biomLECY[3] * denLECY3)
#nLECYpost <- ifelse(nLECYpost < 1, 0, 1)
#biomLECY <- biomLECY * nLECYpost

} # End of inner loop ####################################################################


CACLoutput.biom.DF <- as.data.frame(CACLoutput.biom) %>%
    rename(S1 = V1, S2 = V2, S3 = V3) %>%
    mutate(rep = row.names(.)) %>%
    gather(stage, g, -rep) %>%
    mutate(spp = 'CACL') 
    
GIROoutput.biom.DF <- as.data.frame(GIROoutput.biom) %>%
    rename(S1 = V1, S2 = V2, S3 = V3) %>%
    mutate(rep = row.names(.)) %>%
    gather(stage, g, -rep) %>%
    mutate(spp = 'GIRO')

LECYoutput.biom.DF <- as.data.frame(LECYoutput.biom) %>%
    rename(S1 = V1, S2 = V2, S3 = V3) %>%
    mutate(rep = row.names(.)) %>%
    gather(stage, g, -rep) %>%
    mutate(spp = 'LECY')
   

ALLoutput.biom.DF <- rbind(CACLoutput.biom.DF, GIROoutput.biom.DF, LECYoutput.biom.DF)

head(ALLoutput.biom.DF)
ALLoutput.biom.DF
ggplot(ALLoutput.biom.DF, aes(as.numeric(rep), g, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(stage~spp, scales = "free")


ggplot(flowdata, aes(Year, SpFloodMag)) +
    geom_point() + geom_path()
ggplot(flowdata, aes(Year, BaseDur)) +
    geom_point() + geom_path()

flowresults <- as.data.frame(cbind(SPhighfloodoutput, medfloodoutput, noneventoutput, droughtoutput)) %>%
    mutate(rep = as.numeric(row.names(.))) %>%
    gather(metric, value, -rep)

ggplot(flowresults, aes(rep, value)) +
    geom_point() + geom_path() +
    facet_wrap(~metric)


ggplot(ALLoutput.biom.DF, aes(as.numeric(rep), g, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(~spp)

                                        #
