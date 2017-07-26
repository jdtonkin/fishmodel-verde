# Verde Fish Model
# Jane Rogosch, Jono Tonkin, et al.
# 01-Mar-17
# 26-Jul-17

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

burnin <- 30 # number of years to discard as burn in, randomly sampled from flow record
count <- burnin + length(flowdata$SpFloodMag) # inner loop - number of years to simulate; starting with natural sequence of flow record
# outerreps <- 1 # number of iterations for outer loop that alters drought/flood frequency 
iterations <- 1000 # number of replicate projections to run (mid loop)

# Modifiers
modifiers <- read.csv('data/modifiers-all-spp.csv')

# adding 'Modifier' value from csv to 'Code' in csv
for(j in 1:length(modifiers[,1])) {
    nam <- paste(modifiers[j,4])
    assign(nam, modifiers[j,5]) # SHOULD BE [j,5] FOR REAL VALUES, [j,6] for null values!!!!
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
# CAIN: sonora sucker
# MIDO: smallmouth bass
# CYLU: red shiner
# AMNA: yellow bullhead

# Average total volume of water per 100 m reach in m3: 307
# Average total fish biomass per 100 m reach in g: 4766
# Average total biomass Bonar 2004 in g/100m2: 606

K = 47660 # avg. for 1-km reach across 6 replicate reaches... in g/m3 this is 155 g/m3

# Bunch of functions -----------------------------------------

### :CLARIFY: This rescue function needs to be sorted out. Is it realistic?
# Random rescue function - calculates whether rescue occurs if less than 2 adults are present (or pop extinction occurs)
rescue <- rbinom(1, 1, .5)

# Checks to see if at least two adults are present to allow spawning
# If there are, it carries on, if not, there's a potential for population rescue, via the above random no. generator. i.e. if rescued, 2 individuals are present, if not, 0. 
adult_func <- function(x) { 
    ifelse(x > 1.99999, x, 2)
}

# # rescue with 2 adults
# adult_res <- function(x) {
#   ifelse(x > 1, x, 2*rbinom(1, 1, 0.5))
# }

# rescue biomass for 2 individuals based on density (indiv/g) of species
biom_res <- function(x, spp) { # x is biom object for a given species (e.g. biomCACL), spp is name of species in quotes (e.g. "CACL")
  s.x <- x[3]
  x[3] <- ifelse(s.x > 2*(1/get(noquote(paste0("den", spp, "3")))), s.x, 2*(1/get(noquote(paste0("den",spp,"3"))))*rbinom(1, 1, 0.25))
  print(x)
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

spawnwindow.CAIN_func <- function(x) {
  ifelse(flowdata$SpFloodDate >= 107 & flowdata$SpFloodDate <= 244, 1, 0) # flood can occur within 15 days of starting spawn (water day 124)
}

wipeoutwindow.LECY_func <- function(x) {
    ifelse(flowdata$SuFloodDate >= 214 & flowdata$SuFloodDate <= 305, 1, 0) # late floods stop spawning
}

wipeoutwindow.MIDO_func <- function(x) {
  ifelse(flowdata$SuFloodDate >= 214 & flowdata$SuFloodDate <= 244, 1, 0) # late floods stop spawning
}

    # Nothing stops CYLU, they spawn April to October so they can take advantage of any low flow window

wipeoutwindow.AMNA_func <- function(x) {
  ifelse(flowdata$SuFloodDate >= 184 & flowdata$SuFloodDate <= 274, 1, 0) # late floods stop spawning
}

spawnwindow.CACL <-
    spawnwindow.CACL_func(flowdata$SpFloodDate) # Spring flood date

spawnwindow.GIRO <-
    spawnwindow.GIRO_func(flowdata$SpFloodDate) # Spring flood date

spawnwindow.CAIN <-
  spawnwindow.CAIN_func(flowdata$SpFloodDate) # Spring flood date

wipeoutwindow.LECY <-
    wipeoutwindow.LECY_func(flowdata$SuFloodDate) # Summer flood date

wipeoutwindow.MIDO <-
  wipeoutwindow.MIDO_func(flowdata$SuFloodDate) # Summer flood date

wipeoutwindow.AMNA <-
  wipeoutwindow.AMNA_func(flowdata$SuFloodDate) # Summer flood date

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
#flood <- SP_highflood + medflood # simply whether it's a flood year or not
  ## Four years with SP_highflood and SU_highflood

# ITERATION PARAMETERS -------------------------------------------------------------------
# Setting up arrays/vectors to fill with data from loops

# Mid loop details -----------------------------------------------------------------------
# 'iterations' - number of replicate flow sequences to run for averaging, SE, etc.

years <- as.character(seq(1994, 2008, by = 1))
CACLrep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
GIROrep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
LECYrep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
CAINrep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
MIDOrep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
CYLUrep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
AMNArep <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))  
Total.N <- array(0, dim = c(15, iterations), dimnames = list(years, 1:iterations))

# Inner loop details ---------------------------------------------------------------------
# 'count' - number of years to project simulations (inner loop)

# Output of biomass and no. ind. for each age class for each year projected  
# An array with 3 columns (each stage class) and however many rows there are years projected 

CACLoutput.N <- array(0, dim = c(count, 3))
CACLoutput.biom <- array(0, dim = c(count, 3))
CACL.lambda <- array(0, dim = c(count,1))

GIROoutput.N <- array(0, dim = c(count, 3))
GIROoutput.biom <- array(0, dim = c(count, 3))
GIRO.lambda <- array(0, dim = c(count,1))

LECYoutput.N <- array(0, dim = c(count, 3))
LECYoutput.biom <- array(0, dim = c(count, 3))
LECY.lambda <- array(0, dim = c(count,1))

CAINoutput.N <- array(0, dim = c(count, 3))
CAINoutput.biom <- array(0, dim = c(count, 3))
CAIN.lambda <- array(0, dim = c(count,1))

MIDOoutput.N <- array(0, dim = c(count, 3))
MIDOoutput.biom <- array(0, dim = c(count, 3))
MIDO.lambda <- array(0, dim = c(count,1))

CYLUoutput.N <- array(0, dim = c(count, 3))
CYLUoutput.biom <- array(0, dim = c(count, 3))
CYLU.lambda <- array(0, dim = c(count,1))

AMNAoutput.N <- array(0, dim = c(count, 3))
AMNAoutput.biom <- array(0, dim = c(count, 3))
AMNA.lambda <- array(0, dim = c(count,1))

# Total biomass as % of K 
CACLbiomoutput <- numeric(length = count)
GIRObiomoutput <- numeric(length = count)
LECYbiomoutput <- numeric(length = count)
CAINbiomoutput <- numeric(length = count)
MIDObiomoutput <- numeric(length = count)
CYLUbiomoutput <- numeric(length = count)
AMNAbiomoutput <- numeric(length = count)

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
FCAINoutput <- numeric(length = count)
FMIDOoutput <- numeric(length = count)
FCYLUoutput <- numeric(length = count)
FAMNAoutput <- numeric(length = count)



# Mid loop ###############################################################################
# Middle loop uses iterator "iter" to get "iterations" for suming S2 and S3
for(iter in 1:iterations) {
  
  z <- c(sample(nrow(flowdata), burnin), 1:54) # mid loop - annual flow sequence over which fish model will run
 
# Need to read in initial biom every time so starting biomass is reset to be equal across spp and stages
  # N gives the total number of individuals for each age class.
  # Initially here, this is found by multiplying the number of g occupied by a given class
  # by the density per g
  # biom = g/m3
  # den = indiv/g
  biomCACL <- c(biomCACL1,
                biomCACL2,
                biomCACL3) 
  
  biomGIRO <- c(biomGIRO1,
                biomGIRO2,
                biomGIRO3)
  
  biomLECY <- c(biomLECY1,
                biomLECY2,
                biomLECY3) 
  
  biomCAIN <- c(biomCAIN1,
                biomCAIN2,
                biomCAIN3) 
  
  biomMIDO <- c(biomMIDO1,
                biomMIDO2,
                biomMIDO3) 
  
  biomCYLU <- c(biomCYLU1,
                biomCYLU2,
                biomCYLU3) 
  
  biomAMNA <- c(biomAMNA1,
                biomAMNA2,
                biomAMNA3) 
  
# Inner loop #############################################################################
for(i in 1:count) {

# CHANGE WHAT 'y' IS TO SIMULATE DIFFERENT FLOW REGIMES ACROSS THE 54 YEARS (+ BURNIN YEARS)
 #y = sample(nrow(flowdata), 1) 
 #y = i
 y = z[i] # sample burnin + 54 years

# y was a random number within the length of the flow data to randomly select a year from 
# the 'flood' and 'drought' vector. now it's directly taken from flow vector. 

# VITAL RATE DEFINITIONS: Desert Sucker  -----------------------------------------------------
# G is prob. of transition to next stage
# P is prob. of remaining in that stage

# flood is a vector of either medium or high magnitude
# Summer mortality vital rate is multiplied by multiplier based on yeartype as specified above
# Different from plant model as we have built in multiplier of severity of droughts or non-events on transition probabilities 

# Natives
# Desert Sucker - 
# YOY (GCACL1) survival and recruitment depends on spring flows (Sp is modifier for stage 1, Su is modifier for stage 2&3)
    GCACL1 <- aCACL1 * denCACLJ * (1/denCACL2) *
        (1 - (SP_highflood[y] * S2MortCACL * CACL_Sp_HF)) *
        (1 - (medflood[y] * S2MortCACL * CACL_Sp_MF)) *
        (1 - (nonevent[y] * S2MortCACL* CACL_Sp_NE))  *
        (1 - (drought[y] * S2MortCACL* CACL_Sp_DR)) 
         
    GCACL2 <- aCACL2 * denCACL2 * (1/denCACL3) *
        (1 - (SP_highflood[y] * S2MortCACL* CACL_Su_HF)) *
        (1 - (medflood[y] * S2MortCACL * CACL_Su_MF))*
        (1 - (nonevent[y] * S2MortCACL * CACL_Su_NE))  *
        (1 - (drought[y] * S2MortCACL * CACL_Su_DR))  

    PCACL3 <- (1 - aCACL3) *
        (1 - (SP_highflood[y] * S3MortCACL * CACL_Su_HF)) *
        (1 - (medflood[y] * S3MortCACL * CACL_Su_MF))  *
        (1 - (nonevent[y] * S3MortCACL * CACL_Su_NE))  *
        (1 - (drought[y] * S3MortCACL * CACL_Su_DR)) 
    

# Chub
     GGIRO1 <- aGIRO1 * denGIROJ * (1/denGIRO2) *
        (1 - (SP_highflood[y] * S2MortGIRO * GIRO_Sp_HF)) *
        (1 - (medflood[y] * S2MortGIRO * GIRO_Sp_MF)) *
        (1 - (nonevent[y] * S2MortGIRO * GIRO_Sp_NE))  *
        (1 - (drought[y] * S2MortGIRO * GIRO_Sp_DR)) 
        
    GGIRO2 <- aGIRO2 * denGIRO2 * (1/denGIRO3) *
        (1 - (SP_highflood[y] * S2MortGIRO)) * GIRO_Su_HF *
        (1 - (medflood[y] * S2MortGIRO * GIRO_Su_MF)) *
        (1 - (nonevent[y] * S2MortGIRO * GIRO_Su_NE)) *
        (1 - (drought[y] * S2MortGIRO * GIRO_Su_DR)) 

    PGIRO3 <- (1 - aGIRO3) *
        (1 - (SP_highflood[y] * S3MortGIRO * GIRO_Su_HF)) *
        (1 - (medflood[y] * S3MortGIRO * GIRO_Su_MF)) *
        (1 - (nonevent[y] * S3MortGIRO * GIRO_Su_NE)) *
        (1 - (drought[y] * S3MortGIRO * GIRO_Su_DR))
    
# Sonora sucker
    GCAIN1 <- aCAIN1 * denCAINJ * (1/denCAIN2) *
      (1 - (SP_highflood[y] * S2MortCAIN * CAIN_Sp_HF)) *
      (1 - (medflood[y] * S2MortCAIN * CAIN_Sp_MF)) *
      (1 - (nonevent[y] * S2MortCAIN * CAIN_Sp_NE)) *
      (1 - (drought[y] * S2MortCAIN * CAIN_Sp_DR))
    
    GCAIN2 <- aCAIN2 * denCAIN2 * (1/denCAIN3) *
      (1 - (SP_highflood[y] * S2MortCAIN * CAIN_Su_HF)) *
      (1 - (medflood[y] * S2MortCAIN * CAIN_Su_MF)) *
      (1 - (nonevent[y] * S2MortCAIN * CAIN_Su_NE)) *
      (1 - (drought[y] * S2MortCAIN * CAIN_Su_DR))
    
    PCAIN3 <- (1 - aCAIN3) *
      (1 - (SP_highflood[y] * S3MortCAIN * CAIN_Su_HF)) *
      (1 - (medflood[y] * S3MortCAIN * CAIN_Su_MF)) *
      (1 - (nonevent[y] * S3MortCAIN * CAIN_Su_NE)) *
      (1 - (drought[y] * S3MortCAIN * CAIN_Su_DR))
    
# Non-natives
    # note: in years with a spring and summer flood they get 2 x mortality
# Green sunfish - note the "NN"
 GLECY1 <- aLECY1 * denLECYJ * (1/denLECY2) *
        (1 - (SU_highflood[y] * S2MortLECY * LECY_Sp_HF)) *
        (1 - (medflood[y] * S2MortLECY * LECY_Sp_MF)) *
        (1 - (nonevent[y] * S2MortLECY * LECY_Sp_NE)) *
        (1 - (drought[y] * S2MortLECY * LECY_Sp_DR))
        
    GLECY2 <- aLECY2 * denLECY2 * (1/denLECY3) *
        (1 - (SU_highflood[y] * S2MortLECY * LECY_Su_HF)) *
        (1 - (medflood[y] * S2MortLECY * LECY_Su_MF)) *
        (1 - (nonevent[y] * S2MortLECY * LECY_Su_NE)) *
        (1 - (drought[y] * S2MortLECY * LECY_Su_DR))

    PLECY3 <- (1 - aLECY3) *
        (1 - (SU_highflood[y] * S3MortLECY * LECY_Su_HF)) *
        (1 - (medflood[y] * S3MortLECY * LECY_Su_MF)) *
        (1 - (nonevent[y] * S3MortLECY * LECY_Su_NE)) *
        (1 - (drought[y] * S3MortLECY * LECY_Su_DR))
    
# Smallmouth bass
    GMIDO1 <- aMIDO1 * denMIDOJ * (1/denMIDO2) *
      (1 - (SU_highflood[y] * S2MortMIDO * MIDO_Sp_HF)) *
      (1 - (medflood[y] * S2MortMIDO * MIDO_Sp_MF)) *
      (1 - (nonevent[y] * S2MortMIDO * MIDO_Sp_NE)) *
      (1 - (drought[y] * S2MortMIDO * MIDO_Sp_DR))
    
    GMIDO2 <- aMIDO2 * denMIDO2 * (1/denMIDO3) *
      (1 - (SU_highflood[y] * S2MortMIDO * MIDO_Su_HF)) *
      (1 - (medflood[y] * S2MortMIDO * MIDO_Su_MF)) *
      (1 - (nonevent[y] * S2MortMIDO * MIDO_Su_NE)) *
      (1 - (drought[y] * S2MortMIDO * MIDO_Su_DR))
    
    PMIDO3 <- (1 - aMIDO3) *
      (1 - (SU_highflood[y] * S3MortMIDO * MIDO_Su_HF)) *
      (1 - (medflood[y] * S3MortMIDO * MIDO_Su_MF)) *
      (1 - (nonevent[y] * S3MortMIDO * MIDO_Su_NE)) *
      (1 - (drought[y] * S3MortMIDO * MIDO_Su_DR))
    
# Red shiner
    GCYLU1 <- aCYLU1 * denCYLUJ * (1/denCYLU2) *
      (1 - (SU_highflood[y] * S2MortCYLU * CYLU_Sp_HF)) *
      (1 - (medflood[y] * S2MortCYLU * CYLU_Sp_MF)) *
      (1 - (nonevent[y] * S2MortCYLU * CYLU_Sp_NE)) *
      (1 - (drought[y] * S2MortCYLU * CYLU_Sp_DR))
    
    GCYLU2 <- aCYLU2 * denCYLU2 * (1/denCYLU3) *
      (1 - (SU_highflood[y] * S2MortCYLU * CYLU_Su_HF)) *
      (1 - (medflood[y] * S2MortCYLU * CYLU_Su_MF)) *
      (1 - (nonevent[y] * S2MortCYLU * CYLU_Su_NE)) *
      (1 - (drought[y] * S2MortCYLU * CYLU_Su_DR))
    
    PCYLU3 <- (1 - aCYLU3) *
      (1 - (SU_highflood[y] * S3MortCYLU * CYLU_Su_HF)) *
      (1 - (medflood[y] * S3MortCYLU * CYLU_Su_MF)) *
      (1 - (nonevent[y] * S3MortCYLU * CYLU_Su_NE)) *
      (1 - (drought[y] * S3MortCYLU * CYLU_Su_DR))
    
# Yellow bullhead
    GAMNA1 <- aAMNA1 * denAMNAJ * (1/denAMNA2) *
      (1 - (SU_highflood[y] * S2MortAMNA * AMNA_Sp_HF)) *
      (1 - (medflood[y] * S2MortAMNA * AMNA_Sp_MF)) *
      (1 - (nonevent[y] * S2MortAMNA * AMNA_Sp_NE)) *
      (1 - (drought[y] * S2MortAMNA * AMNA_Sp_DR))
    
    GAMNA2 <- aAMNA2 * denAMNA2 * (1/denAMNA3) *
      (1 - (SU_highflood[y] * S2MortAMNA * AMNA_Su_HF)) *
      (1 - (medflood[y] * S2MortAMNA * AMNA_Su_MF)) *
      (1 - (nonevent[y] * S2MortAMNA * AMNA_Su_NE)) *
      (1 - (drought[y] * S2MortAMNA * AMNA_Su_DR))
    
    PAMNA3 <- (1 - aAMNA3) *
      (1 - (SU_highflood[y] * S3MortAMNA * AMNA_Su_HF)) *
      (1 - (medflood[y] * S3MortAMNA * AMNA_Su_MF)) *
      (1 - (nonevent[y] * S3MortAMNA * AMNA_Su_NE)) *
      (1 - (drought[y] * S3MortAMNA * AMNA_Su_DR))
    
 
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

totbiom.CAIN <- 
  biomCAIN[1] + #
  biomCAIN[2] + 
  biomCAIN[3] 

totbiom.MIDO <- 
  biomMIDO[1] + #
  biomMIDO[2] + 
  biomMIDO[3] 

totbiom.CYLU <- 
  biomCYLU[1] + #
  biomCYLU[2] + 
  biomCYLU[3] 

totbiom.AMNA <- 
  biomAMNA[1] + #
  biomAMNA[2] + 
  biomAMNA[3] 

### :CLARIFY: at the moment carrying capacity is limiting spawning of all species based on the total biomass occupied at the end of the previous year. i.e. if above K, no spp spawn in that year. If spawning occurs, they can all spawn and there is no sequence, so overseeding COULD be massive.
b <- 1/47660
# POTENTIAL CACL FECUNDITY ---------------------------------------------------------
FCACL2 <- ((0.5*GSI.CACL*(1-S1MortCACL))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denCACL1*(1/denCACLJ)
FCACL3 <- ((0.5*GSI.CACL*(1-S1MortCACL))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denCACL1*(1/denCACLJ)


# POTENTIAL GIRO FECUNDITY ---------------------------------------------------------
FGIRO2 <- ((0.5*GSI.GIRO*(1-S1MortGIRO))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denGIRO1*(1/denGIROJ)
FGIRO3 <- ((0.5*GSI.GIRO*(1-S1MortGIRO))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denGIRO1*(1/denGIROJ)

# POTENTIAL CAIN FECUNDITY ---------------------------------------------------------

FCAIN3 <- ((0.5*GSI.CAIN*(1-S1MortCAIN))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denCAIN1*(1/denCAINJ)

# POTENTIAL LECY FECUNDITY ---------------------------------------------------------
FLECY2 <- ((0.5*GSI.LECY*(1-S1MortLECY))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denLECY1*(1/denLECYJ)
FLECY3 <- ((0.5*GSI.LECY*(1-S1MortLECY))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denLECY1*(1/denLECYJ)
# POTENTIAL MIDO FECUNDITY ---------------------------------------------------------

FMIDO3 <- ((0.5*GSI.MIDO*(1-S1MortMIDO))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denMIDO1*(1/denMIDOJ)

# POTENTIAL CYLU FECUNDITY ---------------------------------------------------------
# because they are serial spawners, they are allowed to spawn twice a season
FCYLUJ <- ((0.5*GSI.CYLU*(1-S1MortCYLU))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denCYLU1*(1/denCYLUJ)
FCYLU2 <- ((0.5*2*GSI.CYLU*(1-S1MortCYLU))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denCYLU1*(1/denCYLUJ)
FCYLU3 <- ((0.5*2*GSI.CYLU*(1-S1MortCYLU))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denCYLU1*(1/denCYLUJ)

# POTENTIAL AMNA FECUNDITY ---------------------------------------------------------

FAMNA3 <- ((0.5*GSI.AMNA*(1-S1MortAMNA))/
             (1+(b*sum(totbiom.CACL, totbiom.GIRO, totbiom.LECY, totbiom.CAIN, totbiom.MIDO, totbiom.CYLU, totbiom.AMNA))))*
  denAMNA1*(1/denAMNAJ)


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

# CAIN
# gives total CAIN biomass as a percentage of K; 
KCAIN <- 100 * (biomCAIN[1] + 
                  biomCAIN[2] + 
                  biomCAIN[3])/K

# MIDO
# gives total MIDO biomass as a percentage of K; 
KMIDO <- 100 * (biomMIDO[1] + 
                  biomMIDO[2] + 
                  biomMIDO[3])/K

# CYLU
# gives total CYLU biomass as a percentage of K; 
KCYLU <- 100 * (biomCYLU[1] + 
                  biomCYLU[2] + 
                  biomCYLU[3])/K

# AMNA
# gives total CYLU biomass as a percentage of K; 
KAMNA <- 100 * (biomAMNA[1] + 
                  biomAMNA[2] + 
                  biomAMNA[3])/K
    
# TRANSITION MATRICES --------------------------------------------------------------------

# TRANSITION MATRIX FOR CACL -------------------------------------------------------
ACACL1 <- c(0, FCACL2, FCACL3)
ACACL2 <- c(GCACL1, 0, 0)
ACACL3 <- c(0, GCACL2, PCACL3)
# Matrix
ACACL <- rbind(ACACL1, ACACL2, ACACL3)
lambda(ACACL) # Checking population growth rate
# TRANSITION MATRIX FOR GIRO -------------------------------------------------------
AGIRO1 <- c(0, FGIRO2, FGIRO3)
AGIRO2 <- c(GGIRO1, 0, 0)
AGIRO3 <- c(0, GGIRO2, PGIRO3)
# Matrix
AGIRO <- rbind(AGIRO1, AGIRO2, AGIRO3)
lambda(AGIRO) # Checking population growth rate

# TRANSITION MATRIX FOR LECY -------------------------------------------------------
ALECY1 <- c(0, FLECY2, FLECY3)
ALECY2 <- c(GLECY1, 0, 0)
ALECY3 <- c(0, GLECY2, PLECY3)
# Matrix
ALECY <- rbind(ALECY1, ALECY2, ALECY3)
lambda(ALECY) # Checking population growth rate

# TRANSITION MATRIX FOR CAIN -------------------------------------------------------
ACAIN1 <- c(0, 0, FCAIN3)
ACAIN2 <- c(GCAIN1, 0, 0)
ACAIN3 <- c(0, GCAIN2, PCAIN3)
# Matrix
ACAIN <- rbind(ACAIN1, ACAIN2, ACAIN3)
lambda(ACAIN) # Checking population growth rate

# TRANSITION MATRIX FOR MIDO -------------------------------------------------------
AMIDO1 <- c(0, 0, FMIDO3)
AMIDO2 <- c(GMIDO1, 0, 0)
AMIDO3 <- c(0, GMIDO2, PMIDO3)
# Matrix
AMIDO <- rbind(AMIDO1, AMIDO2, AMIDO3)
lambda(AMIDO) # Checking population growth rate

# TRANSITION MATRIX FOR CYLU -------------------------------------------------------
ACYLU1 <- c(FCYLUJ, FCYLU2, FCYLU3)
ACYLU2 <- c(GCYLU1, 0, 0)
ACYLU3 <- c(0, GCYLU2, PCYLU3)
# Matrix
ACYLU <- rbind(ACYLU1, ACYLU2, ACYLU3)
lambda(ACYLU) # Checking population growth rate

# TRANSITION MATRIX FOR AMNA -------------------------------------------------------
AAMNA1 <- c(0, 0, FAMNA3)
AAMNA2 <- c(GAMNA1, 0, 0)
AAMNA3 <- c(0, GAMNA2, PAMNA3)
# Matrix
AAMNA <- rbind(AAMNA1, AAMNA2, AAMNA3)
lambda(AAMNA) # Checking population growth rate

# COMPILING OUTPUTS ----------------------------------------------------------------------
# CACL
CACLoutput.biom[i,1:3] <- biomCACL # array of biomass of each age class for each yr projected. biomCACL = total biomass for each age class
#CACLbiomoutput[i] <- KCACL # total biomass as % of K; this is in g
FCACLoutput[i] <- FCACL3 + FCACL2
CACLoutput.N [i, 1:3] <- 
  c(biomCACL[1] * denCACLJ,
    biomCACL[2] * denCACL2,
    biomCACL[3] * denCACL3)
CACL.lambda[i] <- lambda(ACACL)
    
# GIRO
GIROoutput.biom[i,1:3] <- biomGIRO # array of biomass of each age class for each yr projected. biomGIRO = total biomass for each age class
#GIRObiomoutput[i] <- KGIRO # total biomass as % of K; this is in g
FGIROoutput[i] <- FGIRO3 + FGIRO2
GIROoutput.N [i, 1:3] <- 
  c(biomGIRO[1] * denGIROJ,
    biomGIRO[2] * denGIRO2,
    biomGIRO[3] * denGIRO3)
GIRO.lambda[i] <- lambda(AGIRO)   

# LECY
LECYoutput.biom[i,1:3] <- biomLECY # array of biomass of each age class for each yr projected. biomLECY = total biomass for each age class
#LECYbiomoutput[i] <- KLECY # total biomass as % of K; this is in g 
FLECYoutput[i] <- FLECY3 + FLECY2
LECYoutput.N [i, 1:3] <- 
  c(biomLECY[1] * denLECYJ,
    biomLECY[2] * denLECY2,
    biomLECY[3] * denLECY3)
LECY.lambda[i] <- lambda(ALECY)

# CAIN
CAINoutput.biom[i,1:3] <- biomCAIN # array of biomass of each age class for each yr projected. biomCAIN = total biomass for each age class
#CAINbiomoutput[i] <- KCAIN # total biomass as % of K; this is in g 
FCAINoutput[i] <- FCAIN3
CAINoutput.N [i, 1:3] <- 
  c(biomCAIN[1] * denCAINJ,
    biomCAIN[2] * denCAIN2,
    biomCAIN[3] * denCAIN3)
CAIN.lambda[i] <- lambda(ACAIN)

# MIDO
MIDOoutput.biom[i,1:3] <- biomMIDO # array of biomass of each age class for each yr projected. biomMIDO = total biomass for each age class
#MIDObiomoutput[i] <- KMIDO # total biomass as % of K; this is in g 
FMIDOoutput[i] <- FMIDO3
MIDOoutput.N [i, 1:3] <- 
  c(biomMIDO[1] * denMIDOJ,
    biomMIDO[2] * denMIDO2,
    biomMIDO[3] * denMIDO3)
MIDO.lambda[i] <- lambda(AMIDO)

# CYLU
CYLUoutput.biom[i,1:3] <- biomCYLU # array of biomass of each age class for each yr projected. biomCYLU = total biomass for each age class
#CYLUbiomoutput[i] <- KCYLU # total biomass as % of K; this is in g 
FCYLUoutput[i] <- FCYLU3 + FCYLU2 + FCYLUJ
CYLUoutput.N [i, 1:3] <- 
  c(biomCYLU[1] * denCYLUJ,
    biomCYLU[2] * denCYLU2,
    biomCYLU[3] * denCYLU3)
CYLU.lambda[i] <- lambda(ACYLU)

# AMNA
AMNAoutput.biom[i,1:3] <- biomAMNA # array of biomass of each age class for each yr projected. biomAMNA = total biomass for each age class
#AMNAbiomoutput[i] <- KAMNA # total biomass as % of K; this is in g 
FAMNAoutput[i] <- FAMNA3
AMNAoutput.N [i, 1:3] <- 
  c(biomAMNA[1] * denAMNAJ,
    biomAMNA[2] * denAMNA2,
    biomAMNA[3] * denAMNA3)
AMNA.lambda[i] <- lambda(AAMNA)

# Records flood settings of each particular projected year (0 for nonflood, 1 for flood)
SPhighfloodoutput[i] <- SP_highflood[y]
SUhighfloodoutput[i] <- SU_highflood[y]
medfloodoutput[i] <- medflood[y]
noneventoutput[i] <- nonevent[y]
droughtoutput[i] <- drought[y]

# if we want to add total fecundity in afterwards
# CACLplaceholder <- FCACL3
# GIROplaceholder <- FGIRO3
# LECYplaceholder <- FLECY3
    
# MATRIX MULTIPLICATION ------------------------------------------------------------------
# can include rescue function for each with 0.5 chance of reach being colonized by 2 individuals
# CACL
biomCACL <- ACACL %*% biomCACL #biom_res(biomCACL, "CACL") # ACACL is transition matrix, biomCACL = total biomass for each age class,

# GIRO
biomGIRO <- AGIRO %*% biomGIRO #biom_res(biomGIRO, "GIRO") # AGIRO is transition matrix, biomGIRO = total biomass for each age class

# LECY
biomLECY <- ALECY %*% biomLECY #biom_res(biomLECY, "LECY") # ALECY is transition matrix, biomLECY = total biomass for each age class

# CAIN
biomCAIN <- ACAIN %*% biomCAIN #biom_res(biomCAIN, "CAIN") # ACAIN is transition matrix, biomCAIN = total biomass for each age class

# MIDO
biomMIDO <- AMIDO %*% biomMIDO #biom_res(biomMIDO, "MIDO") # AMIDO is transition matrix, biomMIDO = total biomass for each age class

# CYLU
biomCYLU <- ACYLU %*% biomCYLU #biom_res(biomCYLU, "CYLU") # ACYLU is transition matrix, biomCYLU = total biomass for each age class

# AMNA
biomAMNA <- AAMNA %*% biomAMNA #biom_res(biomAMNA, "AMNA") # AAMNA is transition matrix, biomAMNA = total biomass for each age class

} # End of inner loop ####################################################################

# Mean values for each iteration run over each sequence of years
CACLrep[,iter] <- apply(CACLoutput.N[61:75, 2:3], 1, sum)
GIROrep[,iter] <- apply(GIROoutput.N[61:75, 2:3], 1, sum)
LECYrep[,iter] <- apply(LECYoutput.N[61:75, 2:3], 1, sum)
CAINrep[,iter] <- apply(CAINoutput.N[61:75, 2:3], 1, sum)
MIDOrep[,iter] <- apply(MIDOoutput.N[61:75, 2:3], 1, sum)
CYLUrep[,iter] <- apply(CYLUoutput.N[61:75, 2:3], 1, sum)
AMNArep[,iter] <- apply(AMNAoutput.N[61:75, 2:3], 1, sum)
Total.N[,iter] <- apply(
  cbind(CACLoutput.N[61:75, 2:3], GIROoutput.N[61:75, 2:3], LECYoutput.N[61:75, 2:3], CAINoutput.N[61:75, 2:3], 
        MIDOoutput.N[61:75, 2:3], CYLUoutput.N[61:75, 2:3], AMNAoutput.N[61:75, 2:3]), 1, sum)

} # End of mid loop ######################################################################

# 1000 iterations takes less than 9 mins
# ########################################################################################
# OUTPUTS --------------------------------------------------------------------------------
# ########################################################################################

CACLoutput.biom.DF <- as.data.frame(CACLoutput.biom) %>%
    rename(S1 = V1, S2 = V2, S3 = V3) %>%
    mutate(rep = row.names(.)) %>%
    gather(stage, g, -rep) %>%
    mutate(spp = 'CACL') 

CACLoutput.N.DF <- as.data.frame(CACLoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'CACL')
    
GIROoutput.biom.DF <- as.data.frame(GIROoutput.biom) %>%
    rename(S1 = V1, S2 = V2, S3 = V3) %>%
    mutate(rep = row.names(.)) %>%
    gather(stage, g, -rep) %>%
    mutate(spp = 'GIRO')

GIROoutput.N.DF <- as.data.frame(GIROoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'GIRO')

LECYoutput.biom.DF <- as.data.frame(LECYoutput.biom) %>%
    rename(S1 = V1, S2 = V2, S3 = V3) %>%
    mutate(rep = row.names(.)) %>%
    gather(stage, g, -rep) %>%
    mutate(spp = 'LECY')

LECYoutput.N.DF <- as.data.frame(LECYoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'LECY')

CAINoutput.biom.DF <- as.data.frame(CAINoutput.biom) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, g, -rep) %>%
  mutate(spp = 'CAIN')

CAINoutput.N.DF <- as.data.frame(CAINoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'CAIN')

MIDOoutput.biom.DF <- as.data.frame(MIDOoutput.biom) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, g, -rep) %>%
  mutate(spp = 'MIDO')

MIDOoutput.N.DF <- as.data.frame(MIDOoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'MIDO')

CYLUoutput.biom.DF <- as.data.frame(CYLUoutput.biom) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, g, -rep) %>%
  mutate(spp = 'CYLU')

CYLUoutput.N.DF <- as.data.frame(CYLUoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'CYLU')

AMNAoutput.biom.DF <- as.data.frame(AMNAoutput.biom) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, g, -rep) %>%
  mutate(spp = 'AMNA')

AMNAoutput.N.DF <- as.data.frame(AMNAoutput.N) %>%
  rename(S1 = V1, S2 = V2, S3 = V3) %>%
  mutate(rep = row.names(.)) %>%
  gather(stage, N, -rep) %>%
  mutate(spp = 'AMNA')

ALLoutput.N.DF <- rbind(CACLoutput.N.DF, GIROoutput.N.DF, LECYoutput.N.DF, CAINoutput.N.DF, MIDOoutput.N.DF, CYLUoutput.N.DF, AMNAoutput.N.DF)  

ALLoutput.biom.DF <- rbind(CACLoutput.biom.DF, GIROoutput.biom.DF, LECYoutput.biom.DF, CAINoutput.biom.DF, MIDOoutput.biom.DF, CYLUoutput.biom.DF, AMNAoutput.biom.DF)

head(ALLoutput.biom.DF)

# Graph biomass
ggplot(ALLoutput.biom.DF, aes(as.numeric(rep), g, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(stage~spp, scales = "free")

# Graph abundance
ggplot(ALLoutput.N.DF, aes(as.numeric(rep), N, colour = stage)) +
  geom_point() +
  geom_path() +
  facet_grid(stage~spp, scales = "free")

# Graph flows
# ggplot(flowdata, aes(Year, SpFloodMag)) +
#     geom_point() + geom_path()
# ggplot(flowdata, aes(Year, BaseDur)) +
#     geom_point() + geom_path()

flowresults <- as.data.frame(cbind(SPhighfloodoutput, SUhighfloodoutput, medfloodoutput, noneventoutput, droughtoutput)) %>%
    mutate(rep = as.numeric(row.names(.))) %>%
    gather(metric, value, -rep)

ggplot(flowresults, aes(rep, value)) +
    geom_point() + geom_path() +
    facet_wrap(~metric)

# Graph all species together
ggplot(ALLoutput.biom.DF, aes(as.numeric(rep), g, colour = stage)) +
    geom_point() +
    geom_path() +
    facet_grid(~spp)

ggplot(ALLoutput.N.DF, aes(as.numeric(rep), N, colour = stage)) + # [ALLoutput.N.DF$spp != "CYLU",]
  geom_point() +
  geom_path() +
  facet_grid(~spp)

# Relative abundance ------------------------------------------------------------------------------
Verde <- read.csv("data/Rel_Abu_Verde_94-08.csv", header = T)
par(mfrow = c(3,3), mar = c(4,4,1,0), oma=c(1,0,1,1))
plot(Verde$Year[Verde$SppCode == "CACL"], Verde$MeanRelAbu[Verde$SppCode == "CACL"], ylim = c(0, 0.5),
     ylab = "Relative Abundance", main = "CACL", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CACL"] - Verde$SERelAbu[Verde$SppCode == "CACL"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CACL"] + Verde$SERelAbu[Verde$SppCode == "CACL"], code = 3, length = 0.1, angle = 90)

# CACL
CACL.RelAbu <- CACLrep/Total.N
CACL.mean.RA <- apply(CACL.RelAbu, 1, mean)
CACL.sd.RA <- apply(CACL.RelAbu, 1, sd)
CACL.SE.RA <- CACL.sd.RA/sqrt(iterations)
CACL.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CACL"] - CACL.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
CACL.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CACL"] - CACL.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "CACL"]) - min(Verde$MeanRelAbu[Verde$SppCode == "CACL"]))

plot(Verde$Year[Verde$SppCode == "CACL"], Verde$MeanRelAbu[Verde$SppCode == "CACL"], ylim = c(0, 0.5),
     ylab = "Relative Abundance", main = "CACL", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CACL"] - Verde$SERelAbu[Verde$SppCode == "CACL"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CACL"] + Verde$SERelAbu[Verde$SppCode == "CACL"], code = 3, length = 0.1, angle = 90)
points(years, CACL.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = CACL.mean.RA + CACL.SE.RA, 
       x1 = c(1994:2008), y1 = CACL.mean.RA - CACL.SE.RA, code = 3, length = 0.1, angle = 90)

# GIRO
GIRO.RelAbu <- GIROrep/Total.N
GIRO.mean.RA <- apply(GIRO.RelAbu, 1, mean)
GIRO.sd.RA <- apply(GIRO.RelAbu, 1, sd)
GIRO.SE.RA <- GIRO.sd.RA/sqrt(iterations)
GIRO.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "GIRO"] - GIRO.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
GIRO.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "GIRO"] - GIRO.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "GIRO"]) - min(Verde$MeanRelAbu[Verde$SppCode == "GIRO"]))

plot(Verde$Year[Verde$SppCode == "GIRO"], Verde$MeanRelAbu[Verde$SppCode == "GIRO"], ylim = c(0, 0.5),
     ylab = "Relative Abundance", main = "GIRO", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "GIRO"] - Verde$SERelAbu[Verde$SppCode == "GIRO"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "GIRO"] + Verde$SERelAbu[Verde$SppCode == "GIRO"], code = 3, length = 0.1, angle = 90)
points(years, GIRO.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = GIRO.mean.RA + GIRO.SE.RA, 
       x1 = c(1994:2008), y1 = GIRO.mean.RA - GIRO.SE.RA, code = 3, length = 0.1, angle = 90)

# LECY
LECY.RelAbu <- LECYrep/Total.N
LECY.mean.RA <- apply(LECY.RelAbu, 1, mean)
LECY.sd.RA <- apply(LECY.RelAbu, 1, sd)
LECY.SE.RA <- LECY.sd.RA/sqrt(iterations)
LECY.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "LECY"] - LECY.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
LECY.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "LECY"] - LECY.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "LECY"]) - min(Verde$MeanRelAbu[Verde$SppCode == "LECY"]))

plot(Verde$Year[Verde$SppCode == "LECY"], Verde$MeanRelAbu[Verde$SppCode == "LECY"], ylim = c(0, 0.5),
     ylab = "Relative Abundance", main = "LECY", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "LECY"] - Verde$SERelAbu[Verde$SppCode == "LECY"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "LECY"] + Verde$SERelAbu[Verde$SppCode == "LECY"], code = 3, length = 0.1, angle = 90)
points(years, LECY.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = LECY.mean.RA + LECY.SE.RA, 
       x1 = c(1994:2008), y1 = LECY.mean.RA - LECY.SE.RA, code = 3, length = 0.1, angle = 90)

# CAIN
CAIN.RelAbu <- CAINrep/Total.N
CAIN.mean.RA <- apply(CAIN.RelAbu, 1, mean)
CAIN.sd.RA <- apply(CAIN.RelAbu, 1, sd)
CAIN.SE.RA <- CAIN.sd.RA/sqrt(iterations)
CAIN.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CAIN"] - CAIN.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
CAIN.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CAIN"] - CAIN.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "CAIN"]) - min(Verde$MeanRelAbu[Verde$SppCode == "CAIN"]))

plot(Verde$Year[Verde$SppCode == "CAIN"], Verde$MeanRelAbu[Verde$SppCode == "CAIN"], ylim = c(0, 0.5),
     ylab = "Relative Abundance", main = "CAIN", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CAIN"] - Verde$SERelAbu[Verde$SppCode == "CAIN"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CAIN"] + Verde$SERelAbu[Verde$SppCode == "CAIN"], code = 3, length = 0.1, angle = 90)
points(years, CAIN.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = CAIN.mean.RA + CAIN.SE.RA, 
       x1 = c(1994:2008), y1 = CAIN.mean.RA - CAIN.SE.RA, code = 3, length = 0.1, angle = 90)

# MIDO
MIDO.RelAbu <- MIDOrep/Total.N
MIDO.mean.RA <- apply(MIDO.RelAbu, 1, mean)
MIDO.sd.RA <- apply(MIDO.RelAbu, 1, sd)
MIDO.SE.RA <- MIDO.sd.RA/sqrt(iterations)
MIDO.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "MIDO"] - MIDO.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
MIDO.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "MIDO"] - MIDO.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "MIDO"]) - min(Verde$MeanRelAbu[Verde$SppCode == "MIDO"]))

plot(Verde$Year[Verde$SppCode == "MIDO"], Verde$MeanRelAbu[Verde$SppCode == "MIDO"], ylim = c(0, 1),
     ylab = "Relative Abundance", main = "MIDO", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "MIDO"] - Verde$SERelAbu[Verde$SppCode == "MIDO"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "MIDO"] + Verde$SERelAbu[Verde$SppCode == "MIDO"], code = 3, length = 0.1, angle = 90)
points(years, MIDO.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = MIDO.mean.RA + MIDO.SE.RA, 
       x1 = c(1994:2008), y1 = MIDO.mean.RA - MIDO.SE.RA, code = 3, length = 0.1, angle = 90)

# CYLU
CYLU.RelAbu <- CYLUrep/Total.N
CYLU.mean.RA <- apply(CYLU.RelAbu, 1, mean)
CYLU.sd.RA <- apply(CYLU.RelAbu, 1, sd)
CYLU.SE.RA <- CYLU.sd.RA/sqrt(iterations)
CYLU.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CYLU"] - CYLU.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
CYLU.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "CYLU"] - CYLU.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "CYLU"]) - min(Verde$MeanRelAbu[Verde$SppCode == "CYLU"]))

plot(Verde$Year[Verde$SppCode == "CYLU"], Verde$MeanRelAbu[Verde$SppCode == "CYLU"], ylim = c(0, 1),
     ylab = "Relative Abundance", main = "CYLU", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "CYLU"] - Verde$SERelAbu[Verde$SppCode == "CYLU"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "CYLU"] + Verde$SERelAbu[Verde$SppCode == "CYLU"], code = 3, length = 0.1, angle = 90)
points(years, CYLU.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = CYLU.mean.RA + CYLU.SE.RA, 
       x1 = c(1994:2008), y1 = CYLU.mean.RA - CYLU.SE.RA, code = 3, length = 0.1, angle = 90)

# AMNA
AMNA.RelAbu <- AMNArep/Total.N
AMNA.mean.RA <- apply(AMNA.RelAbu, 1, mean)
AMNA.sd.RA <- apply(AMNA.RelAbu, 1, sd)
AMNA.SE.RA <- AMNA.sd.RA/sqrt(iterations)
AMNA.RMSE <- sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "AMNA"] - AMNA.mean.RA)^2)) # RMSE: sqrt(mean((y-y_pred)^2))
AMNA.NRMSE <- (sqrt(mean((Verde$MeanRelAbu[Verde$SppCode == "AMNA"] - AMNA.mean.RA)^2)))/(max(Verde$MeanRelAbu[Verde$SppCode == "AMNA"]) - min(Verde$MeanRelAbu[Verde$SppCode == "AMNA"]))

plot(Verde$Year[Verde$SppCode == "AMNA"], Verde$MeanRelAbu[Verde$SppCode == "AMNA"], ylim = c(0, 0.5),
     ylab = "Relative Abundance", main = "AMNA", xlab = "Year")
arrows(x0 = c(1994:2008), y0 = Verde$MeanRelAbu[Verde$SppCode == "AMNA"] - Verde$SERelAbu[Verde$SppCode == "AMNA"], 
       x1 = c(1994:2008), y1 = Verde$MeanRelAbu[Verde$SppCode == "AMNA"] + Verde$SERelAbu[Verde$SppCode == "AMNA"], code = 3, length = 0.1, angle = 90)
points(years, AMNA.mean.RA, pch = 19)
arrows(x0 = c(1994:2008), y0 = AMNA.mean.RA + AMNA.SE.RA, 
       x1 = c(1994:2008), y1 = AMNA.mean.RA - AMNA.SE.RA, code = 3, length = 0.1, angle = 90)

# Correlation tests
Model.RelAbu <- rbind(AMNA.mean.RA, CACL.mean.RA, CAIN.mean.RA, CYLU.mean.RA, GIRO.mean.RA, LECY.mean.RA, MIDO.mean.RA)

#filter(Verde, Year == 1994)

spear.cor <- c(
cor(filter(Verde, Year == 1994)$MeanRelAbu, Model.RelAbu[,"1994"], method = "spear"),
cor(filter(Verde, Year == 1995)$MeanRelAbu, Model.RelAbu[,"1995"], method = "spear"),
cor(filter(Verde, Year == 1996)$MeanRelAbu, Model.RelAbu[,"1996"], method = "spear"),
cor(filter(Verde, Year == 1997)$MeanRelAbu, Model.RelAbu[,"1997"], method = "spear"),
cor(filter(Verde, Year == 1998)$MeanRelAbu, Model.RelAbu[,"1998"], method = "spear"),
cor(filter(Verde, Year == 1999)$MeanRelAbu, Model.RelAbu[,"1999"], method = "spear"),
cor(filter(Verde, Year == 2000)$MeanRelAbu, Model.RelAbu[,"2000"], method = "spear"),
cor(filter(Verde, Year == 2001)$MeanRelAbu, Model.RelAbu[,"2001"], method = "spear"),
cor(filter(Verde, Year == 2002)$MeanRelAbu, Model.RelAbu[,"2002"], method = "spear"),
cor(filter(Verde, Year == 2003)$MeanRelAbu, Model.RelAbu[,"2003"], method = "spear"),
cor(filter(Verde, Year == 2004)$MeanRelAbu, Model.RelAbu[,"2004"], method = "spear"),
cor(filter(Verde, Year == 2005)$MeanRelAbu, Model.RelAbu[,"2005"], method = "spear"),
cor(filter(Verde, Year == 2006)$MeanRelAbu, Model.RelAbu[,"2006"], method = "spear"),
cor(filter(Verde, Year == 2007)$MeanRelAbu, Model.RelAbu[,"2007"], method = "spear"),
cor(filter(Verde, Year == 2008)$MeanRelAbu, Model.RelAbu[,"2008"], method = "spear"))
spear <- cbind(Verde$Year[1:15], spear.cor)
round(spear, digits = 2)
pears <- cbind(Verde$Year[1:15], pearson.cor)
round(pears, digits = 2)


tail(AMNAoutput.N.DF)
tail(MIDOoutput.N.DF)
tail(LECYoutput.N.DF)
tail(GIROoutput.N.DF)
tail(CYLUoutput.N.DF)

head(AMNAoutput.N.DF)
head(MIDOoutput.N.DF)
head(LECYoutput.N.DF)
head(GIROoutput.N.DF)
head(CYLUoutput.N.DF)



