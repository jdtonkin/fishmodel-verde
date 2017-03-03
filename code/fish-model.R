
# Required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# SETUP ----------------------------------------------------------------------------------

rm(list = ls()) # clearing the workspace 

# bringing in flow data
# Verde flow data at Paulden 1984-2013, 30 years continuous
flowdata <- read.csv("flowdata_Verde.csv") 
str(flowdata)
head(flowdata)
# floodmag - magnitude of flood in cfs
# basedur - baseflow duration in days
# flooddate is peak dates of all floods (Oct 1 = 1)

count <- length(flowdata$floodmag) # inner loop - number of years to simulate; starting with natural sequence of flow record
# burnin = 100 # number of years to discard as burn in during long term mean estimation 
# outerreps <- 1 # number of iterations for outer loop that alters drought/flood frequency 
# replicates <- 100 # number of replicate projections to run (mid loop)

# Modifiers
modifiers <- read.csv('modifiers.csv')

# adding 'Modifier' value from csv to 'Code' in csv
for(j in 1:length(modifiers[,1])) {
    nam <- paste(modifiers[j,4])
    assign(nam, modifiers[j,6]) # SHOULD BE 5 FOR REAL VALUES!!!!
}

# Key ------------------------------------------------------------
# CACL: desert sucker
# GIRO: chub
# LECY: green sunfish

# Total volume of water in m3: 317
# Total fish biomass in g (with main spp): 3809
# 
#

K = 3809 # avg. for 100-m reach across 6 replicate reaches... in g/m3 this is 15 g/m3

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
    ifelse(flowdata$flooddate > 124 & flowdata$flooddate < 213, 1, 0)
}

spawnwindow.GIRO_func <- function(x) {
    ifelse(flowdata$flooddate > 153 & flowdata$flooddate < 244, 1, 0)
}

wipeoutwindow.LECY_func <- function(x) {
    ifelse(flowdata$flooddate > 198, 1, 0) # late floods stop spawning
}


spawnwindow.CACL <-
    spawnwindow.CACL_func(flowdata$flooddate)

spawnwindow.GIRO <-
    spawnwindow.GIRO_func(flowdata$flooddate)

wipeoutwindow.LECY <-
    wipeoutwindow.LECY_func(flowdata$flooddate)


# FLOOD THRESHOLD FUNCTIONS --------------------------------------------------------------
# flowdata$floodmag - vector containing peak flood magnitude 

# Magnitude of peak flow over which is considered a large flood event
highfloodcutoff = 700 # This is in CFS, as are the Maybell data points 
medfloodcutoff = 100 # non-event is 
# Same for drought (threshold below rather than above)
mindroughtlength = 30 # length of drought - a vector specified in csv. 

# Convert peak discharge values into a vector of floods/no floods, 
highflood_func <- function(x) {
    ifelse(x > highfloodcutoff, 1, 0)
}

medflood_func <- function(x) {
    ifelse(x > medfloodcutoff & x <= highfloodcutoff, 1, 0)
}

nonevent_func <- function(x,y) {
    ifelse(x <= medfloodcutoff & y <= mindroughtlength, 1, 0)
}

drought_func <- function(x,y) {
    ifelse(x <= medfloodcutoff & y > mindroughtlength, 1, 0)
}

# Defining year types here for now. May need to put into outer loop if simulating
highflood <- highflood_func(flowdata$floodmag)
medflood <- medflood_func(flowdata$floodmag)
drought <- drought_func(flowdata$floodmag, flowdata$basedur)
nonevent <- nonevent_func(flowdata$floodmag, flowdata$basedur)             
flood <- highflood + medflood # simply whether it's a flood year or not

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
highfloodoutput <- numeric(length = count) 
medfloodoutput <- numeric(length = count) 
droughtoutput <- numeric(length = count) 
noneventoutput <- numeric(length = count) 

# Fecundities
FCACLoutput <- numeric(length = count)
FGIROoutput <- numeric(length = count)
FLECYoutput <- numeric(length = count)



# VITAL RATES ------------------------------------------------

# Desert Sucker ----------------------------------------

# "Self thinning" rates, or equivalency rules, for stage transitions 
#bCACL1 <- denCACL2/denCACL1
#bCACL2 <- denCACL3/denCACL2

# Baseline maturation probabiliity, aCACL3 (adult senescence rate)
aCACL1 <- 1
aCACL2 <- 1 
aCACL3 <- .167

# Background mortality
S1MortCACL <- .8
S2MortCACL <- .5
S3MortCACL <- .1

# Initial volume in grams in 100-m reach
biomCACL1 <- 333
biomCACL2 <- 333
biomCACL3 <- 333

# Fecundity based on yeartype
# these are raw numbers taken from survey data: max, mean and min fecundities
maxfec.CACL <- 2772 
meanfec.CACL <- 1140
minfec.CACL <- 450

# Stage specific densities (ind./g)
denCACL3 <- .008 # real number based on avg. biomass of individuals on site
denCACL2 <- .077 # this will be calculated based on size of age 1
denCACL1 <- meanfec.CACL / ((380) * 0.05) # eggs; based on 5% of body weight going into ovaries (literature derived); 380 is avg. female weight (from bruder 2006)

# Chub ----------------------------------------

# "Self thinning" rates, or equivalency rules, for stage transitions 
#bGIRO1 <- denGIRO2/denGIRO1
#bGIRO2 <- denGIRO3/denGIRO2

# Baseline maturation probabiliity, aGIRO3 (adult senescence rate)
aGIRO1 <- 1
aGIRO2 <- 1 
aGIRO3 <- .125

# Background mortality
S1MortGIRO <- .8
S2MortGIRO <- .5
S3MortGIRO <- .1

# Initial volume in grams in 100-m reach
biomGIRO1 <- 333
biomGIRO2 <- 333
biomGIRO3 <- 333

# Fecundity based on yeartype
maxfec.GIRO <- 26903 
meanfec.GIRO <- 18699
minfec.GIRO <- 13816

# Stage specific densities (ind./g)
denGIRO3 <- .0069
denGIRO2 <- .05
denGIRO1 <- 1000  


# Green sunfish ----------------------------------------

# "Self thinning" rates, or equivalency rules, for stage transitions 
#bLECY1 <- denLECY2/denLECY1
#bLECY2 <- denLECY3/denLECY2

# Baseline maturation probabiliity, aLECY3 (adult senescence rate)
aLECY1 <- 1
aLECY2 <- 1 
aLECY3 <- .33

# Background mortality
S1MortLECY <- .8
S2MortLECY <- .5
S3MortLECY <- .1

# Initial weight in grams in 100-m reach 
biomLECY1 <- 333
biomLECY2 <- 333
biomLECY3 <- 333

# Fecundity based on yeartype
# THESE ARE MADE UP
maxfec.LECY <- 20000
meanfec.LECY <- 10000
minfec.LECY <- 5000

# Stage specific densities (ind./g)
denLECY3 <- .038
denLECY2 <- .5
denLECY1 <- 1000 

# N gives the total number of individuals for each age class.
# Initially here, this is found by multiplying the number of g occupied by a given class
# by the the density per g 

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

# Desert Sucker
#GCACL1 <- aCACL1 * bCACL1 * (1 - flood[y] * SpCACL1) * (1 - drought[y] * (SuCACL1 * droughtmult.N)) * (1 - nonevent[y] * (SuCACL1 * noneventmult.N))
#GCACL2 <- aCACL2 * bCACL2 * (1 - flood[y] * SpCACL2) * (1 - drought[y] * (SuCACL2 * droughtmult.N)) * (1 - nonevent[y] * (SuCACL2 * noneventmult.N))
#PCACL3 <- (1 - aCACL3) * (1 - flood[y] * SpCACL3) * (1 - drought[y] * (SuCACL3 * droughtmult.N)) * (1 - nonevent[y] * (SuCACL3 * noneventmult.N))

    GCACL1 <- aCACL1 *
        ((1 - (highflood[y] * S1MortCACL)) * CACL_Sp_HF * CACL_Su_HF) *
        ((1 - (medflood[y] * S1MortCACL)) * CACL_Sp_MF * CACL_Su_MF) *
        ((1 - (nonevent[y] * S1MortCACL)) * CACL_Sp_NE * CACL_Su_NE) *
        ((1 - (drought[y] * S1MortCACL)) * CACL_Sp_DR * CACL_Su_DR)
        
    GCACL2 <- aCACL2 *
        ((1 - (highflood[y] * S2MortCACL)) * CACL_Sp_HF * CACL_Su_HF) *
        ((1 - (medflood[y] * S2MortCACL)) * CACL_Sp_MF * CACL_Su_MF) *
        ((1 - (nonevent[y] * S2MortCACL)) * CACL_Sp_NE * CACL_Su_NE) *
        ((1 - (drought[y] * S2MortCACL)) * CACL_Sp_DR * CACL_Su_DR)

    PCACL3 <- (1 - aCACL3) *
        ((1 - (highflood[y] * S3MortCACL)) * CACL_Sp_HF * CACL_Su_HF) *
        ((1 - (medflood[y] * S3MortCACL)) * CACL_Sp_MF * CACL_Su_MF) *
        ((1 - (nonevent[y] * S3MortCACL)) * CACL_Sp_NE * CACL_Su_NE) *
        ((1 - (drought[y] * S3MortCACL)) * CACL_Sp_DR * CACL_Su_DR)
    
 

# Chub
#GGIRO1 <- aGIRO1 * bGIRO1 * (1 - flood[y] * SpGIRO1) * (1 - drought[y] * (SuGIRO1 * droughtmult.N)) * (1 - nonevent[y] * (SuGIRO1 * noneventmult.N))
#GGIRO2 <- aGIRO2 * bGIRO2 * (1 - flood[y] * SpGIRO2) * (1 - drought[y] * (SuGIRO2 * droughtmult.N)) * (1 - nonevent[y] * (SuGIRO2 * noneventmult.N))
#PGIRO3 <- (1 - aGIRO3) * (1 - flood[y] * SpGIRO3) * (1 - drought[y] * (SuGIRO3 * droughtmult.N)) * (1 - nonevent[y] * (SuGIRO3 * noneventmult.N))

     GGIRO1 <- aGIRO1 *
        ((1 - (highflood[y] * S1MortGIRO)) * GIRO_Sp_HF * GIRO_Su_HF) *
        ((1 - (medflood[y] * S1MortGIRO)) * GIRO_Sp_MF * GIRO_Su_MF) *
        ((1 - (nonevent[y] * S1MortGIRO)) * GIRO_Sp_NE * GIRO_Su_NE) *
        ((1 - (drought[y] * S1MortGIRO)) * GIRO_Sp_DR * GIRO_Su_DR)
        
    GGIRO2 <- aGIRO2 *
        ((1 - (highflood[y] * S2MortGIRO)) * GIRO_Sp_HF * GIRO_Su_HF) *
        ((1 - (medflood[y] * S2MortGIRO)) * GIRO_Sp_MF * GIRO_Su_MF) *
        ((1 - (nonevent[y] * S2MortGIRO)) * GIRO_Sp_NE * GIRO_Su_NE) *
        ((1 - (drought[y] * S2MortGIRO)) * GIRO_Sp_DR * GIRO_Su_DR)

    PGIRO3 <- (1 - aGIRO3) *
        ((1 - (highflood[y] * S3MortGIRO)) * GIRO_Sp_HF * GIRO_Su_HF) *
        ((1 - (medflood[y] * S3MortGIRO)) * GIRO_Sp_MF * GIRO_Su_MF) *
        ((1 - (nonevent[y] * S3MortGIRO)) * GIRO_Sp_NE * GIRO_Su_NE) *
        ((1 - (drought[y] * S3MortGIRO)) * GIRO_Sp_DR * GIRO_Su_DR)
    
 


# Green sunfish - note the "NN"
#GLECY1 <- aLECY1 * bLECY1 * (1 - flood[y] * SpLECY1) * (1 - drought[y] * (SuLECY1 * droughtmult.NN)) * (1 - nonevent[y] * SuLECY1)
#GLECY2 <- aLECY2 * bLECY2 * (1 - flood[y] * SpLECY2) * (1 - drought[y] * (SuLECY2 * droughtmult.NN)) * (1 - nonevent[y] * SuLECY2)
#PLECY3 <- (1 - aLECY3) * (1 - flood[y] * SpLECY3) * (1 - drought[y] * (SuLECY3 * droughtmult.NN)) * (1 - nonevent[y] * SuLECY3)

 GLECY1 <- aLECY1 *
        ((1 - (highflood[y] * S1MortLECY)) * LECY_Sp_HF * LECY_Su_HF) *
        ((1 - (medflood[y] * S1MortLECY)) * LECY_Sp_MF * LECY_Su_MF) *
        ((1 - (nonevent[y] * S1MortLECY)) * LECY_Sp_NE * LECY_Su_NE) *
        ((1 - (drought[y] * S1MortLECY)) * LECY_Sp_DR * LECY_Su_DR)
        
    GLECY2 <- aLECY2 *
        ((1 - (highflood[y] * S2MortLECY)) * LECY_Sp_HF * LECY_Su_HF) *
        ((1 - (medflood[y] * S2MortLECY)) * LECY_Sp_MF * LECY_Su_MF) *
        ((1 - (nonevent[y] * S2MortLECY)) * LECY_Sp_NE * LECY_Su_NE) *
        ((1 - (drought[y] * S2MortLECY)) * LECY_Sp_DR * LECY_Su_DR)

    PLECY3 <- (1 - aLECY3) *
        ((1 - (highflood[y] * S3MortLECY)) * LECY_Sp_HF * LECY_Su_HF) *
        ((1 - (medflood[y] * S3MortLECY)) * LECY_Sp_MF * LECY_Su_MF) *
        ((1 - (nonevent[y] * S3MortLECY)) * LECY_Sp_NE * LECY_Su_NE) *
        ((1 - (drought[y] * S3MortLECY)) * LECY_Sp_DR * LECY_Su_DR)
    
 
# Total grams occupied after year -----------------------------------------------
#totbiom.CACL <- 
#    biomCACL[1] * GCACL1/bCACL1 + 
#    biomCACL[2] * GCACL2/(bCACL2 * bCACL1) + 
#    biomCACL[3] * PCACL3/(bCACL2 * bCACL1) 

totbiom.CACL <- 
#    biomCACL[1] + 
    biomCACL[2] + 
    biomCACL[3] 

#totbiom.GIRO <- 
#    biomGIRO[1] * GGIRO1/bGIRO1 + 
#    biomGIRO[2] * GGIRO2/(bGIRO2 * bGIRO1) + 
#    biomGIRO[3] * PGIRO3/(bGIRO2 * bGIRO1) 

totbiom.GIRO <- 
#    biomGIRO[1] + 
    biomGIRO[2] + 
    biomGIRO[3] 

#totbiom.LECY <- 
#    biomLECY[1] * GLECY1/bLECY1 + 
#    biomLECY[2] * GLECY2/(bLECY2 * bLECY1) + 
#    biomLECY[3] * PLECY3/(bLECY2 * bLECY1) 

totbiom.LECY <- 
#    biomLECY[1] + 
    biomLECY[2] + 
    biomLECY[3] 



### :CLARIFY: at the moment carrying capacity is limiting spawning of all species based on the total biomass occupied at the end of the previous year. i.e. if above K, no spp spawn in that year. If spawning occurs, they can all spawn and there is no sequence, so overseeding COULD be massive.

# POTENTIAL CACL FECUNDITY ---------------------------------------------------------
    FCACL3 <- checkpos((adult_func(biomCACL[3] * denCACL3)) * # checks to see if at least 2 adults are present.
                       #(1/nonind(biomCACL[3])) *
    ifelse(highflood[y] == 1 & spawnwindow.CACL[y] == 1,
    (floor((biomCACL[3] * denCACL3) / 2) * maxfec.CACL) / denCACL1, # converting maxfec into g 
    ifelse(medflood[y] == 1 & spawnwindow.CACL[y] == 1,
           (floor((biomCACL[3] * denCACL3) / 2) * meanfec.CACL) / denCACL1,
           (floor((biomCACL[3] * denCACL3) / 2) * minfec.CACL) / denCACL1)) * # fecundity function of yeartype and whether flood occurred during spawning window
                   ifelse((totbiom.CACL + totbiom.GIRO + totbiom.LECY) < K, 1, 0)) # if K is already occupied, then no fecundity
                     
# '(1/nonind(NCACL[3]))' keeps FCACL3 from dividing by zero by substituting an arbitrary non-0 
# number that will be multiplied by 0 later anyway during matrix multiplication
# This means it's fecund per individual, rather than having overall fecundity double multipled by overall biomass (i.e. here and in matrix mult)     


# POTENTIAL GIRO FECUNDITY ---------------------------------------------------------
    FGIRO3 <- checkpos((adult_func(biomGIRO[3] * denGIRO3)) *  # checks to see if at least 2 adults are present. 
                       #(1/nonind(biomGIRO[3])) *
                   ifelse(highflood[y] == 1 & spawnwindow.GIRO[y] == 1,
                          (floor((biomGIRO[3] * denGIRO3) / 2) * maxfec.GIRO) / denGIRO1, # converting maxfec into g 
                   ifelse(medflood[y] == 1 & spawnwindow.GIRO[y] == 1,
                           (floor((biomGIRO[3] * denGIRO3) / 2) * meanfec.GIRO) / denGIRO1,
                           (floor((biomGIRO[3] * denGIRO3) / 2) * minfec.GIRO) / denGIRO1)) * # fecundity function of yeartype
                   ifelse((totbiom.CACL + totbiom.GIRO + totbiom.LECY) < K, 1, 0)) # if K is already occupied, then no fecundity

# POTENTIAL LECY FECUNDITY ---------------------------------------------------------
    FLECY3 <- checkpos((adult_func(biomLECY[3] * denLECY3)) * # checks to see if at least 2 adults are present.
                       #(1/nonind(biomLECY[3])) *
                   ifelse(highflood[y] == 1 & wipeoutwindow.LECY[y] == 1,
                          0, # converting maxfec into g
                   ifelse(highflood[y] == 1 & wipeoutwindow.LECY[y] == 0,
                           (floor((biomLECY[3] * denLECY3) / 2) * minfec.LECY) / denLECY1,
                   ifelse(medflood[y] == 1 & wipeoutwindow.LECY[y] == 1,
                          (floor((biomLECY[3] * denLECY3) / 2) * meanfec.LECY) / denLECY1,
                          (floor((biomLECY[3] * denLECY3) / 2) * meanfec.LECY) / denLECY1))) * # fecundity function of yeartype
                   ifelse((totbiom.CACL + totbiom.GIRO + totbiom.LECY) < K, 1, 0)) # if K is already occupied, then no fecundity






# K --------------------------------------------------------------------------------------
# CACL
# gives total CACL biomass as a percentage of K; 
#KCACL <- 100 * (biomCACL[1] + 
#             biomCACL[2]/(bCACL1) + 
#             biomCACL[3]/(bCACL2 * bCACL1))/K

KCACL <- 100 * (biomCACL[1] + 
             biomCACL[2] + 
             biomCACL[3])/K

# GIRO
# gives total GIRO biomass as a percentage of K; 
#KGIRO <- 100 * (biomGIRO[1] + 
#             biomGIRO[2]/(bGIRO1) + 
#             biomGIRO[3]/(bGIRO2 * bGIRO1))/K

KGIRO <- 100 * (biomGIRO[1] + 
             biomGIRO[2] + 
             biomGIRO[3])/K

# LECY
# gives total LECY biomass as a percentage of K; 
#KLECY <- 100 * (biomLECY[1] + 
#             biomLECY[2]/(bLECY1) + 
#             biomLECY[3]/(bLECY2 * bLECY1))/K 

KLECY <- 100 * (biomLECY[1] + 
             biomLECY[2] + 
             biomLECY[3])/K
    
# TRANSITION MATRICES --------------------------------------------------------------------

# TRANSITION MATRIX FOR CACL -------------------------------------------------------
ACACL1 <- c(0, 0, 0)
ACACL2 <- c(GCACL1, 0, 0)
ACACL3 <- c(0, GCACL2, PCACL3)
# Matrix
ACACL <- rbind(ACACL1, ACACL2, ACACL3)

# TRANSITION MATRIX FOR GIRO -------------------------------------------------------
AGIRO1 <- c(0, 0, 0)
AGIRO2 <- c(GGIRO1, 0, 0)
AGIRO3 <- c(0, GGIRO2, PGIRO3)
# Matrix
AGIRO <- rbind(AGIRO1, AGIRO2, AGIRO3)

# TRANSITION MATRIX FOR LECY -------------------------------------------------------
ALECY1 <- c(0, 0, 0)
ALECY2 <- c(GLECY1, 0, 0)
ALECY3 <- c(0, GLECY2, PLECY3)
# Matrix
ALECY <- rbind(ALECY1, ALECY2, ALECY3)



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
highfloodoutput[i] <- highflood[y]
medfloodoutput[i] <- medflood[y]
noneventoutput[i] <- nonevent[y]
droughtoutput[i] <- drought[y]


CACLplaceholder <- FCACL3
GIROplaceholder <- FGIRO3
LECYplaceholder <- FLECY3
    
# MATRIX MULTIPLICATION ------------------------------------------------------------------
# CACL
biomCACL <- ACACL %*% biomCACL # ACACL is transition matrix, biomCACL = total biomass for each age class
    biomCACL[1] <- CACLplaceholder
    
#nCACLpost <- c(biomCACL[1] * denCACL1, biomCACL[2] * denCACL2, biomCACL[3] * denCACL3)
#nCACLpost <- ifelse(nCACLpost < 1, 0, 1)
#biomCACL <- biomCACL * nCACLpost

# GIRO
biomGIRO <- AGIRO %*% biomGIRO # AGIRO is transition matrix, biomGIRO = total biomass for each age class
    biomGIRO[1] <- GIROplaceholder
    
#nGIROpost <- c(biomGIRO[1] * denGIRO1, biomGIRO[2] * denGIRO2, biomGIRO[3] * denGIRO3)
#nGIROpost <- ifelse(nGIROpost < 1, 0, 1)
#biomGIRO <- biomGIRO * nGIROpost

# LECY
biomLECY <- ALECY %*% biomLECY # ALECY is transition matrix, biomLECY = total biomass for each age class
    biomLECY[1] <- LECYplaceholder
    
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


ggplot(flowdata, aes(year, floodmag)) +
    geom_point() + geom_path()
ggplot(flowdata, aes(year, basedur)) +
    geom_point() + geom_path()

flowresults <- as.data.frame(cbind(highfloodoutput, medfloodoutput, noneventoutput, droughtoutput)) %>%
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
