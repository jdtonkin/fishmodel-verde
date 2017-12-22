## Verde Fish Model
## Jane Rogosch, Jono Tonkin, et al.

## Required libraries
library(ggplot2)
library(plyr)
library(tidyverse)
library(popbio)

## * SETUP ---------------------------------------------------------------------

rm(list = ls()) # clearing the workspace 

## Bringing in flow data
all.scenarios.list <- readRDS('data/all_scenarios_list.rds')

## Making a drought and nonevent scenario list
## Includes natural too

## names of dataframes to match partially
toMatch <- c('drought', 'nonevent', 'natural')

## selecting dataframes in the list with partial strings matching those in toMatch
dr.ne.scenarios.list <- all.scenarios.list[grep(paste(toMatch, collapse = "|"), names(all.scenarios.list), value = TRUE)]


## When using natural flow data, just pull it out from the list here
flowdata <- all.scenarios.list$natural.flow

count <- 10 # how many years to project ## Prob go for 250
iterations <- 1 # number of replicate projections to run (mid loop) ## Prob go for 500

## Modifiers
modifiers <- read.csv('data/modifiers-all-spp.csv')

## adding 'Modifier' value from csv to 'Code' in csv
for(j in 1:length(modifiers[,1])) {
    nam <- paste(modifiers[j,4])
    assign(nam, modifiers[j,5]) 
}

## Vital rates
## Baseline maturation probability, aCACL3 (adult senescence rate)
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

## * Key -----------------------------------------------------------------------
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

## Create a list of all possible combinations of species
## e.g. combn(sppnames, 4, simplify = T)

## Make list of all combose - a nested list
allcombos <- lapply(seq_along(sppnames), function(y)
    combn(sppnames, y, simplify = F))

## Renaming the outer list so these are kept when unlisting
names(allcombos) <- c('one_', 'two_', 'three_', 'four_', 'five_', 'six_', 'seven')

## Unlisting to bring up sub-lists into a single list.
## recursive = FALSE is key here to just do the one level of unlisting.
## Names are now combined. e.g. one_1 
allcombos <- unlist(allcombos, recursive = FALSE)
allcombos[[29]]
allcombos$three_1
allcombos[['three_1']]
names(allcombos)




### ----------------------------------------------------------------------------
## outer spp selection loop goes here
### ----------------------------------------------------------------------------

## set this to 'sppnames' if you want to run on all seven 
## make something here that has all 7 spp but a 0 or 1 next to it based on
## whether it was included in list or not. 
sppvec <- allcombos[['six_1']]

## This has two columns, the sppnames and the value to multiply biomass with
spp.val <- as.data.frame(sppnames) %>%
    mutate(val = ifelse(match(sppnames, sppvec, nomatch = 0) > 0, 1, 0))

rownames(spp.val) <- spp.val$sppnames


### ----------------------------------------------------------------------------
## outer flow modification loop goes here
### ----------------------------------------------------------------------------


K = 47660 # mean for 1-km reach across 6 replicate reaches

## Loading functions from functions.R file -------------------------------------
source('code/functions.R')

## * ITERATION PARAMETERS ------------------------------------------------------
## Setting up arrays/vectors to fill with data from loops

## Mid loop details ------------------------------------------------------------
## 'iterations' - number of replicate flow sequences to run for averaging

stages <- as.character(c("S1", "S2", "S3"))

## Total. N of stages 2 and 3 each year ----------------------------------------
## Takes all stages 2 and 3 and sums them for each year and iteration
Total.N <- array(0,
                 dim = c(count, iterations),
                 dimnames = list(1:count, 1:iterations)
                 )

## replist. List of arrays w/ abundance data for each spp ----------------------
## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
reparray <- array(0,
                  dim = c(count, 3, iterations),
                  dimnames = list(1:count, stages, 1:iterations)
                  )

## Repeating the array 7 times 
replist <- rep(list(reparray), 7)

## Assigning names to each array from sppnames vector
names(replist) <- sppnames

## Inner loop details ----------------------------------------------------------
## 'count' - number of years to project simulations (inner loop)

## N ---------------------------------------------------------------------------
## Output of biomass and no. ind. for each age class for each year projected  
## Array w/ 3 cols (stage classes) and however many rows there are yrs projected
## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
output.N.array <- array(0, dim = c(count, 3))

## Repeating the array 7 times 
output.N.list <- rep(list(output.N.array), 7)

## Assigning names to each array from sppnames vector
names(output.N.list) <- sppnames

## Biomass ---------------------------------------------------------------------
## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
output.biom.array <- array(0, dim = c(count, 3))

## Repeating the array 7 times 
output.biom.list <- rep(list(output.biom.array), 7)

## Assigning names to each array from sppnames vector
names(output.biom.list) <- sppnames

## Total biomass as % of K -----------------------------------------------------
## Creating a list of 7 vectors to fill in. One for each spp. 
## Create a vector to be repeated
biomoutput.vector <- numeric(length = count)

## Repeating the vector 7 times 
biomoutput.list <- rep(list(biomoutput.vector), 7)

## Assigning names to each vector from sppnames vector
names(biomoutput.list) <- sppnames

## Flow results ----------------------------------------------------------------
## Flood and drought settings for each yr projected into future (i.e. 0 or 1) 
## Create data frame with 5 cols and 'count' rows to fill in with flow results
flowresults <- data.frame(matrix(ncol = 5, nrow = count))
names(flowresults) <- c('SPhighflood',
                        'SUhighflood',
                        'medflood',
                        'drought',
                        'nonevent')
flowresults

## Fecundities -----------------------------------------------------------------
## Creating a list of 7 vectors to fill in. One for each spp. 
## Create a vector to be repeated
fec.vector <- numeric(length = count)

## Repeating the vector 7 times 
fec.list <- rep(list(fec.vector), 7)

## Assigning names to each vector from sppnames vector
names(fec.list) <- sppnames

### ----------------------------------------------------------------------------
## * Mid loop ##################################################################
### ----------------------------------------------------------------------------
## Middle loop uses iterator "iter" to get "iterations" for suming S2 and S3
for(iter in 1:iterations) {

    ## USE THIS to examine different flow year types +++++++++++++++++++++++++++ 
    ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
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
    
    ## Need to read in initial biom every time so starting biomass is reset each
    ## iteration
    ## N gives the total number of individuals for each age class.
    ## Initially here, this is found by multiplying the number of g occupied by
    ## a given class by the density per g
    ## biom = g/m3
    ## den = indiv/g
    
    ## To have different initial starting population sizes for each iteration,
    ## taking biom of stage 3 from negative binomial distribution, where the
    ## parameter (lambda = mean) and K (dispersion) is calculated from mean and
    ## variance in abundance across seven sites in Verde River from 94-08, and
    ## scaled to biomass from Gibson 2012 survey in file:
    ## "Rinne Verde River Data 1994-2008-.xlsx"
    
    biomCACL <- c(biomCACL1 * spp.val['CACL',]$val,
                  biomCACL2 * spp.val['CACL',]$val,
                  rnbinom(1, size = 1.52, mu = 5284) * spp.val['CACL',]$val) 
    
    biomGIRO <- c(biomGIRO1 * spp.val['GIRO',]$val,
                  biomGIRO2 * spp.val['GIRO',]$val,
                  rnbinom(1, 0.44, mu = 2376) * spp.val['GIRO',]$val) 
    
    biomLECY <- c(biomLECY1 * spp.val['LECY',]$val,
                  biomLECY2 * spp.val['LECY',]$val,
                  rnbinom(1, 0.34, mu = 164) * spp.val['LECY',]$val) 
    
    biomCAIN <- c(biomCAIN1 * spp.val['CAIN',]$val,
                  biomCAIN2 * spp.val['CAIN',]$val,
                  rnbinom(1, 1.33, mu = 34068) * spp.val['CAIN',]$val) 
    
    biomMIDO <- c(biomMIDO1 * spp.val['MIDO',]$val,
                  biomMIDO2 * spp.val['MIDO',]$val,
                  rnbinom(1, 0.66, mu = 4202) * spp.val['MIDO',]$val) 
    
    biomCYLU <- c(biomCYLU1 * spp.val['CYLU',]$val,
                  biomCYLU2 * spp.val['CYLU',]$val,
                  rnbinom(1, 1.78, mu = 238) * spp.val['CYLU',]$val) 
    
    biomAMNA <- c(biomAMNA1 * spp.val['AMNA',]$val,
                  biomAMNA2 * spp.val['AMNA',]$val,
                  rnbinom(1, 0.36, mu = 1306) * spp.val['AMNA',]$val)

### ----------------------------------------------------------------------------
## * Inner loop ################################################################
### ----------------------------------------------------------------------------
    for(i in 1:count) {

        ## CHANGE WHAT 'y' IS TO SIMULATE DIFFERENT FLOW REGIMES ACROSS THE 54 Y 
        ## y is directly taken from flow vector. 
        #y = i # follow flow record sequence

        ## Sampling randomly from the flow record
        y = sample(nrow(flowdata), 1)


        ## Transition probabilities  -------------------------------------------
        ## G is prob. of transition to next stage
        ## P is prob. of remaining in that stage
        ## Baseline mortality vital rate (from file object: 'vital-rates.csv') 
        ## 'STmort...' is multiplied by modifier (from file object:
        ## 'modifiers-all-species.csv') based on yeartype as specified above
        ## So we have SP_highflood, SU_highflood, medflood, nonevent, drought
        ## and '..._J/2/3_SUHF/SPHF/MF/NE/DR'

        ## Stage 1 - G
        for(nm in sppnames) {
            assign(paste0('G', nm, '1'),
                   
                   get(paste0('a', nm, '1')) *
                   get(paste0('den', nm, 'J')) *
                   (1 /
                    get(paste0('den', nm, '2'))
                   )  *
                   (1 -
                    (flowdata$SU_highflood[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_J_SUHF')))
                   ) *
                   (1 -
                    (flowdata$SP_highflood[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_J_SPHF')))
                   ) *
                   (1 -
                    (flowdata$medflood[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_J_MF')))
                   ) *
                   (1 -
                    (flowdata$nonevent[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_J_NE')))
                   ) *
                   (1 -
                    (flowdata$drought[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_J_DR')))
                   )
                   )          
        }        

        ## Stage 2 - G
        for(nm in sppnames) {
            assign(paste0('G', nm, '2'),

                   get(paste0('a', nm, '2')) *
                   get(paste0('den', nm, '2')) *
                   (1 /
                    get(paste0('den', nm, '3'))
                   ) *
                   (1 -
                    (flowdata$SU_highflood[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_A_SUHF')))
                   ) *
                   (1 -
                    (flowdata$SP_highflood[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_A_SPHF')))
                   ) *
                   (1 -
                    (flowdata$medflood[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_A_MF')))
                   ) *
                   (1 -
                    (flowdata$nonevent[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_A_NE')))
                   ) *
                   (1 -
                    (flowdata$drought[y] *
                     get(paste0('STMort', nm)) *
                     get(paste0(nm, '_A_DR')))
                   )
                   )          
        }        

        ## Stage 3 - P
        for(nm in sppnames) {
            assign(paste0('P', nm, '3'),

            (1 -
             get(paste0('a', nm, '3'))
            ) *
            (1 -
             (flowdata$SU_highflood[y] *
              get(paste0('STMort', nm)) *
              get(paste0(nm, '_A_SUHF')))
            ) *
            (1 -
             (flowdata$SP_highflood[y] *
              get(paste0('STMort', nm)) *
              get(paste0(nm, '_A_SPHF')))
            ) *
            (1 -
             (flowdata$medflood[y] *
              get(paste0('STMort', nm)) *
              get(paste0(nm, '_A_MF')))
            ) *
            (1 -
             (flowdata$nonevent[y] *
              get(paste0('STMort', nm)) *
              get(paste0(nm, '_A_NE')))
            ) *
            (1 -
             (flowdata$drought[y] *
              get(paste0('STMort', nm)) *
              get(paste0(nm, '_A_DR')))
            )
            )          
        }          
        
        ## POTENTIAL FECUNDITY -------------------------------------------------

        ## 1st calculate total grams occupied after year 
        totbiom <- ldply(sppnames, function(x)
            get(paste0('biom', x))[1] +
            get(paste0('biom', x))[2] +
            get(paste0('biom', x))[3]) %>%
            sum
        
        ## Carrying capacity (K) is limiting spawning of all species based on
        ## the total biomass occupied at the end of the previous year. i.e. if
        ## above K, no spp spawn in that year. If spawning occurs, they all do.
        ## Some slight differences in fecund btwn spp so can't loop/lapply

        ## CACL stage 2
        FCACL2 <- ((0.5 * GSI.CACL * (1 - S0MortCACL)) *
                   checkpos((K - totbiom)/K)) *
            denCACL1 *
            (1/denCACLJ)

        ## CACL stage 3
        FCACL3 <- ((0.5 * GSI.CACL * (1 - S0MortCACL)) *
                   checkpos((K - totbiom)/K)) *
            denCACL1 *
            (1/denCACLJ)

        ## GIRO stage 2
        FGIRO2 <- ((0.5 * GSI.GIRO * (1 - S0MortGIRO)) *
                   checkpos((K - totbiom)/K)) * 
            denGIRO1 *
            (1/denGIROJ)

        ## GIRO stage 3
        FGIRO3 <- ((0.5 * GSI.GIRO * (1 - S0MortGIRO)) *
                   checkpos((K - totbiom)/K)) *
            denGIRO1 *
            (1/denGIROJ)

        ## CAIN stage 3
        FCAIN3 <- ((0.5 * GSI.CAIN * (1 - S0MortCAIN)) *
                   checkpos((K - totbiom)/K)) *
            denCAIN1 *
            (1/denCAINJ)

        ## LECY stage 2
        FLECY2 <- ((0.5 * GSI.LECY * (1 - S0MortLECY)) *
                   checkpos((K - totbiom)/K)) *
            denLECY1 *
            (1/denLECYJ)

        ## LECY stage 3
        FLECY3 <- ((0.5 * GSI.LECY * (1 - S0MortLECY)) *
                   checkpos((K - totbiom)/K)) * 
            denLECY1 *
            (1/denLECYJ)
        
        ## MIDO stage 3
        FMIDO3 <- ((0.5 * GSI.MIDO * (1 - S0MortMIDO)) *
                   checkpos((K - totbiom)/K)) *
            denMIDO1 *
            (1/denMIDOJ)

        ## CYLU 
        ## because they are serial spawners, they are allowed to spawn twice a
        ## season in stage 2 and 3
        ## CYLU stage 1
        FCYLUJ <- ((0.5 * GSI.CYLU * (1 - S0MortCYLU)) *
                   checkpos((K - totbiom)/K)) * 
            denCYLU1 *
            (1/denCYLUJ)

        ## CYLU stage 2
        FCYLU2 <- ((0.5 * 2 * GSI.CYLU * (1 - S0MortCYLU)) *
                   checkpos((K - totbiom)/K)) *
            denCYLU1 *
            (1/denCYLUJ)

        ## CYLU stage 3
        FCYLU3 <- ((0.5 * 2 * GSI.CYLU * (1 - S0MortCYLU)) *
                   checkpos((K - totbiom)/K)) *
            denCYLU1 *
            (1/denCYLUJ)

        ## AMNA 
        FAMNA3 <- ((0.5 * GSI.AMNA * (1 - S0MortAMNA)) *
                   checkpos((K - totbiom)/K)) *
            denAMNA1 *
            (1/denAMNAJ)

        ## K -------------------------------------------------------------------
        ## Calculating the percentage of K occupied and adding to KCACL etc
        for(nm in sppnames) {
            assign(paste0('K', nm), 100 * (
                get(paste0('biom', nm))[1] +   
                get(paste0('biom', nm))[2] +
                get(paste0('biom', nm))[3])/K)                                  
        }
        
        ## TRANSITION MATRICES -------------------------------------------------

        ## CACL 
        ACACL1 <- c(0, FCACL2, FCACL3)
        ACACL2 <- c(GCACL1, 0, 0)
        ACACL3 <- c(0, GCACL2, PCACL3)
        
        ## GIRO 
        AGIRO1 <- c(0, FGIRO2, FGIRO3)
        AGIRO2 <- c(GGIRO1, 0, 0)
        AGIRO3 <- c(0, GGIRO2, PGIRO3)
       
        ## LECY 
        ALECY1 <- c(0, FLECY2, FLECY3)
        ALECY2 <- c(GLECY1, 0, 0)
        ALECY3 <- c(0, GLECY2, PLECY3)
       
        ## CAIN 
        ACAIN1 <- c(0, 0, FCAIN3)
        ACAIN2 <- c(GCAIN1, 0, 0)
        ACAIN3 <- c(0, GCAIN2, PCAIN3)
        
        ## MIDO 
        AMIDO1 <- c(0, 0, FMIDO3)
        AMIDO2 <- c(GMIDO1, 0, 0)
        AMIDO3 <- c(0, GMIDO2, PMIDO3)
        
        ## CYLU 
        ACYLU1 <- c(FCYLUJ, FCYLU2, FCYLU3)
        ACYLU2 <- c(GCYLU1, 0, 0)
        ACYLU3 <- c(0, GCYLU2, PCYLU3)
        
        ## AMNA 
        AAMNA1 <- c(0, 0, FAMNA3)
        AAMNA2 <- c(GAMNA1, 0, 0)
        AAMNA3 <- c(0, GAMNA2, PAMNA3)
       
        ## rbinding the vectors from above into transition matrices
        ## Makes ACACL, AGIRO etc.
        for(nm in sppnames) {
            assign(paste0('A', nm), rbind(
                                        get(paste0('A', nm, '1')),
                                        get(paste0('A', nm, '2')),
                                        get(paste0('A', nm, '3'))
                                    ))
        }
        
        ## COMPILING OUTPUTS ---------------------------------------------------

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
            output.N.list[[nm]][i,1:3] <- c(get(paste0('biom', nm))[1] *
                                            get(paste0('den', nm, 'J')
                                                ),
                                            get(paste0('biom', nm))[2] *
                                            get(paste0('den', nm, '2')
                                                ),
                                            get(paste0('biom', nm))[3] *
                                            get(paste0('den', nm, '3')
                                                ))
        }

        ## Flow results --------------------------------------------------------
        ## Records flood settings of each particular projected year
        ## (1 = yes, 0 = no)
        flowresults$SPhighflood[i] <- flowdata$SP_highflood[y]
        flowresults$SUhighflood[i] <- flowdata$SU_highflood[y]
        flowresults$medflood[i] <- flowdata$medflood[y]
        flowresults$nonevent[i] <- flowdata$nonevent[y]
        flowresults$drought[i] <- flowdata$drought[y]
        
        ## MATRIX MULTIPLICATION -----------------------------------------------
        ## can include rescue function for each with 0.5 chance of reach being
        ## colonized by 2 individuals
        ## Loop essentially == biomAMNA <- AAMNA %*% biomAMNA
        ## AAMNA is transit. matrix, biomAMNA = total biomass for each age class
        for(nm in sppnames) {
            assign(paste0('biom', nm),
                   get(paste0('A', nm)) %*% get(paste0('biom', nm)))
        }
               
    }
### ----------------------------------------------------------------------------
### End of inner loop ##########################################################
### ----------------------------------------------------------------------------

    ## Mean values for each iteration run over each sequence of years ----------
    for(nm in sppnames) {
        replist[[nm]][,,iter] <- output.N.list[[nm]]    
    }
    
    ## Total.N -----------------------------------------------------------------
    ## Caculating Total.N for each year, and adding it to total.N data frame
    ## with however many iterations run.
    ## Total does not incl. juveniles.
    ## map is purrr version of lapply. Can pass fn using ~ and .x instead of
    ## function(x) x
    ## Gets list output of stages 2:3 for ea spp, then cbinds them all together,
    ## then calcs sum. 
    Total.N[,iter] <- map(output.N.list, ~ .x[,2:3]) %>%
        do.call('cbind', .) %>%
        apply(1, sum)

}
### ----------------------------------------------------------------------------
### End of mid loop ############################################################
### ----------------------------------------------------------------------------

## Saving image here as been having trouble in next steps
save.image()

## * OUTPUTS -------------------------------------------------------------------
################################################################################


## turning replist into a df
repdf <- ldply(replist, function(x) {
    adply(x, c(1,2,3))
})

head(repdf, 70)
names(repdf) <- c('species', 'year', 'stage', 'rep', 'abund')
repdf <- filter(repdf, stage != 'S1')
repdf$year <- as.numeric(as.character(repdf$year))

totn <- adply(Total.N, c(1,2))
names(totn) <- c('year', 'rep', 'tot.abund')
totn$year <- as.numeric(as.character(totn$year))

## joining totn and repdf together
repdf <- left_join(totn, repdf)

## calculating relative abundance
repdf <- mutate(repdf, rel.abund = abund/tot.abund)

## Taking mean results to cf w/ observed data
means <- repdf %>%
    select(-tot.abund) %>%
    group_by(year, rep, species) %>% # combining stages
    summarise(abund = sum(abund),
              rel.abund = sum(rel.abund)) %>%
    ungroup() %>%
    group_by(species, year) %>%
    summarise(mean.abund = mean(abund),
              sd.abund = sd(abund),
              se.abund = sd(abund)/sqrt(iterations),
              mean.rel.abund = mean(rel.abund),
              sd.rel.abund = sd(rel.abund),
              se.rel.abund = sd(rel.abund)/sqrt(iterations)) %>%
    ungroup()
means

## Seeing if there are 0s
## means[means$mean.rel.abund == 0,]

## Plotting
rel.abund.trends <- ggplot(means, aes(year,
                                      mean.rel.abund,
                                      colour = species,
                                      fill = species)) +
    geom_ribbon(aes(ymin = mean.rel.abund - 2 * se.rel.abund,
                    ymax = mean.rel.abund + 2 * se.rel.abund),
                colour = 'transparent',
                alpha = .5,
                show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    facet_wrap(~species, ncol = 2) +
    theme_classic_facet() +
    coord_cartesian(ylim = c(0,1)) +
    ylab('Relative abundance') +
    xlab('Year')
## adding observed data
rel.abund.trends 
ggsave('export/multi-spp-all-year.pdf', width = 6, height = 6)




### Local Variables:
### eval: (orgstruct-mode 1)
### orgstruct-heading-prefix-regexp: "## "
### End:
