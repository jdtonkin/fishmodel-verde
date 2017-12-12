## Functions used in 'fish-model-all-spp.R'

## checkpos makes sure that the K-occupied term is positive, assigns 0 if not
checkpos <- function(x) {
    ifelse(x < 0, 0, x)
}

## Flow threshold functions ----------------------------------------------------

## Spring
SP_highflood_func <- function(x) {
    ifelse(x >= SP_highfloodcutoff, 1, 0)
}

## Summer
SU_highflood_func <- function(x) {
    ifelse(x >= SU_highfloodcutoff, 1, 0)
}

## Medium flood
medflood_func <- function(x) {
    ifelse(x >= medfloodcutoff &
           x < SP_highfloodcutoff, 1, 0)
}

## Nonevent
nonevent_func <- function(Spfl, BD, Sufl) {
    ifelse(Spfl < medfloodcutoff &
           BD < mindroughtlength &
           Sufl < SU_highfloodcutoff, 1, 0)
}

## Drought
drought_func <- function(Spfl, BD, Sufl) {
    ifelse(Spfl < medfloodcutoff &
           BD >= mindroughtlength &
           Sufl < SU_highfloodcutoff, 1, 0)
}

## Plotting theme
theme_classic_facet <- function() {
    theme_classic() +
        theme(strip.background = element_rect(colour = NA, fill = NA))
}



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
