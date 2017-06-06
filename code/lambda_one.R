# Jane Rogosch
# 5 June 2017
# Verde Fish model

?floor


# Vital rates sans flow -----------------------------------------------------
# CACL
biom.CACL3 <- 333
den.CACL3 <- 0.019
meanfec.CACL <- 1140
den.CACL1 <- 1000
a.CACL3 <- 0.08

F.CACL3 <- floor((biom.CACL3 * den.CACL3) / 2) * meanfec.CACL / den.CACL1
S.CACL1 <- 0.10
S.CACL2 <- 0.8
S.CACL3 <- 0.8
P.CACL3 <- (1 - a.CACL3)*S.CACL3

# GIRO
biom.GIRO3 <- 333
den.GIRO3 <- 0.0024
meanfec.GIRO <- 16324
den.GIRO1 <- 1313
a.GIRO3 <- 0.125

F.GIRO3 <- floor((biom.GIRO3 * den.GIRO3) / 2 * meanfec.GIRO / den.GIRO1)
S.GIRO1 <- 0.09
S.GIRO2 <- 0.81 
S.GIRO3 <- 0.81
P.GIRO3 <- (1 - a.GIRO3)*S.GIRO3  

# LECY
biom.LECY3 <- 333
den.LECY3 <- 0.015
meanfec.LECY <- 1000
den.LECY1 <- 513
a.LECY3 <- 0.14

F.LECY3 <- floor((biom.LECY3 * den.LECY3)) / 2 * meanfec.LECY / den.LECY1
S.LECY1 <- 0.24
S.LECY2 <- 0.56
S.LECY3 <- 0.56
P.LECY3 <- (1 - a.LECY3)*S.LECY3  
  
# TRANSITION MATRIX FOR CACL -------------------------------------------------------
A.CACL1 <- c(0, 0, F.CACL3)
A.CACL2 <- c(S.CACL1, 0, 0)
A.CACL3 <- c(0, S.CACL2, P.CACL3)
# Matrix
A.CACL <- rbind(A.CACL1, A.CACL2, A.CACL3)
lambda(A.CACL) # Checking population growth rate

gCACL1vec <- c(0.05, 0.06, 0.08, 0.10, 0.09, 0.095)
lambda.gCACL.vec <- c(0.9035605, 0.9270232, 0.9690734, 1.006225, 0.9881707, 0.9973192) # taken from result of corresponding value to S.CACL1

# TRANSITION MATRIX FOR GIRO -------------------------------------------------------
A.GIRO1 <- c(0, 0, F.GIRO3)
A.GIRO2 <- c(S.GIRO1, 0, 0)
A.GIRO3 <- c(0, S.GIRO2, P.GIRO3)
# Matrix
A.GIRO <- rbind(A.GIRO1, A.GIRO2, A.GIRO3)
lambda(A.GIRO) # Checking population growth rate

gGIRO1vec <- c(0.05, 0.10, 0.09)
lambda.gGIRO.vec <- c(0.9060765, 1.020105, 1.000221)
# TRANSITION MATRIX FOR LECY -------------------------------------------------------
A.LECY1 <- c(0, 0, F.LECY3)
A.LECY2 <- c(S.LECY1, 0, 0)
A.LECY3 <- c(0, S.LECY2, P.LECY3)
# Matrix
A.LECY <- rbind(A.LECY1, A.LECY2, A.LECY3)
lambda(A.LECY) # Checking population growth rate

gLECY1vec <- c(0.05, 0.10, 0.2, 0.3, 0.25, 0.24)
lambda.gLECY.vec <- c(0.7026819, 0.812398, 0.9576854, 1.062157, 1.013239, 1.002729)
