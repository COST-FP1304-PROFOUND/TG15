library(rjags)
library(BayesianTools)
set.seed(123)
ndays                 <- 1095
PAR                   <- VSEMcreatePAR(1:ndays)

refPars               <- VSEMgetDefaults()
referenceData         <- VSEM(refPars$best, PAR)
obs                   <- referenceData + rnorm(length(referenceData), sd = (abs(referenceData) + 1E-7) * 0.1)
obs[,1]               <- referenceData[,1] + rnorm(length(referenceData[,1]), sd = 0.003)

NEEobs <- obs[,1]
Cvobs  <- obs[,2]
Csobs  <- obs[,3]
days   <- seq(ndays)

VSEM <- " model {
  GAMMA  = 0.4
  LAR    = 1.5
  Cr[1]  = 3.0
  tauR   = 1440.
  Av     = 0.5
  Cv[1]  ~ dunif(  1., 10. )
  Cs[1]  ~ dunif(  1., 20. )
  ## update the states
  for (i in 2:n) {
    Cv[i]  <- Cv[i-1]  + Av       * NPP[i-1] - Cv[i-1]/(tauV*1000.)
    Cr[i]  <- Cr[i-1]  + (1.0-Av) * NPP[i-1] - Cr[i-1]/tauR
    Cs[i]  <- Cs[i-1]  + Cr[i-1]/tauR + Cv[i-1]/(tauV*1000.) - Cs[i-1]/(tauS*10000.)
  }
  for (i in 1:n) {
      ## Fluxes
      G[i]       <- PAR[i] * (LUE/1000.) * (1.0 - exp(-1.0*KEXT*LAR*Cv[i]))
      NPP[i]     <- (1.0-GAMMA)*G[i]
      NEE[i]     <- Cs[i]/(tauS*10000.) - NPP[i]
      ## NEE[i]  <- (Cs[i]/(tauS*10000.) + GAMMA*G[i]) - G[i]

      ## Calculate precisions
      tauNEE[i] <- pow( 0.003 , -2)
      tauCv[i]  <- pow( cvCv   * Cv[i]  , -2 )
      tauCs[i]  <- pow( cvCs   * Cs[i]  , -2 )

      ## Likelihood
      NEEobs[i]    ~ dnorm( NEE[i], tauNEE[i] )
      Cvobs[i]     ~ dnorm( Cv[i], tauCv[i] )
      Csobs[i]     ~ dnorm( Cs[i], tauCs[i] )
  }

  KEXT   ~ dunif( 0.3, 1.0 )
  LUE    ~ dunif( 1.0, 5.0 )
  tauV   ~ dunif( 1.0, 3.0 )
  tauS   ~ dunif( 1.0, 4.0 )
  cvCv   ~ dunif( 0.0, 1.0 )
  cvCs   ~ dunif( 0.0, 1.0 )
}"

truth = refPars$best
names(truth) = row.names(refPars)
ic       = as.list(truth)
ic$LAR   = NULL
ic$GAMMA = NULL
ic$tauR  = NULL
ic$Av    = NULL
ic[5:8]  = NULL
ic$LUE   = ic$LUE*1000
ic$tauV  = ic$tauV/1000
ic$tauS  = ic$tauS/10000

VSEMdata       <- list(n=ndays, NEEobs=NEEobs, Cvobs=Cvobs, Csobs=Csobs, PAR=PAR)
#VSEMic         <- list(Cv.ic=Cvobs[1],Cs.ic=Csobs[1])
VSEMoutputs    <- c( "cvCv", "cvCs", "KEXT", "LUE",
                    "tauV","tauS","NEE","G","PAR",
                     "Cv","Cs","Cr")
## VSEMoutputs    <- c( "cvCv", "cvCs", "KEXT", "LUE",
##                     "tauV","tauS")
nadapt         <- 1000; nbi <- 1000; nit <- 3000
jagsVSEM       <- jags.model ( textConnection(VSEM), data=VSEMdata, n.chains=3, n.adapt=nadapt, inits = ic)
update( jagsVSEM, n.iter=nbi )
codaSamples    <- coda.samples( jagsVSEM, var=VSEMoutputs, n.iter=nit, thin= 9 )
gelman.diag(codaSamples)
mcmcChain      <- as.matrix( codaSamples )

##
## Parameter histograms
##

par(mfrow=c(3,3))
hist( mcmcChain[,"KEXT"], xlab="", ylab="",
      main=paste( "KEXT \n( sd =",signif(sd(mcmcChain[,"KEXT"]),2), ")" )
     )
abline(v=truth["KEXT"],col=2)

hist( mcmcChain[,"LUE"]/1000, xlab="", ylab="",
      main=paste( "LUE \n( sd =",signif(sd(mcmcChain[,"LUE"]),2), ")" )
     )
abline(v=truth["LUE"],col=2)

hist( mcmcChain[,"tauV"]*1000., xlab="", ylab="",
      main=paste( "tauV \n( sd =",signif(sd(mcmcChain[,"tauV"]),2), ")" )
     )
abline(v=truth["tauV"],col=2)

hist( mcmcChain[,"tauS"]*10000, xlab="", ylab="",
      main=paste( "tauS \n( sd =",signif(sd(mcmcChain[,"tauS"]),2), ")" )
     )
abline(v=truth["tauS"],col=2)

hist( mcmcChain[,"Cv[1]"], xlab="", ylab="",
      main=paste( "Cv[1] \n( sd =",signif(sd(mcmcChain[,"Cv[1]"]),2), ")" )
     )
abline(v=truth["Cv"],col=2)

hist( mcmcChain[,"Cs[1]"], xlab="", ylab="",
      main=paste( "Cs[1] \n( sd =",signif(sd(mcmcChain[,"Cs[1]"]),2), ")" )
     )
abline(v=truth["Cs"],col=2)
hist( mcmcChain[,"cvCv"], xlab="", ylab="",
      main=paste( "cvCv \n( sd =",signif(sd(mcmcChain[,"cvCv"]),2), ")" )
     )
hist( mcmcChain[,"cvCs"], xlab="", ylab="",
      main=paste( "cvCs \n( sd =",signif(sd(mcmcChain[,"cvCs"]),2), ")" )
     )
##
## Plot model posterior mean outputs against obs
##
means_post <- apply( mcmcChain, 2, mean )
sd_post    <- apply( mcmcChain, 2, sd   )

i.Cv   <- sapply( 1:ndays, function(i) {
  which( colnames(mcmcChain)==paste("Cv[",i,"]",sep="") ) } )
ts_Cv  <- means_post[ i.Cv ] ; ts_Cv.sd <- sd_post[ i.Cv ]

i.Cs   <- sapply( 1:ndays, function(i) {
  which( colnames(mcmcChain)==paste("Cs[",i,"]",sep="") ) } )
ts_Cs  <- means_post[ i.Cs ] ; ts_Cs.sd <- sd_post[ i.Cs ]

i.NEE   <- sapply( 1:ndays, function(i) {
  which( colnames(mcmcChain)==paste("NEE[",i,"]",sep="") ) } )
ts_NEE  <- means_post[ i.NEE ] ; ts_NEE.sd <- sd_post[ i.NEE ]

## i.G   <- sapply( 1:n, function(i) {
##   which( colnames(mcmcChain)==paste("G[",i,"]",sep="") ) } )
## ts_G  <- means_post[ i.G ] ; ts_G.sd <- sd_post[ i.G ]

## i.PAR   <- sapply( 1:n, function(i) {
##   which( colnames(mcmcChain)==paste("PAR[",i,"]",sep="") ) } )
## ts_PAR  <- means_post[ i.PAR ] ; ts_PAR.sd <- sd_post[ i.PAR ]

par(mfrow=c(2,2))
plot(   days, ts_NEE, type='l', ylab="", main="NEE",
        ylim=c(min(ts_NEE,NEEobs,na.rm=TRUE),max(ts_NEE,NEEobs,na.rm=TRUE)) )
## arrows( days, ts_NEE-2*ts_NEE.sd, days, ts_NEE+2*ts_NEE.sd,
##         length=0.05, angle=90, code=3)
points( days, NEEobs, type='l', col="blue" )

plot(   days, ts_Cv, pch=16, ylab="", main="Cv",
        ylim=c(0,max(ts_Cv,Cvobs,na.rm=TRUE)) )
## arrows( days, ts_Cv-2*ts_Cv.sd, days, ts_Cv+2*ts_Cv.sd,
##         length=0.05, angle=90, code=3)
points( days, Cvobs, pch=1, col="blue" )

plot(   days, ts_Cs, pch=16, ylab="", main="Cs",
        ylim=c(0,max(ts_Cs,Csobs,na.rm=TRUE)) )
## arrows( days, ts_Cs-2*ts_Cs.sd, days, ts_Cs+2*ts_Cs.sd,
##         length=0.05, angle=90, code=3)
points( days, Csobs, pch=1, col="blue" )

## par(mfrow=c(1,1))
## plot(   days, ts_G, pch=16, ylab="", main="G",
##         ylim=c(min(ts_G,na.rm=TRUE),max(ts_G,na.rm=TRUE)) )
## arrows( days, ts_G-ts_G.sd, days, ts_G+ts_G.sd,
##         length=0.05, angle=90, code=3)

## par(mfrow=c(1,1))
## plot(   days, ts_PAR, pch=16, ylab="", main="PAR",
##         ylim=c(min(ts_PAR,na.rm=TRUE),max(ts_PAR,na.rm=TRUE)) )



