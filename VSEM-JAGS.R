library(rjags)

n <- 365
NEEobs <- obs[1:n,1]
Cvobs  <- obs[1:n,2]
Csobs  <- obs[1:n,3]
days   <- seq(n)
PAR    <- (abs (sin(days/365 * pi)+ rnorm(n) *0.25)) *10

VSEM <- " model {
  Cr[1]  ~ dunif(  0, 200. )
  Cv[1]  ~ dunif(  0, 400. )
  Cs[1]  ~ dunif(  0, 1000. )
  ## update the states
  for (i in 2:n) {
    Cv[i]  <- Cv[i-1]  + Av*NPP[i-1] - Cv[i-1]/(tauV*1000.)
    Cr[i]  <- Cr[i-1]  + (1.0-Av)*NPP[i-1]  - Cr[i-1]/(tauR*1000.)
    Cs[i]  <- Cs[i-1]  + Cr[i-1]/(tauR*1000) + Cv[i-1]/(tauV*1000) - Cs[i-1]/(tauS*10000.)
  }
  for (i in 1:n) {
      ## Fluxes
      G[i]   <- PAR[i] * (LUE/1000.)* (1.0 - exp(-1.0*KEXT*LAR*Cv[i]))
      NPP[i]  <- (1.0-GAMMA)*G[i]
      NEE[i]  <- (Cs[i]/(tauS*10000.) + GAMMA*G[i]) - G[i]

      ## Calculate precisions
      tauNEE[i]  <- pow( cvNEE   * NEE[i]  , -2 )
      tauCv[i]  <- pow( cvCv   * Cv[i]  , -2 )
      tauCs[i]  <- pow( cvCs   * Cs[i]  , -2 )

      ## Likelihood
      NEEobs[i]    ~ dnorm( NEE[i], tauNEE[i] )
      Cvobs[i]     ~ dnorm( Cv[i], tauCv[i] )
      Csobs[i]     ~ dnorm( Cs[i], tauCs[i] )
  }
  KEXT   ~ dunif( 0.2, 1. )
  LAR    ~ dunif( 0.2, 3. )
  LUE    ~ dunif( 0.5,4.0)
  GAMMA  ~ dunif( 0.2, 0.6 )
  tauV   ~ dunif( 0.5,3. )
  tauR   ~ dunif( 0.5,3. )
  tauS   ~ dunif( 0.4,5. )
  Av     ~ dunif( 0.2, 1.0 )
  cvNEE  ~ dunif( 0.0, 1.0 )
  cvCv   ~ dunif( 0.0, 1.0 )
  cvCs   ~ dunif( 0.0, 1.0 )
}"

VSEMdata       <- list(n=n, NEEobs=NEEobs, Cvobs=Cvobs, Csobs=Csobs, PAR=PAR)
VSEMoutputs    <- c( "cvNEE", "cvCv", "cvCs", "cvCr", "KEXT", "LAR","LUE", 
                     "GAMMA","tauV","tauR","tauS","Av","NEE","G","PAR",
                     "Cv","Cs","Cr")
nadapt         <- 100; nbi <- 100; nit <- 1000
## nadapt         <- 1000; nbi <- 1000; nit <- 10000
jagsVSEM       <- jags.model ( textConnection(VSEM), data=VSEMdata, n.chains=1, n.adapt=nadapt)
update( jagsVSEM, n.iter=nbi )
codaSamples    <- coda.samples( jagsVSEM, var=VSEMoutputs, n.iter=nit, thin=20 )
mcmcChain      <- as.matrix( codaSamples )

par(mfrow=c(3,3))
hist( mcmcChain[,"KEXT"], xlab="", ylab="",
      main=paste( "KEXT \n( sd =",signif(sd(mcmcChain[,"KEXT"]),2), ")" ) 
     )
hist( mcmcChain[,"LAR"], xlab="", ylab="",
      main=paste( "LAR \n( sd =",signif(sd(mcmcChain[,"LAR"]),2), ")" ) 
     )
hist( mcmcChain[,"LUE"], xlab="", ylab="",
      main=paste( "LUE \n( sd =",signif(sd(mcmcChain[,"LUE"]),2), ")" ) 
     )
hist( mcmcChain[,"GAMMA"], xlab="", ylab="",
      main=paste( "GAMMA \n( sd =",signif(sd(mcmcChain[,"GAMMA"]),2), ")")
     )
hist( mcmcChain[,"tauV"], xlab="", ylab="",
      main=paste( "tauV \n( sd =",signif(sd(mcmcChain[,"tauV"]),2), ")" ) 
     )
hist( mcmcChain[,"tauR"], xlab="", ylab="",
      main=paste( "tauR \n( sd =",signif(sd(mcmcChain[,"tauR"]),2), ")" ) 
     )
hist( mcmcChain[,"tauS"], xlab="", ylab="",
      main=paste( "tauS \n( sd =",signif(sd(mcmcChain[,"tauS"]),2), ")" ) 
     )
hist( mcmcChain[,"Av"], xlab="", ylab="",
      main=paste( "Av \n( sd =",signif(sd(mcmcChain[,"Av"]),2), ")" ) 
     )


par(mfrow=c(3,2))
hist( mcmcChain[,"Cv[1]"], xlab="", ylab="",
      main=paste( "Cv[1] \n( sd =",signif(sd(mcmcChain[,"Cv[1]"]),2), ")" ) 
     )
hist( mcmcChain[,"Cr[1]"], xlab="", ylab="",
      main=paste( "Cr[1] \n( sd =",signif(sd(mcmcChain[,"Cr[1]"]),2), ")" ) 
     )
hist( mcmcChain[,"Cs[1]"], xlab="", ylab="",
      main=paste( "Cs[1] \n( sd =",signif(sd(mcmcChain[,"Cs[1]"]),2), ")" ) 
     )
hist( mcmcChain[,"cvCv"], xlab="", ylab="",
      main=paste( "cvCv \n( sd =",signif(sd(mcmcChain[,"cvCv"]),2), ")" ) 
     )
hist( mcmcChain[,"cvCr"], xlab="", ylab="",
      main=paste( "cvCr \n( sd =",signif(sd(mcmcChain[,"cvCr"]),2), ")" ) 
     )
hist( mcmcChain[,"cvCs"], xlab="", ylab="",
      main=paste( "cvCs \n( sd =",signif(sd(mcmcChain[,"cvCs"]),2), ")" ) 
     )


plot( mcmcChain[,"Av"], xlab="", ylab="",
      main=paste( "Av \n( sd =",signif(sd(mcmcChain[,"Av"]),2), ")")
     ,type='l'  )


means_post <- apply( mcmcChain, 2, mean )
sd_post    <- apply( mcmcChain, 2, sd   )

i.Cv   <- sapply( 1:n, function(i) {
  which( colnames(mcmcChain)==paste("Cv[",i,"]",sep="") ) } )
ts_Cv  <- means_post[ i.Cv ] ; ts_Cv.sd <- sd_post[ i.Cv ]

i.G   <- sapply( 1:n, function(i) {
  which( colnames(mcmcChain)==paste("G[",i,"]",sep="") ) } )
ts_G  <- means_post[ i.G ] ; ts_G.sd <- sd_post[ i.G ]

i.PAR   <- sapply( 1:n, function(i) {
  which( colnames(mcmcChain)==paste("PAR[",i,"]",sep="") ) } )
ts_PAR  <- means_post[ i.PAR ] ; ts_PAR.sd <- sd_post[ i.PAR ]


i.Cs   <- sapply( 1:n, function(i) {
  which( colnames(mcmcChain)==paste("Cs[",i,"]",sep="") ) } )
ts_Cs  <- means_post[ i.Cs ] ; ts_Cs.sd <- sd_post[ i.Cs ]

i.NEE   <- sapply( 1:n, function(i) {
  which( colnames(mcmcChain)==paste("NEE[",i,"]",sep="") ) } )
ts_NEE  <- means_post[ i.NEE ] ; ts_NEE.sd <- sd_post[ i.NEE ]


par(mfrow=c(1,1))
plot(   days, ts_Cv, pch=16, ylab="", main="Cv",
        ylim=c(0,max(ts_Cv,Cvobs,na.rm=TRUE)) )
arrows( days, ts_Cv-ts_Cv.sd, days, ts_Cv+ts_Cv.sd,
        length=0.05, angle=90, code=3)
points( days, Cvobs, pch=1, col="blue" )

par(mfrow=c(1,1))
plot(   days, ts_G, pch=16, ylab="", main="G",
        ylim=c(min(ts_G,na.rm=TRUE),max(ts_G,na.rm=TRUE)) )
arrows( days, ts_G-ts_G.sd, days, ts_G+ts_G.sd,
        length=0.05, angle=90, code=3)

par(mfrow=c(1,1))
plot(   days, ts_PAR, pch=16, ylab="", main="PAR",
        ylim=c(min(ts_PAR,na.rm=TRUE),max(ts_PAR,na.rm=TRUE)) )


par(mfrow=c(1,1))
plot(   days, ts_Cs, pch=16, ylab="", main="Cs",
        ylim=c(0,max(ts_Cs,Csobs,na.rm=TRUE)) )
arrows( days, ts_Cs-ts_Cs.sd, days, ts_Cs+ts_Cs.sd,
        length=0.05, angle=90, code=3)
points( days, Csobs, pch=1, col="blue" )

par(mfrow=c(1,1))
plot(   days, ts_NEE, pch=16, ylab="", main="NEE",
        ylim=c(min(ts_NEE,NEEobs,na.rm=TRUE),max(ts_NEE,NEEobs,na.rm=TRUE)) )
arrows( days, ts_NEE-ts_NEE.sd, days, ts_NEE+ts_NEE.sd,
        length=0.05, angle=90, code=3)
points( days, NEEobs, pch=1, col="blue" )



refPars$best[1:11]
names(VSEMcreatePAR(1:1000))
row.names(refPars)
