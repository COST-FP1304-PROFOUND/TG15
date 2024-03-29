##
## Function to run the VSEM calibration and save the outputs if they don't already exist.
##
fitVSEM <- function(fname,iter=300000,params=refPars){
  if(!file.exists(paste(exptPath,fname,sep=""))|CLEAN.BUILD){
    bayesianSetup <- createBayesianSetup(likelihood, prior,best = newPars[parSel], names = rownames(params)[parSel])
    settings = list(iterations = iter)
    ## out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
    out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
    save(out,file=paste(exptPath,fname,sep=""))
  } else {
    load(paste(exptPath,fname,sep=""))
  }
  invisible(out)
}

##
## function to run the VSEM model needed for plotting model timeseries outputs from a posterior sample. 
##
runModel <- function(par,pool){
  x = createMixWithDefaultsOld(par, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)
  return(predicted[,pool])
}

##
## Alternative function to run the VSEM model needed for plotting model timeseries outputs from a posterior sample, when using the modified Likelihood that includes the systematic error terms 
##
runModelEsys2 <- function(par,pool){
  x = createMixWithDefaultsOld(par, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)
  if (pool == 1) predicted[,1] = predicted[,1]*x[nvar+1]  + x[nvar+4]
  if (pool == 2) predicted[,2] = (predicted[1,2] + (predicted[,2] - predicted[1,2])*x[nvar+3]) + x[nvar+6]
  
  if (pool == 3) predicted[,3] = (predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+5]
  return(predicted[,pool])
}

##
## Older version noiw unused
## runModelEsys <- function(par,pool){
##   x = createMixWithDefaultsOld(par, newPars, parSel)
##   predicted <- VSEM(x[1:(nvar-1)], PAR)
##   if (pool == 1) predicted[,1] = predicted[,1]*x[nvar+1]  + x[nvar+3]
##   if (pool == 3) predicted[,3] = (predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+4]
##   return(predicted[,pool])
## }

##
## Error functions needed to plot timeseries 
## 
errorFunction <- function(mean, par) rnorm(length(mean), mean = mean, sd = abs(mean)*par[length(defParms)+1])
errorFunctionNEE <- function(mean, par) rnorm(length(mean), mean = mean, sd = pmax((abs(referenceData[,1]) + 1E-7) * par[length(defParms)+1],0.0005))

##
## Functions used to plot posterior parameter histograms shown in supplementary material
##
plotPosteriors <- function(out,refPars){
  np = length(parSel)  ## number of parameters fit
  outM <- as.matrix(out$chain)[,1:np]
  colnames(outM) <- rownames(refPars)[parSel]
  for(i in 1:np){
    true = refPars[parSel[i],1]
    d = density(outM[,i])
    xlim = range(c(d$x,true))
    plot(d,xlim=xlim,main=colnames(outM)[i])
    abline(v=true,col=2,lwd=3)
  }
}
plotParameters <- function(out,thePars=     refPars){
    nmc = nrow(out$chain[[1]])
    
    out$chain <- window(out$chain,start=nmc/2,thin=(nmc/2)/5000)

    par(mfrow=c(4,2))
    plotPosteriors(out,thePars)
}
plotParametersEsys <- function(out,thePars=     refPars){
    nmc = nrow(out$chain[[1]])
    
    out$chain <- window(out$chain,start=nmc/2,thin=(nmc/2)/5000)

    par(mfrow=c(4,3))
    plotPosteriors(out,thePars)
}

##
## Function used to create Figure 2 plot
##
secondFig <- function(out,refPars){
    nmc = nrow(out$chain[[1]])

    par(mfrow = c(2,2))
    par(bg = "white")
     ttitles <- c("NEE","vegetative carbon (Cv)","soil carbon (Cs)")

    for(plotfld in 2:3){
        myModel <- function(x) runModel(x,plotfld)
        myObs = obs[,plotfld]
        if(exists("obsSel") & plotfld %in% isLow) myObs[-obsSel] <- NA
        plotTimeSeriesResultsOld(out, model = myModel, reference=referenceData[,plotfld], observed = myObs, error = errorFunction,start=nmc/2,plotResiduals = FALSE)
        lines(referenceData[,plotfld],col=3,lwd=1)
        if(plotfld==3) text(550,0.95*max(myObs),labels ="Pb (Eb, PbB)", cex=1.5)
        title(main = ttitles[plotfld])
    }

   myObs = obs[,1]
        if(exists("obsSel") & 1 %in% isLow) myObs[-obsSel] <- NA

    myModel <- function(x) runModel(x,1)
    plotTimeSeriesResultsOld(out, model = myModel, reference=referenceData[,1],observed = myObs, error = errorFunctionNEE,start=nmc/2, plotResiduals = FALSE)
        lines(referenceData[,1],col=3,lwd=1)
    title(main = ttitles[1])

}

##
## Function to create Supplementary material timeseries plots from a posterior sample 
##
plotOutputs <- function(out,refPars){
    nmc = nrow(out$chain[[1]])

    par(mfrow = c(2,2))
    par(bg = "white")
        myObs = obs[,1]
        if(exists("obsSel") & 1 %in% isLow) myObs[-obsSel] <- NA

    ttitles <- c("NEE","vegetative carbon (Cv)","soil carbon (Cs)")

    myModel <- function(x) runModel(x,1)
    plotTimeSeriesResultsOld(out, model = myModel, reference=referenceData[,1],observed = myObs, error = errorFunctionNEE,start=nmc/2, plotResiduals = FALSE)
        lines(referenceData[,1],col=3,lwd=1)
    title(main = ttitles[1])

    for(plotfld in 2:3){
        myModel <- function(x) runModel(x,plotfld)
        myObs = obs[,plotfld]
        if(exists("obsSel") & plotfld %in% isLow) myObs[-obsSel] <- NA
        plotTimeSeriesResultsOld(out, model = myModel, reference=referenceData[,plotfld], observed = myObs, error = errorFunction,start=nmc/2,plotResiduals = FALSE)
        lines(referenceData[,plotfld],col=3,lwd=1)
        title(main = ttitles[plotfld])
    }

}

##
## Function to create Supplementary material timeseries plots from a posterior sample, when using the modified Likelihood that includes the systematic error terms  
##
plotOutputsEsys2 <- function(out,refPars){
    nmc = nrow(out$chain[[1]])

    par(mfrow = c(2,2))
    par(bg = "white")
    myObs = obs[,1]
    if(exists("obsSel") & 1 %in% isLow) myObs[-obsSel] <- NA

    ttitles <- c("NEE","vegetative carbon (Cv)","soil carbon (Cs)")

    myModel <- function(x) runModelEsys2(x,1)
    plotTimeSeriesResultsOld(out, model = myModel,reference=referenceData[,1],  observed = myObs, error = errorFunctionNEE,start=nmc/2,plotResiduals = FALSE)
    lines(referenceData[,1],col=3,lwd=1)
    title(main = ttitles[1])

    for(plotfld in 2:3){
        myModel <- function(x) runModelEsys2(x,plotfld)
        myObs = obs[,plotfld]
        if(exists("obsSel") & plotfld %in% isLow) myObs[-obsSel] <- NA
        plotTimeSeriesResultsOld(out, model = myModel, reference=referenceData[,plotfld], observed = myObs, error = errorFunction,start=nmc/2,plotResiduals = FALSE)
        lines(referenceData[,plotfld],col=3,lwd=1)
        title(main = ttitles[plotfld])
    }
}
## Old version 
## plotOutputsEsys <- function(out,refPars){
##     nmc = nrow(out$chain[[1]])

##     par(mfrow = c(2,2))
##     par(bg = "white")
##     myObs = obs[,1]
##     if(exists("obsSel") & 1 %in% isLow) myObs[-obsSel] <- NA

##     ttitles <- c("NEE","vegetative carbon (Cv)","soil carbon (Cs)")

##     myModel <- function(x) runModelEsys(x,1)
##     plotTimeSeriesResultsOld(out, model = myModel,reference=referenceData[,1],  observed = myObs, error = errorFunctionNEE,start=nmc/2,plotResiduals = FALSE)
##     lines(referenceData[,1],col=3,lwd=1)
##     title(main = ttitles[1])

##     for(plotfld in 2:3){
##         myModel <- function(x) runModelEsys(x,plotfld)
##         myObs = obs[,plotfld]
##         if(exists("obsSel") & plotfld %in% isLow) myObs[-obsSel] <- NA
##         plotTimeSeriesResultsOld(out, model = myModel, reference=referenceData[,plotfld], observed = myObs, error = errorFunction,start=nmc/2,plotResiduals = FALSE)
##         lines(referenceData[,plotfld],col=3,lwd=1)
##         title(main = ttitles[plotfld])
##     }
## }

## RMS Error calibration plot
## function to run a BC of the VSEM halfing the included data each time and finding the MAP point
csel <- map(2^{0:8},function(i) seq(1,2048,i))

##
## Function to calculate the maximum a posteriori (MAP) parameter vector.
##
calcMAP <- function(fname){
  if(!file.exists(paste(exptPath,fname,sep=""))|CLEAN.BUILD){
    bayesianSetup <- createBayesianSetup(likelihood, prior,best = newPars[parSel], names = rownames(addPars)[parSel])
    settings = list(iterations = 300000)
    MAP <- array(0.0, c(length(csel),7))
    for (istep in 1:length(csel)){
       isel        <<- csel[[istep]]
       out         <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
       out$chain   <- window(out$chain,start=10000)
       MAP[istep,] <- out$chain[[1]][which.max(out$chain[[1]][,9]),1:7]
    }
    save(MAP,file=paste(exptPath,fname,sep=""))
  } else {
    load(paste(exptPath,fname,sep=""))
  }
  invisible(MAP)
}

##
## old unused functions
##
## vsemDiagnostics2 <- function(out,refPars){
##   ## thin
##   nmc = nrow(out$chain[[1]])
##   ## out$chain <- window(out$chain,start=nmc/2,thin=(nmc/2)/5000)
##   ## could replace above with getSample?

##     par(mfrow = c(2,2))
##     par(bg = "white")
##         myObs = obs[,1]
##         if(exists("obsSel") & 1 %in% isLow) myObs[-obsSel] <- NA

##     myModel <- function(x) runModel(x,1)
##     plotTimeSeriesResultsOld(out, model = myModel, observed = myObs, error = errorFunctionNEE,start=nmc/2, plotResiduals = FALSE)
##         lines(referenceData[,1],col=3,lwd=1)

##     for(plotfld in 2:3){
##         myModel <- function(x) runModel(x,plotfld)
##         myObs = obs[,plotfld]
##         if(exists("obsSel") & plotfld %in% isLow) myObs[-obsSel] <- NA
##         plotTimeSeriesResultsOld(out, model = myModel, observed = myObs, error = errorFunction,start=nmc/2,plotResiduals = FALSE)
##         lines(referenceData[,plotfld],col=3,lwd=1)
##     }

##   if(PLOT.DIAGNOSTICS){

##     #plot(out$chain)

##     out$chain <- window(out$chain,start=nmc/2,thin=(nmc/2)/5000)

##     par(mfrow=c(3,2))
##     plotPosteriors(out,refPars)

##     ## par(mfrow=c(1,1))
##     ## correlationPlot(out)
##   }
## }

## vsemDiagnostics2Esys <- function(out,refPars){
##   ## thin
##   nmc = nrow(out$chain[[1]])
##   ## out$chain <- window(out$chain,start=nmc/2,thin=(nmc/2)/5000)
##   ## could replace above with getSample?

##     par(mfrow = c(2,2))
##     par(bg = "white")
##     myObs = obs[,1]
##     if(exists("obsSel") & 1 %in% isLow) myObs[-obsSel] <- NA

##     ttitles <- c("NEE","vegetative carbon (Cv)","soil carbon (Cs)")

##     myModel <- function(x) runModelEsys(x,1)
##     plotTimeSeriesResultsOld(out, model = myModel, observed = obs[,1], error = errorFunctionNEE,start=nmc/2,plotResiduals = FALSE)
##     lines(referenceData[,1],col=3,lwd=1)
##     title(main = ttitles[1])


##     for(plotfld in 2:3){
##         myModel <- function(x) runModelEsys(x,plotfld)
##         myObs = obs[,plotfld]
##         if(exists("obsSel") & plotfld %in% isLow) myObs[-obsSel] <- NA
##         plotTimeSeriesResultsOld(out, model = myModel, observed = myObs, error = errorFunction,start=nmc/2,plotResiduals = FALSE)
##         lines(referenceData[,plotfld],col=3,lwd=1)
##         title(main = ttitles[plotfld])
        
##     }

##   if(PLOT.DIAGNOSTICS){

##     #plot(out$chain)

##     out$chain <- window(out$chain,start=nmc/2,thin=(nmc/2)/5000)

##     par(mfrow=c(3,2))
##     plotPosteriors(out,addPars)

##     ## par(mfrow=c(1,1))
##     ## correlationPlot(out)
##   }
## }
