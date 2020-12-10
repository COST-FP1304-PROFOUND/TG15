## merged timeseries plots
library(purrr)
library(BayesianTools)
source("TG15-BayesianToolsOld.R")
source("helperFunctions.R")

load("obs.RData")

refPars = VSEMgetDefaults()
PAR = BayesianTools::VSEMcreatePAR(days = seq_len(nrow(obs)))
par(bg = "white")


## generate intervals
intervals = list()
experiments = c(1,2,3,5,12,13,17,15,16,18)
expName = c("Pb","Pu","Eb",NA,"Eu",NA,NA,NA,NA,NA,NA,"PbB","PuB",NA,"EuL","PuBL","EuB","EuBL")
for(fl in experiments){
  print(fl)
  
  fname = paste0("run",fl,".RData")
  if(file.exists(fname)) {
    load(fname)
  } else {
    print("skipped")
    next
  }
  
  intervals[[fl]] <- list()
  nmc = nrow(out$chain[[1]])
  
  for(plotfld in 2:3){
    error = errorFunction
    model <- function(x) runModel(x,plotfld)
    
    ## get parameter samples
    if(inherits(out,"bayesianOutput")){ 
      parMatrix = getSample(out, start = nmc/2)
    }else{ 
      if (class(out) == "matrix"){ parMatrix = out
      }else {stop("wrong type give to variable sampler")}
    }
    
    ## get intervals
    start = 1
    thin = 5000 ## number of samples RETAINED
    quantiles = c(0.025, 0.975)
    newPars <- refPars$best
    names(newPars) = row.names(refPars)
    defParms = parSel = which(names(newPars) %in% colnames(parMatrix))
    errSel = grep("error",colnames(parMatrix))
    nvar = nrow(refPars) + 1
    if(grepl("E",expName)[fl]){
      newPars["Av"] = 1
      newPars["Cr"] = 0
    }
    
    #pred <- getPredictiveIntervalsOld(parMatrix = parMatrix, model = model, thin = 1000, quantiles = c(0.025, 0.5, 0.975), error = error)
    
    pred = getPredictiveDistributionOld(parMatrix[,-errSel], model = model, thin = thin)
      ## some predictions are NaN, not sure why
    CI =  apply(pred,2,quantile,quantiles,na.rm=TRUE)
      #getCredibleIntervalsOld(sampleMatrix = pred, quantiles = quantiles)
    
    if(!is.null(error)){
      
      predDistr = pred
      for (i in 1:nrow(predDistr)){
        predDistr[i,] = error(mean = pred[i,], par = parMatrix[i,]) 
      }
      
      predInt = apply(predDistr,2,quantile,quantiles,na.rm=TRUE)
        #getCredibleIntervalsOld(sampleMatrix = predDistr, quantiles = quantiles)   
      CI = rbind(CI, predInt)
    }
    
    intervals[[fl]][[plotfld]] = CI
  } ## end plotfld
      
} ## end fl
    
#pdf("timeseries.pdf")
pdf("timeseriesSel.pdf",height=14)
par(mfrow = c(4,2))
titles =c("NEE","Cv","Cs")
obsUnbal <- c(1,202,390,550,750,920)*2.0
obsBal <- seq_along(PAR)
for(fl in experiments){
  for(plotfld in 2:3){
    myObs = obs[,plotfld]
    if(grepl("u",expName)[fl] && plotfld==2) obsSel = obsUnbal else obsSel = obsBal
    myObs[-obsSel] <- NA
    if(grepl("B",expName)[fl] && plotfld==3) myObs = myObs*2
    reference=referenceData[,plotfld]
    plotTimeSeriesOld(reference = reference, observed=myObs,confidenceBand = intervals[[fl]][[plotfld]][1:2,], predictionBand = intervals[[fl]][[plotfld]][3:4,],main=titles[plotfld])
    lines(referenceData[,plotfld],col=3,lwd=1)
    if(plotfld == 3) text(200,0.95*max(myObs),labels = as.character(expName[fl]),cex=2)
  }
}
dev.off()