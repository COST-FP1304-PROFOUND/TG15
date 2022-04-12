library(BayesianTools)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(ggridges)
library(data.table)
library(ggplot2)
library(purrr)
CLEAN.BUILD = FALSE

## PLOT.DIAGNOSTICS = TRUE
defParms = c(1,3,5,6,9,10) 
exptPath = "RDataWorkingMSILongNew/"

if(!file.exists(paste(exptPath,"obs.RData",sep=""))|CLEAN.BUILD){
  set.seed(123)
  ndays                 <- 2048
  PAR                   <- VSEMcreatePAR(1:ndays)
  refPars               <- VSEMgetDefaults()
  nvar                  <-  nrow(refPars)+1
  ## add SD
  refPars[nvar,]          <- c(0.1, 0.001, 0.5)
  rownames(refPars)[nvar] <- "error-coeffVar"
  ## calculate 'true' output and pseudodata
  referenceData         <- VSEM(refPars$best[1:(nvar-1)], PAR)
  obs                   <- referenceData + rnorm(length(referenceData), sd = (abs(referenceData) + 1E-7) * refPars$best[nvar])
  ## set the minimum sd for NEE to 5 kgC /m^2 /day. This is to avoid very small uncertainty for NEE values close to zero
  obs[,1] <- referenceData[,1] + rnorm(length(referenceData[,1]), sd = pmax((abs(referenceData[,1]) + 1E-7) * refPars$best[nvar],0.0005))
  save(nvar,ndays,PAR,refPars,referenceData,obs,file=paste(exptPath,"obs.RData",sep=""))
  if(PLOT.DIAGNOSTICS){
    par(mfrow=c(3,1))
    for(i in 1:3){
      plot(obs[,i])
      lines(referenceData[,i],col=3,lwd=3)
    }
  }
} else {
  load(paste(exptPath,"obs.RData",sep=""))
  rownames(refPars)[nvar] <- "error-coeffVar"
}

source("TG15-BayesianToolsOld.R")
source("helperFunctions.R")

newPars <- refPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
rm(obsSel)
isLow = NULL
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaultsOld(x, newPars, parSel)
  predicted <- VSEM(x[-nvar], PAR)
  diff       <- c(predicted[,1] - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[,2] - obs[,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[,3] - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))

}

prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])
out1 <- fitVSEM("run1.RData")   

MAP  <- calcMAP("MAP.Rdata")

## Ridge plot

pDistr <- data.table()

runLab <- c("Pb","Pu","Eb",
            "Eu","PbB","PuB",
            "EuB","EuL",
            "PuBL","EuBL")

getSampleRun <- function(out,runX,startX=5e4){
  sampleX <- data.table(getSample(out,start=startX))
  setnames(sampleX, out$setup$names)
  sampleX[,run:=runX]
  return(sampleX)
}

loadPrevBT <- function(fname){
  load(paste0(exptPath,fname))
  invisible(out)
}

out1 <- loadPrevBT("run1.RData")
out2 <- loadPrevBT("run2.RData")
out3 <- loadPrevBT("run3.RData")
out5 <- loadPrevBT("run5.RData")
out12 <- loadPrevBT("run12.RData")
out13 <- loadPrevBT("run13.RData")
out17 <- loadPrevBT("run17.RData")
out15 <- loadPrevBT("run15.RData")
out16 <- loadPrevBT("run16.RData")
out18 <- loadPrevBT("run18.RData")

sample1 <- getSampleRun(out1,runLab[1])
sample2 <- getSampleRun(out2,runLab[2])
sample3 <- getSampleRun(out3,runLab[3])
sample5 <- getSampleRun(out5,runLab[4])
sample12 <- getSampleRun(out12,runLab[5])
sample13 <- getSampleRun(out13,runLab[6])
sample17 <- getSampleRun(out17,runLab[7])
sample15 <- getSampleRun(out15,runLab[8])
sample16 <- getSampleRun(out16,runLab[9])
sample18 <- getSampleRun(out18,runLab[10])

names(sample1)[7] <- "error-coeffVar"
names(sample2)[7] <- "error-coeffVar"
names(sample3)[7] <- "error-coeffVar"
names(sample5)[7] <- "error-coeffVar"

sampleSet1 <- rbind(sample1,sample2,sample3,sample5,sample12,sample13,sample17)

sampleSet2 <- rbind(sample15,sample16,sample18)

sampleAll <- merge(sampleSet1,sampleSet2,all=T)

parNamAll <- out18$setup$names

levelX <- runLab

sampleAll$run <- factor(sampleAll$run,levels = levelX)

addPars                   <- refPars
nvarTmp                   <- nrow(refPars)+1
addPars[nvarTmp+1,]          <- c(1.0, 0.1, 3.0)
addPars[nvarTmp+2,]          <- c(1.0, 0.1, 3.0)
addPars[nvarTmp+3,]          <- c(0.0, -0.01, 0.01)
addPars[nvarTmp+4,]          <- c(0.0, -1.0, 1.0)

parSel = c(defParms, nvarTmp)
pMAP <- addPars$best[parSel]

pRidgefn <- function(i){
  pRidge <- ggplot(sampleAll, 
    aes(x = get(parNamAll[i]), y = run, fill = stat(x))) + 
    xlab(NULL)+ ylab(NULL) + 
    geom_vline(xintercept = pMAP[i],col="green") +
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_gradient2(midpoint = pMAP[i], low = "blue", high = "red",limits=as.numeric(addPars[parSel,2:3][i,])) +
      labs(title = parNamAll[i])
    return(pRidge)
}

pRidgeCollect <- lapply(1:6,pRidgefn)
print(pRidgeCollect[[1]]+pRidgeCollect[[2]]+pRidgeCollect[[3]]+pRidgeCollect[[4]]+pRidgeCollect[[5]]+pRidgeCollect[[6]] + plot_layout(nrow = 3, byrow = FALSE))

addPars<-c()

## Perfect model Perfect data plot

firstFig(out1,refPars)

## Timseries plot

refPars               <- VSEMgetDefaults()
nvar                  <-  nrow(refPars)+1
## add SD
refPars[nvar,]          <- c(0.1, 0.001, 0.5)
rownames(refPars)[nvar] <- "error-sd"


intervals = list()
#experiments = c(1,2,3,5,12,13,17,15,16,18)
#expName = c("Pb","Pu","Eb",NA,"Eu",NA,NA,NA,NA,NA,NA,"PbB","PuB",NA,"EuL","PuBL","EuB","EuBL")
experiments = c(2,5,13,15,18)

expName = c(NA,"Pu",NA,NA,"Eu",NA,NA,NA,NA,NA,NA,NA,"PuB",NA,"EuL",NA,NA,"EuBL")
for(fl in experiments){
  #print(fl)
  if(fl > 12) rownames(refPars)[nvar] <- "error-sd"    
  fname = paste0(exptPath,"run",fl,".RData")
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
    if (grepl("L",expName)[fl]) {
     model <- function(x) runModelEsys2(x,plotfld)
    }else {    
     model <- function(x) runModel(x,plotfld)
    }
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

    if(grepl("L",expName)[fl]){
        addPars                   <- refPars
        addPars[nvar+1,]          <- c(1.0, 0.1, 2.0)
        addPars[nvar+2,]          <- c(1.0, 0.1, 2.0)
        addPars[nvar+3,]          <- c(1.0, 0.1, 2.0)
        addPars[nvar+4,]          <- c(0.0, -0.01, 0.01)
        addPars[nvar+5,]          <- c(0.0, -1.0, 1.0)
        addPars[nvar+6,]          <- c(0.0, -1.0, 1.0)
        rownames(addPars)[nvar]     <- "error-coeffVar"
        rownames(addPars)[nvar+1]   <- "modmultNEE"
        rownames(addPars)[nvar+2]   <- "modmultCs"
        rownames(addPars)[nvar+3]   <- "modmultCv"
        rownames(addPars)[nvar+4]   <- "modaddNEE"
        rownames(addPars)[nvar+5]   <- "modaddCs"
        rownames(addPars)[nvar+6]   <- "modaddCv"
        newPars <- addPars$best
        names(newPars) = row.names(addPars)
     #nvar = nrow(addPars) + 1
    } else{ 
      newPars <- refPars$best
      names(newPars) = row.names(refPars)
      
     #nvar = nrow(refPars) + 1
    }
   #errSel = grep("error",colnames(parMatrix))
    if(expName[fl]=="PuB") names(newPars)[nvar]="error-coeffVar"
    parSel = which(names(newPars) %in% colnames(parMatrix))
    if(grepl("E",expName)[fl]){
      newPars["Av"] = 1
      newPars["Cr"] = 0
    }
    
    #pred <- getPredictiveIntervalsOld(parMatrix = parMatrix, model = model, thin = 1000, quantiles = c(0.025, 0.5, 0.975), error = error)
    
   # pred = getPredictiveDistributionOld(parMatrix[,-errSel], model = model, thin = thin)
    pred = getPredictiveDistributionOld(parMatrix, model = model, thin = thin)
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

addPars<-c()
newPars<-c()
parMatrix<-c()

expName = c(NA,"Pu",NA,NA,"Eu",NA,NA,NA,NA,NA,NA,NA,"PuB",NA,"EuL",NA,NA,"EuBL")


par(mfrow = c(5,2))
par(mar = c(4, 4, 1.0, 1.0))
titles =c("NEE","Cv","Cs")
obsUnbal <- c(1,202,390,550,750,920)*2.0
obsBal <- seq_along(PAR)
for(fl in experiments){
  for(plotfld in 2:3){
    myObs = obs[,plotfld]
    if(grepl("u",expName)[fl] && plotfld==2) obsSel = obsUnbal else obsSel = obsBal
    myObs[-obsSel] <- NA
    if(grepl("B",expName)[fl] && plotfld==3) myObs = myObs*0.8
    reference=referenceData[,plotfld]
    plotTimeSeriesOld(reference = reference, observed=myObs,confidenceBand = intervals[[fl]][[plotfld]][1:2,], predictionBand = intervals[[fl]][[plotfld]][3:4,],main=titles[plotfld])
    lines(referenceData[,plotfld],col=3,lwd=1)
    if(expName[fl]=="PuB"){
        if(plotfld == 3) text(450,0.95*max(myObs),labels = paste0(as.character(expName[fl]),
                                                                  " (EuB)"),cex=1.2)
    } else if (expName[fl]=="EuBL"){
        if(plotfld == 3) text(550,0.95*max(myObs),labels = paste0(as.character(expName[fl]),
                                                                  " (PuBL)"),cex=1.2)
    } else {    
        if(plotfld == 3) text(200,0.95*max(myObs),labels = as.character(expName[fl]),cex=1.2)
    }
  }
}

refPars               <- VSEMgetDefaults()
nvar                  <-  nrow(refPars)+1
## add SD
refPars[nvar,]          <- c(0.1, 0.001, 0.5)
rownames(refPars)[nvar] <- "error-coeffVar"

newPars <- refPars$best
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[-nvar], PAR)

  diff       <- c(predicted[,1] - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[,3] - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

## Prior
prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

out2 <- fitVSEM("run2.RData")

MAPunbal <- calcMAP("MAPunbal.Rdata")

newPars <- refPars$best
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"] <- 0.0
parSel = c(defParms, nvar)
rm(obsSel)
isLow = NULL
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[-nvar], PAR)
  diff       <- c(predicted[,1] - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[,2] - obs[,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[,3] - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))

}

prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

out3 <- fitVSEM("run3.RData")

MAPErr <- calcMAP("MAPErr.Rdata")

MAPErrunbal <- calcMAP("MAPErrunbal.Rdata")

obs.orig     <- obs
obs[,3] <- obs[,3] * 0.8

newPars <- refPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
rm(obsSel)
isLow = NULL
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaultsOld(x, newPars, parSel)
  predicted <- VSEM(x[-nvar], PAR)
  diff       <- c(predicted[,1] - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[,2] - obs[,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[,3] - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

out12 <- fitVSEM("run12.RData")

obs <- obs.orig

obs[,3] <- obs[,3] * 0.8

newPars <- refPars$best
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaultsOld(x, newPars, parSel)
  predicted <- VSEM(x[-nvar], PAR)

  diff       <- c(predicted[,1] - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[,3] - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

out13 <- fitVSEM("run13.RData")

obs <- obs.orig

obs[,3] <- obs[,3] * 0.8

newPars <- refPars$best
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"] <- 0.0
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaultsOld(x, newPars, parSel)
  predicted <- VSEM(x[-nvar], PAR)

  diff       <- c(predicted[,1] - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[,3] - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

out17 <- fitVSEM("run17.RData")

obs <- obs.orig

## Diagnostic tool plot 1

source("gridArrangeSharedLegend.R")

calcRMS <- function(istep,fldno,tselect,fld,MAPfld){
    x          <- createMixWithDefaults(MAPfld[istep,], newPars, parSel)
    predicted  <- VSEM(x[-nvar], PAR)
    ofld       <- fld
    return(sqrt(mean((predicted[tselect,fldno] - ofld[tselect,fldno])**2)))
}

tt <- tibble::tibble ( istep=1:9)

newPars <- refPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
BCMAPRMSObs <- tt %>% mutate(NEERMSErrunbal=map_dbl(istep,function(i) calcRMS(i,fldno=1,tselect=1:ndays,fld=referenceData,MAPfld=MAPErrunbal))) %>%
                      mutate(CvRMSErrunbal =map_dbl(istep,function(i) calcRMS(i,fldno=2,tselect=1:ndays,fld=referenceData,MAPfld=MAPErrunbal))) %>%
                      mutate(CsRMSErrunbal =map_dbl(istep,function(i) calcRMS(i,fldno=3,tselect=1:ndays,fld=referenceData,MAPfld=MAPErrunbal))) %>%
                      mutate(NEERMSErr     =map_dbl(istep,function(i) calcRMS(i,fldno=1,tselect=1:ndays,fld=referenceData,MAPfld=MAPErr     ))) %>%
                      mutate(CvRMSErr      =map_dbl(istep,function(i) calcRMS(i,fldno=2,tselect=1:ndays,fld=referenceData,MAPfld=MAPErr     ))) %>%
                      mutate(CsRMSErr      =map_dbl(istep,function(i) calcRMS(i,fldno=3,tselect=1:ndays,fld=referenceData,MAPfld=MAPErr     ))) %>%
                      mutate(noObs=map_int(istep,function(i) length(csel[[i]])))

newPars["Av"]  <- 0.5
newPars["Cr"]  <- 3.0
BCMAPRMSObs <- BCMAPRMSObs %>%
                      mutate(NEERMSunbal   =map_dbl(istep,function(i) calcRMS(i,fldno=1,tselect=1:ndays,fld=referenceData,MAPfld=MAPunbal   ))) %>%
                      mutate(CvRMSunbal    =map_dbl(istep,function(i) calcRMS(i,fldno=2,tselect=1:ndays,fld=referenceData,MAPfld=MAPunbal   ))) %>%
                      mutate(CsRMSunbal    =map_dbl(istep,function(i) calcRMS(i,fldno=3,tselect=1:ndays,fld=referenceData,MAPfld=MAPunbal   ))) %>%
                      mutate(NEERMS        =map_dbl(istep,function(i) calcRMS(i,fldno=1,tselect=1:ndays,fld=referenceData,MAPfld=MAP        ))) %>%
                      mutate(CvRMS         =map_dbl(istep,function(i) calcRMS(i,fldno=2,tselect=1:ndays,fld=referenceData,MAPfld=MAP        ))) %>%
                      mutate(CsRMS         =map_dbl(istep,function(i) calcRMS(i,fldno=3,tselect=1:ndays,fld=referenceData,MAPfld=MAP        )))

p1 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = NEERMS,colour="Perfect Model")) +
    geom_line (aes(y = NEERMS,colour="Perfect Model")) +
    geom_point(aes(y = NEERMSunbal,colour="Perfect Model UnBal Data")) +
    geom_line (aes(y = NEERMSunbal,colour="Perfect Model UnBal Data")) +
    geom_point(aes(y = NEERMSErr,colour="Model with Error")) +
    geom_line (aes(y = NEERMSErr,colour="Model with Error")) +
    geom_point(aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data"))+
    labs(x = " ", y="RMS NEE")

p2 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CvRMS,colour="Perfect Model")) +
    geom_line (aes(y = CvRMS,colour="Perfect Model")) +
    geom_point(aes(y = CvRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_line (aes(y = CvRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_point(aes(y = CvRMSErr,colour="Model with Error")) +
    geom_line (aes(y = CvRMSErr,colour="Model with Error")) +
    geom_point(aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))+
    geom_line (aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))+
    labs(x = " ", y="RMS vegetative carbon")

p3 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CsRMS,colour="Perfect Model")) +
    geom_line (aes(y = CsRMS,colour="Perfect Model")) +
    geom_point(aes(y = CsRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_line (aes(y = CsRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_point(aes(y = CsRMSErr,colour="Model with Error")) +
    geom_line (aes(y = CsRMSErr,colour="Model with Error")) +
    geom_point(aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data"))+
    labs(x = "Number of observations included in the calibration", y="RMS soil carbon")

grid_arrange_shared_legend(p1, p2, p3, ncol = 1, nrow = 3)

## Diagnostic tool plot 2

tt <- tibble::tibble ( istep=1:9)

obsSel <- c(1,202,390,550,750,920)*2.0

newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
BCMAPRMSObs <- tt %>% mutate(NEERMSErrunbal=map_dbl(istep,function(i) calcRMS(i,fldno=1,tselect=1:ndays,fld=obs,MAPfld=MAPErrunbal))) %>%
                      mutate(CvRMSErrunbal =map_dbl(istep,function(i) calcRMS(i,fldno=2,tselect=obsSel,fld=obs,MAPfld=MAPErrunbal))) %>%
                      mutate(CsRMSErrunbal =map_dbl(istep,function(i) calcRMS(i,fldno=3,tselect=1:ndays,fld=obs,MAPfld=MAPErrunbal))) %>%
                      mutate(noObs=map_int(istep,function(i) length(csel[[i]])))

p1 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data")) +
    labs(x = " ", y="RMS NEE")

p2 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))+
    geom_line (aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))+ 
    labs(x = " ", y="RMS vegetative carbon")

p3 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data")) +
    labs(x = "Number of observations included in the calibration", y="RMS soil carbon")

grid_arrange_shared_legend(p1, p2, p3, ncol = 1, nrow = 3)

addPars                   <- refPars
addPars[nvar+1,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+2,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+3,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+4,]          <- c(0.0, -0.01, 0.01)
addPars[nvar+5,]          <- c(0.0, -1.0, 1.0)
addPars[nvar+6,]          <- c(0.0, -1.0, 1.0)
rownames(addPars)[nvar+1]   <- "modmultNEE"
rownames(addPars)[nvar+2]   <- "modmultCs"
rownames(addPars)[nvar+3]   <- "modmultCv"
rownames(addPars)[nvar+4]   <- "modaddNEE"
rownames(addPars)[nvar+5]   <- "modaddCs"
rownames(addPars)[nvar+6]   <- "modaddCv"
newPars <- addPars$best
names(newPars) = row.names(addPars)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
parSel = c(defParms, nvar,nvar+1,nvar+2,nvar+3,nvar+4,nvar+5,nvar+6)
npar <- length(parSel)
obsSel <- c(1,202,390,550,750,920)*2
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c((predicted[,1]*x[nvar+1] + x[nvar+4]) - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(((predicted[1,2] + (predicted[obsSel,2] - predicted[1,2])*x[nvar+3]) + x[nvar+6]) - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(((predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+5]) - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

out15 <- fitVSEM("run15.RData",iter=1200000, params=addPars)

obs[,3] <- obs[,3] * 0.8

addPars                   <- refPars
addPars[nvar+1,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+2,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+3,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+4,]          <- c(0.0, -0.01, 0.01)
addPars[nvar+5,]          <- c(0.0, -1.0, 1.0)
addPars[nvar+6,]          <- c(0.0, -1.0, 1.0)
rownames(addPars)[nvar+1]   <- "modmultNEE"
rownames(addPars)[nvar+2]   <- "modmultCs"
rownames(addPars)[nvar+3]   <- "modmultCv"
rownames(addPars)[nvar+4]   <- "modaddNEE"
rownames(addPars)[nvar+5]   <- "modaddCs"
rownames(addPars)[nvar+6]   <- "modaddCv"
newPars <- addPars$best
names(newPars) = row.names(addPars)
parSel = c(defParms, nvar,nvar+1,nvar+2,nvar+3,nvar+4,nvar+5,nvar+6)
npar <- length(parSel)
obsSel <- c(1,202,390,550,750,920)*2
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c((predicted[,1]*x[nvar+1] + x[nvar+4]) - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(((predicted[1,2] + (predicted[obsSel,2] - predicted[1,2])*x[nvar+3]) + x[nvar+6]) - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(((predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+5]) - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

out16 <- fitVSEM("run16.RData",iter=1200000,params=addPars)

#pdf('par16.pdf')
#plot(out16)
#dev.off()

#pdf('out16.pdf')
#plotOutputsEsys2(out16,addPars)
#dev.off()

obs <- obs.orig

obs[,3] <- obs[,3] * 0.8

addPars                   <- refPars
addPars[nvar+1,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+2,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+3,]          <- c(1.0, 0.1, 2.0)
addPars[nvar+4,]          <- c(0.0, -0.01, 0.01)
addPars[nvar+5,]          <- c(0.0, -1.0, 1.0)
addPars[nvar+6,]          <- c(0.0, -1.0, 1.0)
rownames(addPars)[nvar+1]   <- "modmultNEE"
rownames(addPars)[nvar+2]   <- "modmultCs"
rownames(addPars)[nvar+3]   <- "modmultCv"
rownames(addPars)[nvar+4]   <- "modaddNEE"
rownames(addPars)[nvar+5]   <- "modaddCs"
rownames(addPars)[nvar+6]   <- "modaddCv"
newPars <- addPars$best
names(newPars) = row.names(addPars)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
parSel = c(defParms, nvar,nvar+1,nvar+2,nvar+3,nvar+4,nvar+5,nvar+6)
npar <- length(parSel)
obsSel <- c(1,202,390,550,750,920)*2
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c((predicted[,1]*x[nvar+1] + x[nvar+4]) - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(((predicted[1,2] + (predicted[obsSel,2] - predicted[1,2])*x[nvar+3]) + x[nvar+6]) - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(((predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+5]) - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)

  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

out18 <- fitVSEM("run18.RData",iter=1200000, params=addPars)
