setwd("C:/Users/minunno/Downloads/TG15-drcWorking/TG15-drcWorking/")
runLab <- c("perModbalData","perModunbalData","errModbalData",
            "errModunbalData","perModbalDataBias","perModunbalDataBias",
            "errModunbalDataBias","errModunbalDatamodLike",
            "perModunbalDataBiasmodLike","errModunbalDataBiasmodLike")
runID <- c("run1", "run2", "run3", "run5",
           "run12", "run13", "run17",
           "run15","run16","run18")

library(BayesianTools)
library(tidyverse)
library(gridExtra)
library(ggcorrplot)
library(corrplot)
CLEAN.BUILD = FALSE

## PLOT.DIAGNOSTICS = TRUE
defParms = c(1,3,5,6,9,10) 
#exptPath = "RDataWorkingDellLong/"
exptPath = "~/research/TG15/RDataWorkingMSILong/"

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





library(BayesianTools)
library(tidyverse)
library(gridExtra)
CLEAN.BUILD = FALSE
## PLOT.DIAGNOSTICS = TRUE
defParms = c(1,3,5,6,9,10) 
#exptPath = "RDataWorkingDellLong/"
exptPath = "~/research/TG15/RDataWorkingMSILong/"


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

prior <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])
out1 <- fitVSEM("run1.RData")

MAP  <- calcMAP("MAP.Rdata")

plotParameters(out1,refPars)
plotOutputs(out1,refPars)

pCor1 <- plotCor(out1,parSel,refPars,runLab[1])
pCor1



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

plotParameters(out2,refPars)
plotOutputs(out2,refPars)

pCor2 <- plotCor(out2,parSel,refPars,runLab[2])
pCor2

## Likelihood: Gaussian 2048 obs for each of NEE, Cv and Cs
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

plotParameters(out3,refPars)
plotOutputs(out3,refPars)
pCor3 <- plotCor(out3,parSel,refPars,runLab[3])
pCor3

newPars <- refPars$best
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"] <- 0.0
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

prior         <- createUniformPrior(lower = refPars$lower[parSel], upper = refPars$upper[parSel])

out5 <- fitVSEM("run5.RData")

MAPErrunbal <- calcMAP("MAPErrunbal.Rdata")

plotParameters(out5,refPars)
plotOutputs(out5,refPars)

pCor5 <- plotCor(out5,parSel,refPars,runLab[4])
pCor5

obs.orig     <- obs
obs[,3] <- obs[,3] * 2.0

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

plotParameters(out12,refPars)
plotOutputs(out12,refPars)
obs <- obs.orig
pCor12 <- plotCor(out12,parSel,refPars,runLab[5])
pCor12


obs[,3] <- obs[,3] * 2.0

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

plotParameters(out13,refPars)
plotOutputs(out13,refPars)
obs <- obs.orig
pCor13 <- plotCor(out13,parSel,refPars,runLab[6])
pCor13


obs[,3] <- obs[,3] * 2.0

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

plotParameters(out17,refPars)
plotOutputs(out17,refPars)
obs <- obs.orig

pCor17 <- plotCor(out17,parSel,refPars,runLab[7])
pCor17



source("gridArrangeSharedLegend.R")

## csel <- map(2^{0:8},function(i) seq(1,2048,i))

## loadMAP <- function(fname){
##   load(paste(exptPath,fname,sep=""))
##   invisible(MAP)
## }

## MAP         <- loadMAP("MAP.Rdata")
## MAPunbal    <- loadMAP("MAPunbal.Rdata")
## MAPErr      <- loadMAP("MAPErr.Rdata")
## MAPErrunbal <- loadMAP("MAPErrunbal.Rdata")

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
addPars[nvar+3,]          <- c(0.0, -0.01, 0.01)
addPars[nvar+4,]          <- c(0.0, -1.0, 1.0)
rownames(addPars)[nvar+1]   <- "modmultNEE"
rownames(addPars)[nvar+2]   <- "modmultCs"
rownames(addPars)[nvar+3]   <- "modaddNEE"
rownames(addPars)[nvar+4]   <- "modaddCs"
newPars <- addPars$best
names(newPars) = row.names(addPars)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
parSel = c(defParms, nvar,nvar+1,nvar+2,nvar+3,nvar+4)
npar <- length(parSel)
obsSel <- c(1,202,390,550,750,920)*2
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)
  
  diff       <- c((predicted[,1]*x[nvar+1] + x[nvar+3]) - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(((predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+4]) - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)
  
  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])
out15 <- fitVSEM("run15.RData",iter=1200000, params=addPars)

plotParametersEsys(out15,addPars)
plotOutputsEsys(out15,addPars)
pCor15 <- plotCor(out15,parSel,addPars,runLab[8])
pCor15



obs[,3] <- obs[,3] * 2.0

addPars                   <- refPars
addPars[nvar+1,]          <- c(1.0, 0.1, 3.0)
addPars[nvar+2,]          <- c(1.0, 0.1, 3.0)
addPars[nvar+3,]          <- c(0.0, -0.01, 0.01)
addPars[nvar+4,]          <- c(0.0, -1.0, 1.0)
rownames(addPars)[nvar+1]   <- "modmultNEE"
rownames(addPars)[nvar+2]   <- "modmultCs"
rownames(addPars)[nvar+3]   <- "modaddNEE"
rownames(addPars)[nvar+4]   <- "modaddCs"
newPars <- addPars$best
names(newPars) = row.names(addPars)
parSel = c(defParms, nvar,nvar+1,nvar+2,nvar+3,nvar+4)
npar <- length(parSel)
obsSel <- c(1,202,390,550,750,920)*2
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaultsOld(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)
  
  diff       <- c((predicted[,1]*x[nvar+1] + x[nvar+3]) - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(((predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+4]) - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)
  
  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

out16 <- fitVSEM("run16.RData",params=addPars)
plotParametersEsys(out16,addPars)
plotOutputsEsys(out16,addPars)
obs <- obs.orig
pCor16 <- plotCor(out16,parSel,addPars,runLab[9])
pCor16



obs[,3] <- obs[,3] * 2.0

addPars                   <- refPars
addPars[nvar+1,]          <- c(1.0, 0.1, 3.0)
addPars[nvar+2,]          <- c(1.0, 0.1, 3.0)
addPars[nvar+3,]          <- c(0.0, -0.01, 0.01)
addPars[nvar+4,]          <- c(0.0, -1.0, 1.0)
rownames(addPars)[nvar+1]   <- "modmultNEE"
rownames(addPars)[nvar+2]   <- "modmultCs"
rownames(addPars)[nvar+3]   <- "modaddNEE"
rownames(addPars)[nvar+4]   <- "modaddCs"
newPars <- addPars$best
names(newPars) = row.names(addPars)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
parSel = c(defParms, nvar,nvar+1,nvar+2,nvar+3,nvar+4)
npar <- length(parSel)
obsSel <- c(1,202,390,550,750,920)*2
isLow = 2
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaultsOld(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)
  
  diff       <- c((predicted[,1]*x[nvar+1] + x[nvar+3]) - obs[,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(((predicted[1,3] + (predicted[,3] - predicted[1,3])*x[nvar+2]) + x[nvar+4]) - obs[,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[,3])) + 0.0000001) * x[nvar], log = T)
  
  ## if (sum == FALSE) return(llValues)
  ## else return(sum(llValues1,llValues2,llValues3))
  return(sum(llValues1,llValues2,llValues3))
}

prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

out18 <- fitVSEM("run18.RData",iter=1200000,params=addPars)
plotParametersEsys(out18,addPars)
plotOutputsEsys(out18,addPars)
obs <- obs.orig

pCor18 <- plotCor(out18,parSel,addPars,runLab[10])
pCor18

setwd("C:/Users/minunno/GitHub/TG15")

library(data.table)
library(ggplot2)
library(ggridges)

coef2 <- lm(pCor1$cors~pCor2$cors)
coef3 <- lm(pCor1$cors~pCor3$cors)
coef5 <- lm(pCor1$cors~pCor5$cors)
coef12 <- lm(pCor1$cors~pCor12$cors)
coef13 <- lm(pCor1$cors~pCor13$cors)
coef17 <- lm(pCor1$cors~pCor17$cors)

lmCoef <- data.table(
  inter=c(coef2$coefficients[1],coef3$coefficients[1],coef5$coefficients[1],
          coef12$coefficients[1],coef13$coefficients[1],coef17$coefficients[1]),
  slope=c(coef2$coefficients[2],coef3$coefficients[2],coef5$coefficients[2],
          coef12$coefficients[2],coef13$coefficients[2],coef17$coefficients[2]),
  runs = runLab[2:7]
)



corsData1 <- data.table(
  corCoef = c(pCor1$cors, pCor2$cors, pCor3$cors, pCor5$cors,
           pCor12$cors, pCor13$cors, pCor17$cors),
  runs = rep(runLab[1:7],each=length(pCor1$cors))
)
corsData2 <- data.table(
  corCoef = c(pCor15$cors, pCor16$cors, pCor18$cors),
  runs = rep(runLab[8:10],
             each=length(pCor15$cors))
)


pCorAllData <- cbind(corsData1[runs!=runLab[1]],
                     refRun=corsData1[runs==runLab[1]]$corCoef)



corsData1$run <- factor(corsData1$run,levels = runLab[1:7])
corsData2$run <- factor(corsData2$run,levels = runLab[8:10])
pCorAllData$run <- factor(pCorAllData$run,levels = runLab)

pCorAll <- ggplot(pCorAllData,aes(x=corCoef,y=refRun,col=runs))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlim(-1,1) + ylim(-1,1) +xlab("correltions all runs") + ylab("cor perModbalData") 



###correlation distribution
pRidgeCor <- ggplot(corsData1, 
  aes(x = corCoef, y = runs, fill = stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "", option = "C",limits=c(-1,1)) +
  labs(title = 'Correlation Coefficient distributions') 



pRidgeCor

pCorAll

png("pCor/summaryCorr.png",width = 800,height = 400)
grid_arrange_shared_legend(pCorAll,pRidgeCor,
                           ncol = 2, nrow = 1)
dev.off()


png("pCor/CorrFig1.png",width = 1000,height = 1000)
grid_arrange_shared_legend(pCor1$pCor, pCor2$pCor, pCor3$pCor, pCor5$pCor,
                           pCor12$pCor, pCor13$pCor, pCor17$pCor,
                           ncol = 3, nrow = 3)
dev.off()


png("pCor/CorrFig2.png",width = 1500)
grid_arrange_shared_legend(pCor15$pCor, pCor16$pCor, pCor18$pCor, ncol = 3, nrow = 1)
dev.off()






####parameter ditribution

pDistr <- data.table()


getSampleRun <- function(out,runX,startX=5e4){
  sampleX <- data.table(getSample(out,start=startX))
  setnames(sampleX, out$setup$names)
  sampleX[,run:=runX]
  return(sampleX)
}

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

sampleSet1 <- rbind(sample1,sample2,sample3,sample5,sample12,sample13,sample17)
sampleSet2 <- rbind(sample15,sample16,sample18)
sampleAll <- merge(sampleSet1,sampleSet2,all=T)
parNamAll <- out18$setup$names


levelX <- runLab
  # c("errModbalDataPar", "errModunbalDataBiasmodLikePar", "errModunbalDataBiasPar", "errModunbalDatamodLikePar",
  #           "errModunbalDataPar", "perModbalDataBiasOut", "perModbalDataPar", "perModunbalDataBiasmodLikePar",
  #           "perModunbalDataBiasPar","perModunbalDataOut")

sampleAll$run <- factor(sampleAll$run,levels = levelX)

pMAP <- addPars$best[parSel]
# pRidge <- list()
for(i in 1:length(pMAP)){
  # i=8
  pRidge <- ggplot(sampleAll, 
    aes(x = get(parNamAll[i]), y = run, fill = stat(x))) + 
    xlab(NULL)+ ylab(NULL) + 
    geom_vline(xintercept = pMAP[i],col="green") +
    geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    scale_fill_gradient2(midpoint = pMAP[i], low = "blue", high = "red",limits=as.numeric(addPars[parSel,2:3][i,])) +
    labs(title = parNamAll[i]) 
  print(pRidge)
  
  ggsave(paste0("pDistr/",parNamAll[i],".png"))
}

