##
## switches to determine whether to make a complete rebuild and whether to make diagnostic plots
##
CLEAN.BUILD = FALSE
## PLOT.DIAGNOSTICS = TRUE

##
## Define which VSEM parameters to include in the calibration
##
defParms = c(1,3,5,6,9,10)

##
## Folder to store saved outputs from calibrations. These saved outputs avoid the need to rerun calbrations later.
##
exptPath = "RDataWorkingMSILongNew/"

##
## Create virtual observations from VSEM model outputs plus Gaussian noise. 
##
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

##
## Load required functions from an older version of Bayesian Tools 
##
source("TG15-BayesianToolsOld.R")

##
## Load a series of R functions to help with running and plotting the results from the calibrations 
##
source("helperFunctions.R")


##
## run calibrations
##
##
## run1: Perfect model and balanced data Pb
##
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

##
## run2: Perfect model and unbalanced data Pu
##
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

##
## run3: Model with error and balanced data Eb
##
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

##
## run5: Model with error and unbalanced data Eu
##
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

##
## run12: Perfect model and balanced data with a multiplicative bias PbB
##
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

##
## run13: Perfect model and unbalanced data with a multiplicative bias PuB
##
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

##
## run17: Model with error and unbalanced data with a multiplicative bias EuB
##
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

##
## run15: Model with error and unbalanced perfect data with additive and multiplicative parameters to represent model error. EuL
##
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

##
## run16: Perfect model and and unbalanced data with a multiplicative bias and additive and multiplicative parameters to represent the bias. PuBL
##
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

##
## run18: Model with error and unbalanced data with a multiplicative bias and additive and multiplicative parameters to represent model error and the data bias. EuBL
##
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
