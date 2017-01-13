#' ---
#' title: "Analysis of RMS errors of calibrations with increasing unbalanced quantities of data"
#' author: David Cameron
#' date: '`r Sys.Date()`'
#' output: pdf_document
#' ---

library(BayesianTools)
library(dplyr)
library(purrr)
library(ggplot2)
library(gridExtra)

source("gridArrangeSharedLegend.R")
CLEAN.BUILD = FALSE
exptPath = "RDataWorking/"

# Setup VSEM and obs data
#' ## create PAR
set.seed(123)
ndays                 <- 2048
PAR                   <- VSEMcreatePAR(1:ndays)

refPars               <- VSEMgetDefaults()
nvar                  <-  nrow(refPars)+1

#' ## add SD
refPars[nvar,]          <- c(0.1, 0.001, 0.5)
rownames(refPars)[nvar] <- "error-sd"

#' ## calculate 'true' output and pseudodata
#' The "truth" data created from refPars$best
#' Obs are created by adding Gaussian noise to the "truth" data
referenceData         <- VSEM(refPars$best[1:(nvar-1)], PAR)
obs                   <- referenceData + rnorm(length(referenceData), sd = (abs(referenceData) + 1E-7) * 0.1)
#' set the minimum sd for NEE to 5 kgC /m^2 /day. This is to avoid very small uncertainty for NEE values close to zero
obs[,1] <- referenceData[,1] + rnorm(length(referenceData[,1]), sd = pmax((abs(referenceData[,1]) + 1E-7) * 0.1,0.0005))
## obs[,1] <- referenceData[,1] + rnorm(length(referenceData[,1]), sd = 0.0005)

#' ## parameters to be calibrated
defParms = c(1,3,5,6,9,10)

#' # RMS Error calibration plot
#' function to run a BC of the VSEM halfing the included data each time and finding the MAP point
csel <- map(2^{0:8},function(i) seq(1,2048,i))


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

#' ## Perfect Model
#' ### Likelihood
addPars <- refPars
newPars <- addPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c(predicted[isel,1] - obs[isel,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[isel,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[isel,2] - obs[isel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[isel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[isel,3] - obs[isel,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[isel,3])) + 0.0000001) * x[nvar], log = T)

  if (sum == FALSE) return(llValues)
  else return(sum(llValues1,llValues2,llValues3))
}

#' ### Prior
prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

#' ### Calculate MAP points
#+ results="hide"

MAP <- calcMAP("MAP.Rdata")

#' ## Perfect Model unbalanced data
#' ### Likelihood
addPars <- refPars
newPars <- addPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.
isLow = 2 ## which variables has few observations
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c(predicted[isel,1] - obs[isel,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[isel,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[isel,3] - obs[isel,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[isel,3])) + 0.0000001) * x[nvar], log = T)

  if (sum == FALSE) return(llValues)
  else return(sum(llValues1,llValues2,llValues3))
}

#' ### Prior
prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

#' ### Calculate MAP points
#+ results="hide"
MAPunbal <- calcMAP("MAPunbal.Rdata")

#' ## Model with error
#' ### Likelihood
addPars <- refPars
newPars <- addPars$best
## Model Error
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
parSel = c(defParms, nvar)
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c(predicted[isel,1] - obs[isel,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[isel,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[isel,2] - obs[isel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[isel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[isel,3] - obs[isel,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[isel,3])) + 0.0000001) * x[nvar], log = T)


  if (sum == FALSE) return(llValues)
  else return(sum(llValues1,llValues2,llValues3))
}

#' ### Prior
prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

#' ### Calculate MAP points
#+ results="hide"
MAPErr <- calcMAP("MAPErr.Rdata")

#' ## Model with error and unbalanced data
#' ### Likelihood
addPars <- refPars
newPars <- addPars$best
## Model Error
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
qparSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.
isLow = 2 ## which variables has few observations
likelihood <- function(x, sum = TRUE){
  x         <- createMixWithDefaults(x, newPars, parSel)
  predicted <- VSEM(x[1:(nvar-1)], PAR)

  diff       <- c(predicted[isel,1] - obs[isel,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[isel,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[isel,3] - obs[isel,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[isel,3])) + 0.0000001) * x[nvar], log = T)

  if (sum == FALSE) return(llValues)
  else return(sum(llValues1,llValues2,llValues3))
}

#' ### Prior
prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

#' ### Calculate MAP points
#+ results="hide"
MAPErrunbal <- calcMAP("MAPErrunbal.Rdata")

#' ## Create plots
calcRMS <- function(istep,fldno,tselect,fld,MAPfld){
    x          <- createMixWithDefaults(MAPfld[istep,], newPars, parSel)
    predicted  <- VSEM(x[-nvar], PAR)
    ofld       <- fld
    return(sqrt(mean((predicted[tselect,fldno] - ofld[tselect,fldno])**2)))
}

#' ### Plot against "truth"
#+ background='white', fig.height=7, fig.cap="Plot against 'truth'"
tt <- tibble::tibble ( istep=1:9)

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
    geom_line (aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data"))

p2 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CvRMS,colour="Perfect Model")) +
    geom_line (aes(y = CvRMS,colour="Perfect Model")) +
    geom_point(aes(y = CvRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_line (aes(y = CvRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_point(aes(y = CvRMSErr,colour="Model with Error")) +
    geom_line (aes(y = CvRMSErr,colour="Model with Error")) +
    geom_point(aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))+
    geom_line (aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))

p3 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CsRMS,colour="Perfect Model")) +
    geom_line (aes(y = CsRMS,colour="Perfect Model")) +
    geom_point(aes(y = CsRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_line (aes(y = CsRMSunbal,colour="Perfect Model UnBal Data")) +
    geom_point(aes(y = CsRMSErr,colour="Model with Error")) +
    geom_line (aes(y = CsRMSErr,colour="Model with Error")) +
    geom_point(aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data"))

## grid.arrange(p1,p2,p3)

## pdf("RMSErrPlotTruth.pdf")
grid_arrange_shared_legend(p1, p2, p3, ncol = 1, nrow = 3)
## dev.off()

#' ### Plot against obs
#+ background='white', fig.height=7, fig.cap="Plot against obs"
tt <- tibble::tibble ( istep=1:9)

newPars["Av"]  <- 1.0
newPars["Cr"]  <- 0.0
BCMAPRMSObs <- tt %>% mutate(NEERMSErrunbal=map_dbl(istep,function(i) calcRMS(i,fldno=1,tselect=1:ndays,fld=obs,MAPfld=MAPErrunbal))) %>%
                      mutate(CvRMSErrunbal =map_dbl(istep,function(i) calcRMS(i,fldno=2,tselect=obsSel,fld=obs,MAPfld=MAPErrunbal))) %>%
                      mutate(CsRMSErrunbal =map_dbl(istep,function(i) calcRMS(i,fldno=3,tselect=1:ndays,fld=obs,MAPfld=MAPErrunbal))) %>%
                      mutate(noObs=map_int(istep,function(i) length(csel[[i]])))

p1 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = NEERMSErrunbal,colour="Model with Error UnBal Data"))

p2 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))+
    geom_line (aes(y = CvRMSErrunbal,colour="Model with Error UnBal Data"))

p3 <- ggplot(data = BCMAPRMSObs, aes(x=noObs)) +
    geom_point(aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data")) +
    geom_line (aes(y = CsRMSErrunbal,colour="Model with Error UnBal Data"))

## pdf("RMSErrPlotObs.pdf")
grid_arrange_shared_legend(p1, p2, p3, ncol = 1, nrow = 3)
## dev.off()
