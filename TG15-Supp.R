library(BayesianTools)
library(tidyverse)
library(gridExtra)

##
## Run all the calibrations needed for the study
##
source("runCalibrations.R")

# Supplementary material {-}

## Perfect model and balanced data Pb

obs <- obs.orig
newPars <- refPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
rm(obsSel)

plotParameters(out1,refPars)

plotOutputs(out1,refPars)

## Perfect model and unbalanced data Pu

obs <- obs.orig
newPars <- refPars$best
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2

plotParameters(out2,refPars)

plotOutputs(out2,refPars)

## Model with error and balanced data Eb

obs <- obs.orig
## Likelihood: Gaussian 2048 obs for each of NEE, Cv and Cs
newPars <- refPars$best
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"] <- 0.0
parSel = c(defParms, nvar)
rm(obsSel)
isLow = NULL

plotParameters(out3,refPars)

plotOutputs(out3,refPars)

## Model with error and unbalanced data Eu

obs <- obs.orig
newPars <- refPars$best
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"] <- 0.0
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2

plotParameters(out5,refPars)

plotOutputs(out5,refPars)

## Perfect model and balanced data with a multiplicative bias PbB
obs.orig     <- obs
obs[,3] <- obs[,3] * 0.8

newPars <- refPars$best
names(newPars) = row.names(refPars)
parSel = c(defParms, nvar)
rm(obsSel)
isLow = NULL

plotParameters(out12,refPars)

plotOutputs(out12,refPars)
obs <- obs.orig

## Perfect model and unbalanced data with a multiplicative bias PuB
obs[,3] <- obs[,3] * 0.8

newPars <- refPars$best
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2

plotParameters(out13,refPars)

plotOutputs(out13,refPars)
obs <- obs.orig

## Model with error and unbalanced data with a multiplicative bias EuB

obs[,3] <- obs[,3] * 0.8

newPars <- refPars$best
names(newPars) = row.names(refPars)
newPars["Av"]  <- 1.0
newPars["Cr"] <- 0.0
parSel = c(defParms, nvar)
obsSel <- c(1,202,390,550,750,920)*2.0
isLow = 2

plotParameters(out17,refPars)

plotOutputs(out17,refPars)
obs <- obs.orig

## Model with error and unbalanced perfect data with additive and multiplicative parameters to represent model error. EuL

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

plotParametersEsys(out15,addPars)

plotOutputsEsys2(out15,addPars)

## Perfect model and and unbalanced data with a multiplicative bias and additive and multiplicative parameters to represent the bias. PuBL

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

plotParametersEsys(out16,addPars)

plotOutputsEsys2(out16,addPars)
obs <- obs.orig

## Model with error and unbalanced data with a multiplicative bias and additive and multiplicative parameters to represent model error and the data bias. EuBL

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

plotParametersEsys(out18,addPars)

plotOutputsEsys2(out18,addPars)
obs <- obs.orig
