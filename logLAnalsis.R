#' ---
#' title: "log likelihood analysis"
#' author: David Cameron
#' date: '`r Sys.Date()`'
#' output: pdf_document
#' ---

library(BayesianTools)
library(dplyr)
library(coda)
library(purrr)
library(ggplot2)


#'
#' # Aim
#' Th aim of this  analysis is to determine the influence of the observational data on subsequent Bayesian
#' calibrations.
#' The idea is to quantify how the likelihood will change when increasing quantities of
#' observations are added to the calibration. This is assessed by sampling from the posterior created below with
#' just a few obs (6) calculating the "inter-quartile range" of the resultant likelihood function as progressively
#' more data are added. The log IQR gives a measure of the range of the likelihood function and hence a metric to
#' assess the influence of the obs on the likelihood. This analysis in effect quantifies the influence of the obs
#' on the calibration.
#'

CLEAN.BUILD = FALSE
exptPath = "RDataWorking/"

#' # Setup VSEM and obs data
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

#' # Calibration with just 6 obs (obsSel)
#' The idea here is to make an initial calbration of VSEM with a balanced number of obs for each output (NEE, Cv and Cs).
#' ## Likelihood function
addPars        <- refPars
newPars        <- addPars$best
names(newPars) <- row.names(refPars)
parSel         <- c(defParms, nvar)
obsSel         <- c(1,202,390,550,750,920)*2.
isLow          <-  2 ## which variables has few observations
likelihood <- function(x, sum = TRUE){
  x          <- createMixWithDefaults(x, newPars, parSel)
  predicted  <- VSEM(x[1:(nvar-1)], PAR)
  diff       <- c(predicted[obsSel,1] - obs[obsSel,1])
  llValues1  <- dnorm(diff, sd = pmax((abs(c(predicted[obsSel,1])) + 0.0000001) * x[nvar],0.0005), log = T)
  diff       <- c(predicted[obsSel,2] - obs[obsSel,2])
  llValues2  <- dnorm(diff, sd = (abs(c(predicted[obsSel,2])) + 0.0000001) * x[nvar], log = T)
  diff       <- c(predicted[obsSel,3] - obs[obsSel,3])
  llValues3  <- dnorm(diff, sd = (abs(c(predicted[obsSel,3])) + 0.0000001) * x[nvar], log = T)
  if (sum == FALSE) return(llValues)
  else return(sum(llValues1,llValues2,llValues3))
}

#' ## prior
prior         <- createUniformPrior(lower = addPars$lower[parSel], upper = addPars$upper[parSel])

#' ## run MCMC
#+ results="hide"
bayesianSetup <- createBayesianSetup(likelihood, prior,best = newPars[parSel], names = rownames(addPars)[parSel])
settings      <- list(iterations = 300000)
fname         <- paste(exptPath,"balancedSixObs.Rdata",sep="")
if(!file.exists(fname)|CLEAN.BUILD){
  balencedSixObs <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
  save(balencedSixObs, file=fname)
} else {
  load(fname)
}

#' ## Remove burn-in check gelman and take a sample
balencedSixObs$chain <- window(balencedSixObs$chain,start=60000)
gelman.diag(balencedSixObs$chain[,1:7])
Sample    <- getSample(balencedSixObs,thin=40,coda=T)
nSample   <- nrow(Sample[[1]])

#' # Likelihood analysis
#'
#+ background='white', fig.height=7, fig.cap=" "
csel <- map(2^{0:8},function(i) seq(1,2048,i))

calcLogLIQR <- function(istep,iob,fld){
    llout <- vector("numeric", nSample)
    isel  <- csel[[istep]]
    for (i in 1:nSample){
       samplePars <- Sample[[1]][i,]
       x          <- createMixWithDefaults(samplePars, newPars, parSel)
       predicted  <- VSEM(x[-nvar], PAR)
       diff       <- c(predicted[isel,iob] - obs[isel,iob])
       if (fld=="NEE"){
          llValues   <- dnorm(diff, sd = pmax((abs(c(predicted[isel,iob])) + 0.0000001) * x[nvar],0.0005), log = T)
          ## llValues   <- dnorm(diff, sd = (abs(c(predicted[isel,iob])) + 0.0000001) * x[nvar], log = T)
       } else {
          llValues   <- dnorm(diff, sd = (abs(c(predicted[isel,iob])) + 0.0000001) * x[nvar], log = T)
       }
       llout[i]   <- sum(llValues)
    }
    return(IQR(llout)/1.349)
}

tt   <- tibble::tibble ( istep=1:9,NEE = 1, Cv = 2, Cs = 3)
logLIRQ <- tt %>%
        mutate(CvlogL =map2_dbl(istep,Cv, function(istep,iob) calcLogLIQR(istep,iob,fld="Cv")))  %>%
        mutate(CslogL =map2_dbl(istep,Cs, function(istep,iob) calcLogLIQR(istep,iob,fld="Cs")))  %>%
        mutate(NEElogL=map2_dbl(istep,NEE,function(istep,iob) calcLogLIQR(istep,iob,fld="NEE"))) %>%
        mutate(noObs=map_int(istep,function(i) length(csel[[i]])))

ggplot(data = logLIRQ, aes(x=noObs)) +
    geom_point(aes(y = NEElogL,colour="NEElogL")) +
    geom_line(aes(y = NEElogL,colour="NEElogL")) +
    geom_point(aes(y = CvlogL ,colour="CvlogL"))  +
    geom_line(aes(y = CvlogL ,colour="CvlogL"))  +
    geom_point(aes(y = CslogL ,colour="CslogL"))  +
    geom_line(aes(y = CslogL ,colour="CslogL"))  +
    xlab("Number of obs inc. in logL calc")       +
    ylab("inter-quartile range logL")             +
    theme(legend.title=element_blank())
