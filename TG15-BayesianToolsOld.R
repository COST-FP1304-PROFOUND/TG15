createMixWithDefaultsOld <- function (pars, defaults, locations) 
{
    out = defaults
    out[locations] = pars
    return(out)
}

getCredibleIntervalsOld <- function(sampleMatrix, quantiles = c(0.025, 0.975)){

  x = matrix (ncol = ncol(sampleMatrix), nrow = length(quantiles))
  rownames(x) = quantiles
    
  for (i in 1:length(quantiles)){
    x[i,] = apply(sampleMatrix,2,function(x)quantile(x,probs=quantiles[i]))  
  } 
  return(x)
}

getPredictiveDistributionOld<-function(parMatrix, model, thin = 1000){
  
  # Do thinning if wanted and neccessary
  if (thin != F & nrow(parMatrix) > 2*thin){
    sel = round(seq(1,nrow(parMatrix), len = thin ))
    parMatrixSel = parMatrix[sel,]
  }else{
    parMatrixSel = parMatrix
  }
  
  # calculate predictions
  
  run1 = model(parMatrixSel[1,])
  
  out = matrix(NA, ncol = length(run1), 
               nrow = nrow(parMatrixSel))
  
  out[1,] = run1
  
  for (i in 2:nrow(parMatrixSel)){
    out[i,] = model(parMatrixSel[i,])
  }
  return(out)
}

getPredictiveIntervalsOld<-function(parMatrix, model, thin = 1000, quantiles = c(0.025, 0.975), error = NULL){
  pred = getPredictiveDistributionOld(parMatrix, model = model, thin = thin)
  out = getCredibleIntervalsOld(sampleMatrix = pred, quantiles = quantiles)
  
  if(!is.null(error)){
    
    predDistr = pred
    for (i in 1:nrow(predDistr)){
      predDistr[i,] = error(mean = pred[i,], par = parMatrix[i,]) 
    }

    predInt = getCredibleIntervalsOld(sampleMatrix = predDistr, quantiles = quantiles)   
    out = rbind(out, predInt)
  }
  
  return(out)
}

plotTimeSeriesOld <- function(reference=NULL, observed = NULL, predicted = NULL, x = NULL, confidenceBand = NULL, predictionBand = NULL, xlab = "Time", ylab = "Observed / predicted", ...){
  
  ylim = range(reference, observed, predicted, confidenceBand, predictionBand,na.rm=T)
  
  if (is.null(x)){
    if(!is.null(observed)) x = 1:length(observed)
    else if(!is.null(predicted)) x = 1:length(predicted)
    else stop("either observed or predicted must be supplied")
  }
  
  len = length(x)
  
  plot(x, ylim = ylim, type = "n", xlab = xlab, ylab = ylab, ...)
  
  if(!is.null(predictionBand)) polygon(c(1:len,len:1),c(predictionBand[1,],predictionBand[2,len:1]),col="moccasin",border=NA)
  
  if(!is.null(confidenceBand)) polygon(c(1:len,len:1),c(confidenceBand[1,],confidenceBand[2,len:1]),col="#99333380",border=NA)    
    
  if(!is.null(predicted)) lines(predicted, col = "red")
  if(!is.null(observed)) points(observed, col = "black", pch = 3, cex = 0.6)
  
}

plotTimeSeriesResidualsOld <- function(residuals, x = NULL, main = "residuals"){
  
  ylim = range(residuals)
  
  if (is.null(x)){
    x  = 1:length(residuals)
  }
  barplot(residuals)
}

plotTimeSeriesResultsOld <- function(sampler, model, reference,observed, error = NULL, plotResiduals = T, start = 1){
  
  if(inherits(sampler,"bayesianOutput")) parMatrix = getSample(sampler, start = start)
  else if (class(sampler) == "matrix") parMatrix = sampler
  else stop("wrong type give to variable sampler")
  
  pred <- getPredictiveIntervalsOld(parMatrix = parMatrix, model = model, thin = 1000, quantiles = c(0.025, 0.5, 0.975), error = error)
  
  plotTimeSeriesOld(reference=reference, observed = observed, predicted = pred[2,], confidenceBand = pred[c(1,3),], predictionBand = pred[c(4,6),] )
}
