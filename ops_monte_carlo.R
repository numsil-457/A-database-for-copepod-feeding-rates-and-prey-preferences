#### Script to the prey preferences of copepods
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Change path to R script directory

library(plyr)
library(nloptr)
library(parallel)

## Functions to extract the OPS from the feeding rates database
# Selectivity model from Wirtz (2012)
sel = function(lesd, s, lops){
  x = exp( -s* (lesd - lops)**2 )
  return(x)
}

# >> Using NLOPT
RMSE = function( par, ESD, log_y, n_ESD ){
  #OPS=par[1] ; s=par[2] ; p = par[3]
  y.mod = sel(ESD, par[2], par[1]) * par[3]
  
  #rmse = sqrt( sum( y.mod - y )**2 / length(ESD) ) # RMSE
  #rmse = sum( abs(y.mod - y) ) / length(ESD) # MAE
  rmse = sum( abs( log(y.mod+1e-6) - log_y ) ) / n_ESD # MAE
  
  #rmse = sum( abs(y.mod - y ) / (y + 1e-5) ) / length(ESD) * 100  # MAPE
  
  #rmse = sum( abs(y.mod - y)**2 ) # SSE
  
  return( rmse )
} 

# RMSE for bimodal
RMSE2 = function( par, ESD, log_y, n_ESD ){
  
  y.mod = sel(ESD, par[2], par[1]) * par[3] +
    sel(ESD, par[5], par[4]) * par[6] # Bimodal kernel
  
  #rmse = sqrt( sum( y.mod - y )**2 / length(ESD) ) # RMSE
  rmse = sum( abs( log(y.mod+1e-6) - log_y ), na.rm=T ) / n_ESD # MAE
  return( rmse )
} 

FWHM = function( peak, res, low ){
  
  # Searching for the indexes at half maximum
  hm = res$y[peak] / 2 # Half maximum
  ind.hm = which( res$y - hm > 0 )
  
  #### Select the indices that have a direct relation to "peak" (1 increment)
  # These belong to the same peak
  # -> We do this by comparing with a series with 1 indices increment
  ind.peak = which( ind.hm == peak )
  sercomp = c( (1 - ind.peak):0, 1:(length(ind.hm) - ind.peak) ) + peak
  ind.hm = ind.hm[ which(sercomp == ind.hm) ] # Selecting indices that follow the 1-inc order
  
  ind.hm = c( ind.hm[1], tail(ind.hm, 1) ) # First indice and last 
  
  ## Check if the full peak is present
  # If yes, ok
  # If not, we assume that the peak is gaussian, and multiply the width between the peak and
  # the existing side by 2. If the estimated width is larger than what was computed with
  # the indices, then we take the estimation
  
  fwestim = 0
  if( ind.hm[1] == 0 | ind.hm[2] == length(res$y) ){      # If the peak is not full
    side.missing = which(ind.hm %in% c(0, length(res$y))) # Missing side
    halfwidth = ind.hm[-side.missing]:peak
    fwestim = abs( res$x[ tail(halfwidth, 1) ] - res$x[halfwidth[1]] ) * 2   # Estimated half, assuming a gaussian
  }
  
  fw = res$x[ind.hm[2]] - res$x[ind.hm[1]] # Full width
  
  fw = max(fw, fwestim) # Take the largest width
  
  return( fw / 2.355 )
} # Full Width at Half Maximum rule

PIFREQ = function( mu1, mu2, sd1, sd2, h1, h2 ){ # Gaussian parameters
  dx = 1e-2
  xseq = seq(-1e2, 1e2, by = dx)
  
  g1 = h1 * exp( -0.5 * (xseq - mu1)**2 / sd1**2 )
  g1 = sum(g1 * dx) # Frequency for the peak 1
  
  g2 = h2 * exp( -0.5 * (xseq - mu2)**2 / sd2**2 )
  g2 = sum(g2 * dx) # Frequency for the peak 1
  
  # Compute pi.(1 - pi)
  pif = g1 * g2 / (g1 + g2)**2
  return(pif)
} # Compute the contribution in frequency of each peak

check.peak = function(tk){
  res = list( x = tk$relative.size, y = tk$p )
  
  if(nrow(tk) >= 3){ # At least 3 points for the test
    tkdiff = diff(tk$p)
    d1 = tkdiff[-length(tkdiff)]
    d2 = tkdiff[-1]
    cross = d1 * d2 < 0
    
    # Check if the data forms a peak already, it it does there is no need to smooth
    peaks = which(cross & d1 > 0) + 1
    lows  = which(cross & d1 < 0) + 1
    
    if(length(peaks) == 0) peaks = which.max(tk$p) # Last check
  } else {
    peaks = which.max(tk$p)
    lows  = integer(0)
  }
  
  return( list(res, peaks, lows) )} # Find the number of maxima in a distribution

PI.score = function(peaks, lows, res, ip, jp){ # Compute the PI score for bimodality, combining the FWHM and PIfreq functions
  # This function doesn't use a smooth kernel
  
  # Modes
  mu1 = res$x[ peaks[ip] ] 
  mu2 = res$x[ peaks[jp] ]
  
  # FWHM for both peaks
  # Keep the original points to compare
  xout = c(res$x, seq(min(res$x), max(res$x), length.out=1000))
  res.approx = approx(x=res$x, y=res$y, xout=xout[order(xout)] )
  
  # Searching for the indexes at half maximum
  half_max = max(res$y) / 2 # Half maximum
  ind.half_max = which( res.approx$y - half_max > 0 ) # All the spectrum that is above
  
  height1 = res$y[peaks[ip]]
  height2 = res$y[peaks[jp]]
  
  if(height1 <= half_max || height2 <= half_max) return(0)
  
  # Use which.min distance instead of exact float equality
  ind.p1 = which.min(abs(res.approx$y - height1))  
  ind.p2 = which.min(abs(res.approx$y - height2))
  
  # Find the points attached to each peak using the approximation
  get.hm.indices = function(ind.p, ind.half_max){
    ind.peak = which(ind.half_max == ind.p)
    if(length(ind.peak) == 0) return(NULL)
    sercomp = unique(c((1 - ind.peak):0, 
                       (1 - length(ind.half_max) == ind.peak):(length(ind.half_max) - ind.peak)) + ind.p)
    ind.hm = ind.half_max[which(sercomp == ind.half_max)]
    c(ind.hm[1], tail(ind.hm, 1))
  }
  
  ind.hm1 = get.hm.indices(ind.p1, ind.half_max)
  ind.hm2 = get.hm.indices(ind.p2, ind.half_max)
  
  if(is.null(ind.hm1) || is.null(ind.hm2)) return(NaN)
  
  ## Compute the FWHM
  sd1 = (res.approx$x[ind.hm1[2]] - res.approx$x[ind.hm1[1]]) / 2.355 # Full width
  sd2 = (res.approx$x[ind.hm2[2]] - res.approx$x[ind.hm2[1]]) / 2.355
  
  if(sd1 <= 0 || sd2 <= 0) return(NaN)
  
  # PIFREQ
  dx = 1e-2
  xseq = seq(-1e2, 1e2, by = dx)
  
  # Closed-form Gaussian integral: height * sd * sqrt(2π)
  g1 = height1 * sd1 * sqrt(2 * pi) # Frequency for the peak 1
  g2 = height2 * sd2 * sqrt(2 * pi) # Frequency for the peak 2
  
  # Compute pi.(1 - pi)
  pif = g1 * g2 / (g1 + g2)**2
  
  # Bimodality score
  bi = sqrt(pif)*abs(mu1/sd1 - mu2/sd2)
  
  return(bi)
}

init.par = function(lesd, p){
  opsini = lesd[ which.max(p) ]
  sini = 1/sd(lesd, na.rm=T)
  suggestpar = c(opsini, sini, max(p))
  return(suggestpar)
}

fit.model = function(lesd, p, suggestpar, func, n_restarts = 4, maxeval = 1e5L){
  maxeval = maxeval / n_restarts # Adapt the evaluation time on the restarts
  #print(maxeval)
  # Algorithms crs2lm (fatser) or isrers then lbfgs
  
  # Set boundaries for the parameter search
  n_groups = length(suggestpar) %/% 3
  lb = ub = suggestpar
  
  for (g in seq_len(n_groups)) {
    i = (g - 1) * 3 + 1:3
    lb[i[1]] = suggestpar[i[1]] * 0.8 # OPS
    ub[i[1]] = suggestpar[i[1]] * 1.2
    lb[i[2]] = 0.1 # s
    ub[i[2]] = 50
    lb[i[3]] = suggestpar[i[3]] * 0.5 # Height
    ub[i[3]] = suggestpar[i[3]] * 2
  }
  
  # --- Parallelized ISRES restarts ---
  
  # Generate perturbed starting points for each restart
  start_points <- lapply(seq_len(n_restarts), function(i) {
    if (i == 1) return(suggestpar)          # keep one clean start
    suggestpar * runif(length(suggestpar), 0.85, 1.15)  # perturb others
  })
  
  n_cores = min(n_restarts, parallel::detectCores() - 1)
  
  results = parallel::mclapply(
    start_points,
    function(x0) {
      tryCatch(
        isres(x0    = x0,
              fn    = func,
              lower = lb,
              upper = ub,
              maxeval = maxeval,
              ESD    = lesd,
              log_y  = log(p+1e-6),
              n_ESD  = length(lesd)),
        error = function(e) NULL   # skip failed restarts
      )
    },
    mc.cores = n_cores
  )
  
  # Filter failed runs, pick best result by objective value
  results  = Filter( Negate(is.null), results)
  best     = results[[ which.min(sapply(results, `[[`, "value")) ]]
  
  s  = sign(best$par)
  lb = best$par * 0.9 * (1 + (1 - s) * 0.1)
  ub = best$par * 1.1 * (1 - (1 - s) * 0.1)
  
  mod = lbfgs(x0 = best$par, fn = func,
              lower = lb, upper = ub,
              ESD=lesd, log_y=log(p+1e-6), n_ESD = length(lesd) )
  
  return( mod$par )
}

#### TO CHECK
# Check the parameter initialization for the model, maybe it can be further optimized

## Get the feeding rates
pref = read.csv(file = "copepod_feeding_rates.csv")

## Homogenize the stage variable

# Look for experiments with adult copepods
adult_groups = c("F", "M", "NA", "CVI", "A")
ind.adults = unique( c(grep( paste(adult_groups, collapse='|'), pref$stage ),
                       which(pref$stage=='' | is.na(pref$stage))) )
pref$stage[ind.adults] = "A"

## Classification of life stages in broad groups
pref$stage.group=NA
indN=grep("N", pref$stage) ; indC=grep("C", pref$stage)
pref$stage.group[ indN ] = "N"
pref$stage.group[ indC ] = "C"
pref$stage.group[ intersect(indN, indC) ] = "NC"
pref$stage.group[ pref$stage == 'A' ] = "A"

## Set up the groups of data points
pref$unit = NA

# Use first the Fmax, clearance rate
# Most studies report the clearance rate. If the Imax is calculated with a linear model, it is usually
# biased towards large prey.
# Then, use the Imax if the Fmax is not available (happens for some rare studies)
# Then, use the selectivity if available (main grazing variable in the article)

pref$unit[!is.na(pref$Fmax.at.15.degreeC..ml.mgC.1.h.1.)] = 'Fmax.at.15.degreeC..ml.mgC.1.h.1.'
pref$unit[is.na(pref$Fmax.at.15.degreeC..ml.mgC.1.h.1.)] = 'Imax.at.15.degreeC..mugC.mugC.1.h.1.'
pref$unit[!is.na(pref$selectivity)] = 'selectivity'
pref$Ukern = paste(pref$species, pref$stage, pref$unit, sep="_") # Unique species/stage

## Remove the NAs in the points
pref = pref[-which( is.na(pref$prey.esd) | is.na(pref$unit) | is.na(pref$pred.esd) ),]

## Should treat all the data from Hansen et al (1994) separately - The OPS is already extracted
pref.hansen = pref[ which(pref$secondary.reference == 'Hansen et al (1994)'), ]
pref = pref[ which(pref$secondary.reference != 'Hansen et al (1994)' |
                     is.na(pref$secondary.reference)), ]
# Some points are alone, but the OPS is known. Should not analyze this in the algorithm, just do a mean
# per precise life stage

pref.hansen$ref = paste(pref.hansen$primary.reference, pref.hansen$secondary.reference, sep='; ')
hansen.data = ddply( pref.hansen, .(Ukern, species, stage, ref), summarize,
                                    phylum='Copepod',
                                    pred.esd = mean(pred.esd, na.rm = T),
                                    pred.esd.sd = sd(pred.esd, na.rm = T),
                                    ops = mean(prey.esd[which(selectivity==1)], na.rm = T), # Data
                                    prey.esd.min = min(prey.esd), 
                                    prey.esd.max = max(prey.esd))

## Check which species of a similar stage are present (Nauplius, Copepodite, Adult)
# Then, aggregate by size, if not too far (e.g. if the group variance is bigger that the distance?)
Ukern.discarded = c()
for(ui in unique(pref$Ukern)){
  t = pref[ pref$Ukern == ui, ]
  if(length( unique(t$prey.esd) ) <= 2){ # Check which points are not in a group
    Ukern.discarded = c(Ukern.discarded, unique(t$Ukern) )
  }
} 

## Pull the alone data points in other groups, if not too far 
# Needs to be the same species, same life group (N/C/A) and a similar size
ui.seq = c() ; others.seq = c() # Check if it is the only point, if so not useful to keep
for(ui in Ukern.discarded){
  ind.alone = which(pref$Ukern == ui)
  t.alone = pref[ind.alone,]
  species.i = unique(t.alone$species)
  stage.i = unique(t.alone$stage.group)
  unit.i = unique(t.alone$unit)
  
  t.others = pref[ pref$stage.group == stage.i & pref$species == species.i & 
                     pref$unit == unit.i & pref$Ukern != ui, ]
  ui.seq = c(ui.seq, ui) ; others.seq = c(others.seq, nrow(t.others))
  
  t.variance = ddply(t.others, .(Ukern), summarize, m.esd = mean(pred.esd), sd.esd = sd(pred.esd))
  
  ## Check the distance of every point and assign to the lowest distance
  # And to the group where |mu - i| < std
  for(ii in 1:nrow(t.alone)){ 
    dist = abs(t.variance$m.esd - t.alone$pred.esd[ii])/(t.variance$sd.esd+1e-4)
    ind = intersect( which( dist < 1 ), which.min(dist) )
    
    if(length(ind)>0){pref$Ukern[ ind.alone[ii] ] = t.variance$Ukern[ind]}
    # Assign the point to this group
  }
}
# Maybe search for another criteria; sd is sometimes very low or equal to 0 but the size is close

#data.frame(ui=ui.seq, others=others.seq) # See which species did not find a group

## OPS extraction and Monte Carlo simulation on the Q10 correction

# Clean up the temperatures that are not numeric
pref$temperature[pref$temperature == "23-27"] = 25 # Mean

cols_to_convert = c("temperature")
pref[cols_to_convert] = lapply(pref[cols_to_convert], as.numeric)

# Q10 ranges from 1 to 4.9 in Tyrrell et al 2019 (JPR)
# Q10  = 2.8 (Hansen et al, 1997)
q10 = function( r1, t1, t2 = 15, Q10 = 2.8 ){
  return( r1 * Q10**( (t2 - t1)/10 ) )
}

MC.perturb = function(t){
  #Q10i = rnorm(1, mu.q10, sd.q10) # Normal distribution; test other distributions and see the effect
  Q10i = runif(1, 1.9, 3.8) # Unimodal
  #print(Q10i)
  
  # Modify the rates using the new Q10
  t$Imax.at.15.degreeC..mugC.mugC.1.h.1. = q10(t$Imax.ref, 
                                               t1=as.numeric(t$temperature), t2=15,
                                               Q10=Q10i) # New Q10
  
  t$Fmax.at.15.degreeC..ml.mgC.1.h.1. = q10(t$Fmax.ref, 
                                            t1=as.numeric(t$temperature), t2=15,
                                            Q10=Q10i) # New Q10
  
  # Check which variable to use for the OPS compilation
  var.t = names(t)[ which(names(t) == unique(t$unit)) ]
  t$p = t[[var.t]]
  t$p = t$p/max(t$p, na.rm=T)
  
  return(t)
}

fit.selectivity.model = function(t, output.df){
  # For plotting
  x.seq = seq( log(0.1), 8, length.out = 500) # lESD
  
  # Checking the number of potential modes in the data
  # Do a mean of preference on duplicated prey - otherwise the kernel smoothing doesn't work?
  t = t[ which( !duplicated(t[c('relative.size', 'p')]) ), ] # Remove duplicates, if any
  
  t$row_id = seq_len(nrow(t))
  tk = ddply( t, .(relative.size), summarize, i.max = row_id[which(p == max(p, na.rm = T) )] ) 

  t = t[tk$i.max,] # Avoid the repeated points at one prey size
  t = t[!duplicated(t$relative.size),] # In case different prey of same size have the same p, very rare
  
  # The dataset needs to be ordered by ESD for the peak detection, etc.
  res = check.peak(t)
  lows = res[[3]] ; peaks = res[[2]] ; res = res[[1]]
  
  ## Check for bimodality
  if(length(peaks)>1){
  combos = combn(seq_along(peaks), 2, simplify = FALSE) # Trying all combinations of peaks
  BI = do.call(rbind, lapply(combos, function(idx) {
    bi = PI.score(peaks, lows, res, idx[1], idx[2])
    data.frame(bi = bi, p1 = peaks[idx[1]], p2 = peaks[idx[2]])
  }))
  
  # Remove NaNs, happens when a peak is not really one
  BI = BI[!is.nan(BI$bi), ]
  # BI was defined in Wang et al (2009) to determine the "visually bimodal" distributions
  # The BI limit is 1.1
  
  }else{
    BI = data.frame()
  }
  
  ## Define the densities to fit with the model
  t$gp   = NA # Group affectation for the model
  t$peak = NA # Used for estimating the initial OPS
  if( length(BI$bi) >= 1 & suppressWarnings( max(BI$bi) )>= 1.1 ){ # Bimodal
    
    # The points are split using the lows
    # The peaks with the smallest BI are grouped together
    if( nrow(BI) > 1 ){
      bip = BI[ rev( order(BI$bi) ), ][2,]
    }else{
      bip = BI
    }
    
    bip = BI[which.max(BI$bi),]                 # Highest bimodality index
    low = lows[ lows > bip$p1 & lows < bip$p2 ] # Find the lows in the distribution
    
    # If there are 2 lows in between, take the deepest one
    low = low[ which.min(t$p[low]) ]
    
    t$gp[ t$relative.size <= res$x[low] ] = 1 # Peaks separated by the low
    t$gp[ t$relative.size >  res$x[low] ] = 2
    
    t$peak[ t$gp == 1 ] = res$x[ bip$p1 ] # Detected peaks
    t$peak[ t$gp == 2 ] = res$x[ bip$p2 ]
    bimod = T
    
  }else{ # Unimodal
    t$gp = 1 # No group
    t$peak = res$x[ peaks[ which.max( res$y[ peaks ] ) ] ] # Highest detected peak
    
    bimod = F # Is bimodality detected ?
  } 
  
  # Computing the model on each group
  RSS = 0 # AIC estimation
  model.output = data.frame() # Parameter estimation of the model
  aic.bi = 1e6
  aic.uni = 1e6
  
  # Fit the model of selectivity
  # For bimodal
  if( max(t$gp)>1 ){ 
    t1 = t[ t$gp == 1, ]
    t2 = t[ t$gp == 2, ]
    init.par1 = init.par(t1$relative.size, t1$p)
    init.par2 = init.par(t2$relative.size, t2$p)
    par.bi = c(init.par1, init.par2)
    
    params = fit.model(t$relative.size, t$p, par.bi, RMSE2)
    # Fix the order of the parameters so they are the same order as the rest
    ind.pars = order(params[c(1,4)])
    opsi1 = params[1+(ind.pars[1]-1)*3] ; si1 = params[2+(ind.pars[1]-1)*3] ; pi1 = params[3+(ind.pars[1]-1)*3]
    opsi2 = params[1+(ind.pars[2]-1)*3] ; si2 = params[2+(ind.pars[2]-1)*3] ; pi2 = params[3+(ind.pars[2]-1)*3]
    
    sel.bi = sel(x.seq, si1, opsi1)*pi1 + sel(x.seq, si2, opsi2)*pi2
    
    # RMSE
    rmse.bi = RMSE2( params, t$relative.size, log(t$p+1e-6), nrow(t) )
    
    # AIC
    resid = t$p - ( sel(t$relative.size, si1, opsi1)*pi1 + sel(t$relative.size, si2, opsi2)*pi2 )
    RSS = sum(resid**2)
    
    # AIC estimation for bimodal
    n = nrow(t)
    Log_L = -n/2 * log(2*pi) - n/2 * log( sum(RSS) /n) - n/2 # Log Likelihood
    aic.bi = 2*( 6 ) - 2*Log_L 
  }else{
    
    # For Unimodal
    par.uni = init.par(t$relative.size, t$p)
    
    params = fit.model(t$relative.size, t$p, par.uni, RMSE)
    opsi = params[1] ; si = params[2] ; pi = params[3]
    
    sel.uni = sel(x.seq, si, opsi)*pi
    
    # RMSE
    rmse.uni = RMSE( params, t$relative.size, log(t$p+1e-6), nrow(t)  )
    
    # AIC
    resid = t$p - sel(t$relative.size, si, opsi)
    RSS = sum(resid**2)
    
    # AIC estimation for bimodal
    n = nrow(t)
    Log_L = -n/2 * log(2*pi) - n/2 * log( sum(RSS)/n ) - n/2 # Log Likelihood
    aic.uni = 2*( 2 ) - 2*Log_L 
  }
  
  # Save information
  if( max(t$gp)>1 ){ # Or AIC, but probably not correct Log likelihood
    model.output = data.frame( bi = T,
                               ops1 = exp( par.bi[1] ), # Data
                               ops2 = exp( par.bi[4] ),
                               rmse = rmse.bi
                               ) # Model
    
  }else{
    model.output = data.frame( bi = F,
                               ops1 = exp( par.uni[1] ), # Data
                               ops2 = NA,
                               rmse = rmse.uni
                               ) 
  }
  
  output.df = rbind( output.df, model.output )
}

bs.test = function( pref, nDraws=100, test.type = 'None', n.cores = max(1, detectCores() - 4) ){
  #nDraws = 50 # Draws for Monte Carlo and Bootstrap
  
  ## Q10 = 2.8 c [1.9-3.8] at 95% (Saiz et al 2022)
  mu.q10 = 2.8 ; sd.q10 = (3.8-1.9)/4 # mu+2.sd - (mu-2.sd)/4 = sd
  
  total.df = data.frame() # Store the final dataset
  for(ixp in 1:length(unique(pref$Ukern)) ){
    t = pref[ pref$Ukern == unique(pref$Ukern)[ixp], ]
    #t = pref[ pref$Ukern == "Calanus finmarchicus_A_Fmax.at.15.degreeC..ml.mgC.1.h.1.", ]
    #t = pref[ pref$Ukern == "Tortanus spp_CI-III_selectivity", ]
    
    # Remove any NAs
    unique(t$unit)[!is.na(unique(t$unit))] 
    
    # Searches for data only if Imax at more than 3 prey
    if( nrow(t) > 2 & length( unique(t$prey.esd) ) > 2 ){
      print(t$Ukern[1])
      
      t.refs = paste(unique(t$primary.reference), collapse = '; ')
      
      # Choose the X-axis variable (OPS or OPS:ESD ratio)
      t$relative.size = log( t$prey.esd ) # lESD of prey, simplifies the code in the loops below

      # Variables already know prior simulation
      Ukern = unique(t$Ukern)
      species = unique(t$species)
      stage = unique(t$stage)
      pred = mean(t$pred.esd, na.rm = T)
      pred.sd = sqrt( sum(t$pred.esd.sd**2, na.rm=T)*1/(nrow(t)**2) )
      # Variable for all the references used for the OPS measurement?
      
      # Extract the prey ESD limits being looked at
      esd.min = min(t$prey.esd, na.rm=T)
      esd.max = max(t$prey.esd, na.rm=T)
      
      # If error propagation test
      if( test.type %in% c('Q10', 'BS', 'Q10BC') ){
        #output.df = data.frame() # Store data for a species
        
        ## Q10 MC 
        # Q10 sensitivity loop
        t$Imax.ref = q10(t$Imax.at.15.degreeC..mugC.mugC.1.h.1., 
                         t1=15, t2=as.numeric(t$temperature),
                         Q10=2.8) # Retroacting the Q10
        
        t$Fmax.ref = q10(t$Fmax.at.15.degreeC..ml.mgC.1.h.1., 
                         t1=15, t2=as.numeric(t$temperature),
                         Q10=2.8) # Retroacting the Q10
        
        if( !(unique(t$unit) %in% c("selectivity")) ){
          # Run the loop only if Fmax and Imax have a value
          # Selectivity does not need to be resampled, as all rates are unitless and already at the same temperature
          
          mc.results = mclapply(seq_len(nDraws), function(i) {
            t.mc = MC.perturb(t) # MC draw here
            fit.selectivity.model(t.mc, data.frame())
          }, mc.cores = n.cores)
          
          output.df = do.call(rbind, mc.results)
          
        }else{ # Just one run to find the OPS, but no sd, as there is no influence of Q10 on the data
          # Check which variable to use for the OPS compilation
          var.t = names(t)[ which(names(t) == unique(t$unit)) ]
          t$p = t[[var.t]]
          t$p = t$p/max(t$p, na.rm=T)
          
          output.df = fit.selectivity.model(t, data.frame())
        }
        
        ## Bootstrap resample
        # for(i in 1:nDraws){
        #   # idx = sample(seq_len(nrow(t)), replace = TRUE) # Bootstrap indices
        #   # t =  t[idx,] # New sample
        
        # Check which variable to use for the OPS compilation
            # var.t = names(t)[ which(names(t) == unique(t$unit)) ]
            # t$p = t[[var.t]]
            # t$p = t$p/max(t$p, na.rm=T)
        
        #   if( length(unique(t$prey,size)) > 1){ # Checks that at least two different points
        #     output.df = fit.selectivity.model(t, output.df)}
        #   
        #   if(i%%10 == 0){ print(i) }
        # }
        # 
        # ## Monte Carlo (Q10) and Bootstrap resample
        # for(i in 1:(nDraws**2)){
        #   # Monte Carlo Q10
        #   
        #   # Bootstrapt
        #   # idx = sample(seq_len(nrow(t)), replace = TRUE) # Bootstrap indices
        #   # t =  t[idx,] # New sample
        #   if( length(unique(t$prey)) > 1){ # Checks that at least two different points
        #     output.df = fit.selectivity.model(t, output.df)}
        #   
        #   if(i%%10 == 0){ print(i) }
        # }
      
        ## Summary statistics of the dataset
        output.df = output.df[ c("ops1", "ops2", 
                                 "rmse"
                                 ) ] # Keep the same dataframe format
        
        # Mean and sd
        output.m = apply(output.df, 2, mean, na.rm=T) # Apply the function and correct format + names
        output.m = setNames(as.data.frame(matrix(output.m, nrow=1)), names(output.df))
        
        output.sd = apply(output.df, 2, sd, na.rm=T)
        output.sd = setNames(as.data.frame(matrix(output.sd, nrow=1)), names(output.df))
        
        # Compiling all of the information
        names(output.m) = paste(names(output.m), '.m', sep='') 
        names(output.sd) = paste(names(output.sd), '.sd', sep='') 
        
        # Quantiles
        qt.probs = c(0.05, 0.95) 
        output.qt = apply(output.df, 2, quantile, probs = qt.probs, na.rm=T)
        
        qt.df = data.frame()
        for(i in 1:nrow(output.qt)){ # Split the output.qt and rename
          qt.dfi = setNames(as.data.frame(matrix(output.qt[i,], nrow=1)), names(output.df))
          names(qt.dfi) = paste(names(qt.dfi), qt.probs[i]*100, sep='.') # Store the quantiles
          
          if(i==1){qt.df = qt.dfi
          }else{ qt.df = cbind(qt.df, qt.dfi) }
        }
  
        # Line with all the variables
        output.refs = data.frame(species=species,
                                 phylum='Copepod',
                                 stage=stage, 
                                 pred.esd=pred, 
                                 pred.esd.sd=pred.sd,
                                 prey.esd.min=esd.min, 
                                 prey.esd.max=esd.max)
        spec.line = cbind(output.refs, output.m, output.sd, qt.df)#, refs)
        
        total.df = rbind(total.df, spec.line)
      
      ## Normal simulation, no resample
      }else{
        # Check which variable to use for the OPS compilation
        var.t = names(t)[which( names(t) == unique(t$unit) )]
        t$p = t[[var.t]]
        t$p = t$p/max(t$p, na.rm=T)
        
        output.df = fit.selectivity.model(t, data.frame())
        
        output.refs = data.frame(species=species,
                                 phylum='Copepod',
                                 stage=stage, 
                                 pred.esd=pred, 
                                 pred.esd.sd=pred.sd,
                                 prey.esd.min=esd.min, 
                                 prey.esd.max=esd.max
                                 )
        
        spec.line = cbind( output.refs, output.df )
        total.df = rbind(total.df, spec.line)
      }
    } # Finished fitting the OPSs
    
    if(ixp%%1 == 0 | ixp == length(unique(pref$Ukern))){ # Save iteratively and at the end
      write.csv(total.df, 'copepod_ops_database_wip.csv',
                row.names=F)}
    
    ixp = ixp+1
  }
  
  ## Output a clean dataset
  # Need to reshape the data.frame like it was before, one OPS = one row
  # Easier to work with later
  common.names = c("species", "stage", "phylum", "pred.esd", "pred.esd.sd", "prey.esd.min", "prey.esd.max")
  
  if(test.type %in% c('Q10', 'BS', 'Q10BC')){ # With resamples
    n1 = c("ops1", "rmse")
    n2 = c("ops2", "rmse")
    
    # Split the measurements according to bimodality and homogenize the dataset
    col1 = paste(n1, rep(c('m', 'sd', qt.probs*100), each=length(n1)), sep='.')
    col2 = paste(n2, rep(c('m', 'sd', qt.probs*100), each=length(n2)), sep='.')
    
    # Unique names for the dataset
    unique.names = c(common.names, 
                     paste(gsub('1', '', n1), 
                           rep(c('', '.sd', paste('.', qt.probs*100, sep='')), 
                               each=length(n1)), sep='') )
  }else{ # No resamples
    n1 = c("ops1", "rmse")
    n2 = c("ops2", "rmse")
    
    col1=n1
    col2=n2
    unique.names = c(common.names, gsub('1', '', n1))
  }
  
  # Append the two subsets
  total1 = setNames( total.df[ c(common.names, col1) ], unique.names )
  total2 = setNames( total.df[ c(common.names, col2) ], unique.names )
  
  sf = rbind(total1, total2)
  sf = sf[ which(!is.na(sf$ops)), ] # Remove useless lines (if no bimodality, OPS2==NA)
  
  return(sf)
}

#ops.dt = bs.test(pref) # For extracting the OPS without doing the Monte Carlo simulation
ops.dt = bs.test(pref, nDraws=3, test.type='Q10') # Monte Carlo on Q10 correction

# Append the data by Hansen et al (1994)
names.to.fill = names(ops.dt)[ !names(ops.dt)%in% names(hansen.data)]
hansen.data[names.to.fill] = NA
hansen.data = hansen.data[ names(hansen.data)%in%names(ops.dt) ]

ops.dt = rbind(ops.dt, hansen.data)

## Plot the OPS with the standard deviation from the Monte Carlo simulation
plot.error.bars = function(mean.x.dt, sd.x.dt, mean.y.dt, sd.y.dt, log=F, col.error=rgb(0.5, 0.5, 0.5, alpha = 0.5)){
  
  #col.error = rgb(0.5, 0.5, 0.5, alpha = 0.5)
  
  # Clean the vectors from NAs
  sd.x.dt[ is.na(sd.x.dt) | is.nan(sd.x.dt) ] = 0
  sd.y.dt[ is.na(sd.y.dt) | is.nan(sd.y.dt) ] = 0
  
  for(i in 1:length(mean.x.dt)){
    if(sd.x.dt[i] > 0 | sd.y.dt[i] > 0){
      sd.vec.x = c( mean.x.dt[i] - sd.x.dt[i], mean.x.dt[i] + sd.x.dt[i] ) # Plot the mu +/- sd
      sd.vec.y = c( mean.y.dt[i] - sd.y.dt[i], mean.y.dt[i] + sd.y.dt[i] )
      
      if(log){ # Put in log scale, if necessary
        sd.vec.x = log(sd.vec.x)
        sd.vec.y = log(sd.vec.y)
        mean.x.dt[i] = log(mean.x.dt[i])
        mean.y.dt[i] = log(mean.y.dt[i])
      }
      
      # Plot the error bars
      Arrows(sd.vec.x[1], mean.y.dt[i], sd.vec.x[2], mean.y.dt[i], lty = 1, col = col.error,
             lwd = 1, arr.type = "T", code = 3, arr.length = 0.2)
      Arrows(mean.x.dt[i], sd.vec.y[1], mean.x.dt[i], sd.vec.y[2], lty = 1, col = col.error,
             lwd = 1, arr.type = "T", code = 3, arr.length = 0.2)
    }
  }
}

log10ticks = function(side.ax, tck = 0.02, log = F, axis.lab = T, ...){
  ticks = (-100:100)  # Tick positions
  
  if(axis.lab){
    axis.lab = parse(text = paste0("10^", ticks))
  }
  
  if(log){  # Add the axis ticks
    axis(side.ax, at = log(10^ticks), labels = axis.lab, tck = tck, ...)
  }else{
    axis(side.ax, at = 10^ticks, labels = axis.lab, tck = tck, ...)
  }
  
  # Sub-ticks
  for( xi in 1:9 ){
    if(log){
      axis(side.ax, at = log(10^ticks*xi), labels = F, tck = tck * 0.5, ...)
    }else{
      axis(side.ax, at = 10^ticks*xi, labels = F, tck = tck * 0.5, ...)
    }
  }
} # Plot log10 ticks

x11(height=5, width=8)
par(mgp=c(3, 0.5, 0), cex=1.5, mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04)

plot(ops.dt$pred.esd, ops.dt$ops, xlim=exp( c(4, 8.5) ),
     type='n', lwd=2.5, ann=F, xaxt='n', yaxt='n', log='xy')
plot.error.bars(ops.dt$pred.esd, ops.dt$pred.esd.sd, 
                ops.dt$ops, ops.dt$ops.sd,
                col.error = 'gray30')
points(ops.dt$pred.esd, ops.dt$ops, pch=19, col='gray60')
log10ticks(1, tck=0.04)
log10ticks(2, tck=0.04, las=1)

mtext(side=2, expression('OPS, ESD' *mu *'m'), line=2.1, cex=1.4)
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.6, cex=1.5)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/OPS_error_MC.pdf')
dev.off()

## Save a clean output 
write.csv( ops.dt[ c("species", "phylum", "stage", "pred.esd", "pred.esd.sd", "ops",
                     "prey.esd.min", "prey.esd.max") ], 
           file = "copepod_ops_database.csv", row.names = F )
