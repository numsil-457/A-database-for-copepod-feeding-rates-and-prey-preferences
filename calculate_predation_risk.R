#### Script to calculate the predation risk of copepods based on body size
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Change path to R script directory

library(plyr)
library(dplyr)

source("get_data/get_other_compilations.R")    # Download datasets from Brun et al. (2016) and Pata and Hunt (2023)
source("get_data/get_feeding_behavior_data.R") 
source("get_data/get_zoo_time_series.R") # Extract copepod time series with size > cop.ts

modb = read.csv(file = "copepod_ops_database.csv", header = T) # Copepod OPS dataset
lesd.seq = seq(3.5, 7, length.out=100) # Function of predation risk per size

# Calculate the Imax following the method of Wirtz, 2013 (http://academic.oup.com/plankt/article/35/1/33/1517918/How-fast-can-plankton-feed-Maximum-ingestion-rate)
imax = function(lesd, lops, a0=0.2, v=1.2, n){
  a = n * a0 * (lesd -(-1)**n * lops)
  
  a.lim=2
  a[which(a < 0.0)] = 0
  a[which(a > a.lim)] = a.lim
  
  # Eq.4 in Wirtz JPR 2013
  Imax = 6 * v * exp( (a + (a.lim - a) * lops + (a - a.lim-1) * lesd) ) *24
  
  return( Imax )
}

# Calculate the OPS following the method of Garcia-Oliva and Wirtz, 2025 (https://www.nature.com/articles/s41559-025-02647-1)
OPS.model = function(lesd, s, m){
  PFG.list = log( c(50, 1e2, 1e3, 1e4, 1e5, 1e6) )
  
  # PFGi = PFG.list[ which.min( abs(PFG.list - lesd) ) ] # Corresponding size group
  PFGi = PFG.list[ unlist(sapply(lesd, function(x) which.min(abs(x-PFG.list)))) ]
  
  Stiff = 0.05 * exp(0.26 * PFGi) # Stiffness
  
  mf = 0.12 - 0.23 * PFGi # Ref. feeding mode
  
  smax = 1.9e-2 * PFGi**2 # Max specialization
  
  sf = s * smax # Scaled specialization
  
  lops.mod = -1.73 + exp(-sf**2) * (lesd - PFGi) + PFGi -0.011 * (lesd - (PFGi -3))**2 + m + mf + sf/sqrt(Stiff)
  
  return(lops.mod)
}

## Add the OPS to the gelatinous zooplankton
gel.ts = gel.ts[ which(!is.na(gel.ts$mean.esd)), ]
gel.ts$lops = OPS.model( log(gel.ts$mean.esd), -0.22, 2.45 ) # parameter values were derived using the datasets used by Garcia-Oliva and Wirtz, 2025
gel.ts$imax = imax(log(gel.ts$mean.esd), gel.ts$lops, n=2)

## Extract the fish time series
# fish.juv = read.csv(file = "fish_larvae_wadden_sea.csv", header = T) # Data from Maathuis et al, 2025 (https://linkinghub.elsevier.com/retrieve/pii/S0272771424004311)

# Add the OPS

# Imax from the literature ?

## Add the OPS to the copepod dataset
cop.ts$ops = NA
for(si in unique(cop.ts$Species)){ # Check for the same species first
  ops.sdi = modb[which(modb$species == si),]
  
  if(nrow(ops.sdi)<1){ # Search larger, to the genera level
    si = strsplit(si, ' ')[[1]][1]
    ops.sdi = modb[grep(si, modb$species, fixed=T),]
  }
  
  if(nrow(ops.sdi)>1){
    OPSi = ops.sdi[ unlist(sapply(cop.ts$mean.esd[cop.ts$Species == si], 
                                  function(x) which.min(abs(x-ops.sdi$pred.esd)))), ]
    
    cop.ts$ops[cop.ts$Species == si] = OPSi$ops
  }
}
# length(which(is.na(cop.ts$ops)))/nrow(cop.ts) # 30 % of unknown OPS
# cop.b = add_feeding_trait(cop.ts, feeding.behavior)
# unique(cop.b$detail) # No cruiser detected
# cop.b = cop.b[ which( !is.na(cop.b$ac) ), ] 

# Imax after the model in main.R
imax.metab = function(params, params.activity.1, params.activity.2, params.scale.1, params.scale.2, lesd, lops, n){
  
  activity = function(lesd, params){
    k_a=params[1]; a_f=params[2]; a_shift=params[3]
    
    a = ( a_f + exp(-k_a*(lesd-a_shift)) ) / ( 1 + exp(-k_a*(lesd-a_shift)) )
    
    return( a )
  }
  
  a0 = params[1] # See Wirtz, 2013
  v  = params[2]

    a = n * a0 * (lesd -(-1)**n * lops)
  
  a.lim=2
  a[which(a < 0.0)] = 0
  a[which(a > a.lim)] = a.lim
  
  a.scale = activity(lesd, params.scale.1) * activity(lesd, params.scale.2)
  act  = activity(lesd, params.activity.1) * activity(lesd, params.activity.2)
  vdig = v * act
  
  # Eq.4 in Wirtz JPR 2013 with additional factor e^8: 172h-1*24*e^-8 = 1.4 d-1
  Imax = 6 * vdig * exp( (a + (a.lim - a) * lops + (a - a.lim-1) * lesd) * a.scale ) *24
  
  return( Imax )
}
cop.ts$imax = imax.metab(c(0.24, 1.9/24), 
                         c(2.5, 3.9, 5.5), c(2.5, 3.9, 6.1), 
                         c(2.5, 1.1, 5.5), c(2.5, 1.1, 6.4),
                         log(cop.ts$mean.esd), log(cop.ts$ops), 1)
# Parameter values are taken from Table 1 in the manuscript

## Calculate the overlap of gelatinous, fish, and carnivorous copepods on the copepods size distribution

pred.risk = function(obs.dat, xseq = lesd.seq){
  preference = function(l_prey, l_ops, s=1.5){ # Prey preference model after Wirtz, 2014 (http://www.int-res.com/abstracts/meps/v507/p81-94/)
    exp(-s*(l_prey - l_ops)**2)
  }
  
  pref      = sapply(xseq, preference, l_ops=log(obs.dat$ops) )
  pred.rate = sweep(pref, 1, obs.dat$imax, '*')             # Weighted by imax
  pred.rate = sweep(pred.rate, 1, obs.dat$C..mg.m..3., '*') # Weighted by predator biomass
  
  # Calculate the mean predation risk per event
  obs.dat$UE = paste(obs.dat$Event, obs.dat$date)
  colnames(pred.rate) = xseq
  pred.rate = as.data.frame(pred.rate) # Prepare for merging with obs.dat

  # Merge and calculate the mean predation pressure per unique event
  obs.dat[as.character(xseq)] = pred.rate[as.character(xseq)]

  pred.dat = obs.dat[c('UE', 'C..mg.m..3.', as.character(xseq))] %>% 
    group_by(UE, .add=T) %>% summarise_all('sum', na.rm=T) # Total predation rate and biomass per event
  
  # Plot to see the output
  plot(xseq, apply( pred.dat[as.character(xseq)], 2, 'mean', na.rm=T), 
       ylab='predation risk, micro gC.l-1.d-1',
       xlab='prey log ESD', type='l') # Mean
  
  return(pred.dat)
}

cop.pred = pred.risk(cop.ts)
gel.ts$ops = exp(gel.ts$lops)
gel.pred = pred.risk(gel.ts)

# Merge the datasets per event date, and calculate a total predation risk
pressure_cols = setdiff(names(cop.pred), "UE") # columns to sum
combined      = inner_join(cop.pred, gel.pred, by = "UE", suffix = c("_1", "_2"))

for (col in pressure_cols) {
  combined[[col]] = combined[[paste0(col, "_1")]] + combined[[paste0(col, "_2")]]
}

pred.dat = combined %>% select(UE, all_of(pressure_cols))

plot(lesd.seq, apply( pred.dat[as.character(lesd.seq)], 2, 'mean', na.rm=T), 
     ylab='potential predation, micro gC.l-1.d-1',
     xlab='prey log ESD', type='l') # Mean

# Calculate a mean and quantile of the potential predation risk
pred.mean = exp(apply( log(pred.dat[as.character(lesd.seq)]), 2, 'mean', na.rm=T))
pred.qt   = apply( pred.dat[as.character(lesd.seq)], 2, 'quantile', probs=c(0.1, 0.9))

## Complete plot - predation risk + prey encounter
x11(height=6, width=7)
par(mgp=c(3, 0.3, 0), cex=1.5, mar=c(2.5, 4, 1.5, 0.5), tck=0.02, xpd=F)

plot(lesd.seq, lesd.seq, ylim=log(range(pred.qt)), type='n', ann=F, 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i') # Empty frame

polygon(x=c(lesd.seq, rev(lesd.seq)),
        y=log( c(pred.qt[1,], rev(pred.qt[2,])) ),
        col='gray80', border=NA)
lines(lesd.seq, log(pred.mean), lwd=2)

mtext(side=2, expression('predation risk, ' *mu*'gC.l'^-1*'.d'^-1),  cex=1.5, line=2.5)
mtext(side=1, expression('copepod body size, ESD ' *mu*'m'), cex=1.5, line=1.5)

axis(1, at=log(c(50, 100, 200, 500, 1e3, 2e3)), labels=c(50, 100, 200, 500, 1e3, 2e3))
axis(2, at=log( c(10, 1, 0.1, 0.01, 0.001, 1e-6) ), labels=c(10, 1, 0.1, 0.01, 0.001, 1e-6), las=1)

#mtext(side=3, '(A)', cex=1.5, adj=0.0, font=2)
#mtext(side=3, 'predation pressure in the EFWS', cex=1.5, adj=0.06)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/passive_trade-offs.pdf')
dev.off()
