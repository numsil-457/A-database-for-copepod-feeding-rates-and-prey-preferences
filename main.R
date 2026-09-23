#### Script to fit the metabolic activity model to copepod grazing rates
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Change path to R script directory

library(plyr)
library(dplyr)
library(nloptr)
library(mixtools)
library(parallel)
library(gdata)
# library(readxl)
library(shape)
library(GA)
library(KernSmooth) # Function dpih for binwidths of aggregated data based on an extended rule of thumb by Wand (1995)

## Load other source files
source("get_data/get_other_compilations.R")    # Download datasets from Brun et al. (2016) and Pata and Hunt (2023)
source("get_data/get_respiration_data.R")      # Extract the respiration rates from Brun et al. (2016)
source("get_data/get_feeding_behavior_data.R") # Extract the feeding behavior information

## Count the species in the resp.data 
resp.data$sz.gp = 's'
resp.data$sz.gp[resp.data$esd >=600] = 'l'

resp.count = ddply(resp.data, .(sz.gp,species), summarize,
                   n = length(species))
resp.count[ order(resp.count$n), ] # Species

resp.count = ddply(resp.data, .(sz.gp,reference), summarize,
                   n = length(reference))
resp.count[ order(resp.count$n), ] # Reference

## Colors per feeding mode - @TODO change names
col.filter = rgb(0.2, 0.3, 0.9, alpha = 0.9)
col.ambush = rgb(0.9, 0.1, 0.2, alpha = 0.9)
col.switcher = 'gray50'

## Read the OPS dataset
modb = read.csv(file = "copepod_ops_database.csv", header = T)
# modb = read.csv(file = "modb.csv", header = T)
modb$lops = log(modb$ops)
modb$lesd = log(modb$pred.esd)

# library(mclust)
# plot(Mclust(modb$lops), what='uncertainty')

## Split the OPS values into groups
# Fit a gaussian mixture model with the AIC criteria to find the specialization groups, using the s
wmult = makemultdata(modb$lops, cuts = seq(1, 7, by=0.5)) # Prepared subsets of data to fit mixed gaussians on
multmixmodel.sel(wmult, comps=1:length(seq(1, 7, by=0.5)), epsilon=1e-3) # AIC selects 2 kernels
w1 = normalmixEM(modb$lops, lambda = 0.5, mu = c(1,7), arbvar=F)
#w1 = normalmixEM(modb$lops, lambda=0.5, arbmean=T, arbvar=F)

# Choose the group by using the posterior probabilities of observation
g=rep(NA, nrow(modb))
g[w1$posterior[,1] < w1$posterior[,2]] = 2
g[w1$posterior[,1] >= w1$posterior[,2]] = 1

modb$gp.k = g
plot(modb$lesd, modb$lops, col = modb$gp.k, pch=19)

## Read the Imax data

# Plot the raw Imax dataset
pref = read.csv(file = "copepod_feeding_rates.csv")

## Select the max Ingestion rate per study and species/stage; a bit less evolved than the OPS detection
pref$ind.row = row(pref)[,1]
pref.max = ddply(pref[which(!is.na(pref$Imax.at.15.degreeC..mugC.mugC.1.h.1.)),], 
                 .(species, stage, primary.reference), summarize,
                 i.max = ind.row[ which.max(Imax.at.15.degreeC..mugC.mugC.1.h.1.) ] )

pref.max = pref[pref.max$i.max,]
pref.max = add_feeding_trait(pref.max, feeding.behavior)

## Give an OPS group to Imax

# Give a 'Low' or 'High' criteria for the groups
gp.mops = ddply(modb, .(gp.k), summarize, mops = mean(lops, na.rm=T))
# gp.mops$group = c('low', 'high')[c( which.min(gp.mops$mops), which.max(gp.mops$mops) )]
# 
# modb = merge(modb, gp.mops[c('gp', 'group')], by='gp')

# Based on distance of prey_esd to the mean OPS
pref.max$gp.k = NA
dist.high = abs( log(pref.max$prey.esd) - gp.mops$mops[gp.mops$gp.k==2] ) # Or w1$mu?
dist.low  = abs( log(pref.max$prey.esd) - gp.mops$mops[gp.mops$gp.k==1] )

pref.max$gp.k[dist.high < dist.low]  = 2
pref.max$gp.k[dist.high >= dist.low] = 1

plot(log(pref.max$pred.esd), log(pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.), col = pref.max$gp.k, pch=19)

## Adding the feeding behavior traits to the other datasets
modb      = add_feeding_trait(modb, feeding.behavior)
resp.data = add_feeding_trait(resp.data, feeding.behavior)

## Apply the feeding mode for Acartia and Centropages (active or passive depend on study)
# modb      = add_foraging_switchers(modb)
# pref.max  = add_foraging_switchers(pref.max)
#resp.data = add_foraging_switchers(resp.data) # Cannot say for this, not feeding experiments

## Function to add the plot coding (pch, color) to the dataset
add_graphics_trait = function(dt, is.group.ops=F){
  
  # if(is.group.ops){
  #   dt$colo                   = col.filter
  #   dt$colo[dt$group=='high'] = col.ambush
  # }
  # 
  # dt$pch = 18 # Other
  # dt$pch[dt$ac=='P'] = 1
  # dt$pch[dt$ac=='A'] = 19

  dt$pch = 19 # Other
  dt$colo = 'gray50'
  dt$colo[dt$ac=='P'] = col.ambush
  dt$colo[dt$ac=='S'] = 'darkorange'
  dt$colo[dt$ac=='A'] = col.filter
  
  # + feeding mode detail
  dt$colo.fm = 'gray50'
  dt$colo.fm[dt$detail=='mixed'] = 'darkorange'
  dt$colo.fm[dt$detail=='filter'] = 'blue'
  dt$colo.fm[dt$detail=='ambusher'] = 'red'
  dt$colo.fm[dt$detail=='cruiser'] = 'brown'
  
  return(dt)
}

modb      = add_graphics_trait(modb,     is.group.ops=T)
pref.max  = add_graphics_trait(pref.max, is.group.ops=T)
resp.data = add_graphics_trait(resp.data)

pref.max$lesd  = log(pref.max$pred.esd)
resp.data$lesd = log(resp.data$esd)

plot(modb$lesd, modb$lops, pch=19, col=modb$colo)
plot(pref.max$pred.esd, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., pch=19, 
     col=pref.max$colo, log='xy')
plot(resp.data$esd, resp.data$r, pch=19, col=resp.data$colo, log='xy')

add_OPS_colo = function(dt, is.resp.dt=F, ops.dt=modb){ # Add colors related to the OPS group
  
  if(is.resp.dt){ # Need to check the OPS group in the respiration dataset
    dt$gp.k = NA # Initialize
    
    dt.stage = dt$stage
    dt.stage[dt.stage %in% c('F', 'M')] = 'A' # Adults
    dt.species = gsub("\\s*\\([^\\)]+\\)","",dt$species) # Remove parentheses to compare
    
    for(li in 1:nrow(dt)){
      di = dt[li,]
      
      # print(li)
      # Search for same species and closest stage or size
      ind.species = grep(dt.species[li], ops.dt$species)
      
      if(length(ind.species)>0){
        ind.stage = intersect(ind.species, which(dt.stage[li] == ops.dt$stage))
        
        if(length(ind.stage)<1){
          ind.stage = which.min(abs(dt$lesd[li] - ops.dt$lesd[ind.species]))
          if(length(ind.stage)>0){dt$gp.k[li] = ops.dt$gp.k[ind.species][ind.stage]}
          
        }else{
          dt$gp.k[li] = max(ops.dt$gp.k[ind.stage])
        }
      }
    }
  }
  
  dt$colo = 'gray50'
  dt$colo[dt$gp.k==1] = rgb(0.6, 0.72, 0.3)
  dt$colo[dt$gp.k==2] = rgb(0.3, 0.4, 0.8)
  # dt$colo[dt$gp.k==4] = rgb(0.8, 0.23, 0.2)
  
  return(dt)
}

modb = add_OPS_colo(modb)
plot(modb$lesd, modb$lops, pch=19, col=modb$colo)

pref.max = add_OPS_colo(pref.max)
plot(pref.max$pred.esd, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., pch=19, 
     col=pref.max$colo, log='xy')

resp.data = add_OPS_colo(resp.data, is.resp.dt=T)
plot(resp.data$esd, resp.data$r, pch=19, col=resp.data$colo, log='xy')

## Preparing the datasets and convert rates to .d-1 for better plotting
pref.max$imax = pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1. * 24
resp.data$r   = resp.data$r * 24

## Functions to fit a new specialization model
activity = function(lesd, params){
  k_a=params[1]; a_f=params[2]; a_shift=params[3]
  
  a = ( a_f + exp(-k_a*(lesd-a_shift)) ) / ( 1 + exp(-k_a*(lesd-a_shift)) )
  
  return( a )
}

ops.specialization = function(lesd, params, params.activity.1=NULL, params.activity.2=NULL, linear=F){
  specialization=params[1]; m=params[2]

  if(!linear){ # Bifurcation
    if(is.null(params.activity.2)){
      a = activity(lesd, params.activity.1)
    }else{
      a = activity(lesd, params.activity.1) * activity(lesd, params.activity.2)
    }
  }else{ # Feeding mode constant
    a = 1
  }
  
  lops.estim = exp(-specialization**2) * lesd + m * a
    
  return(lops.estim)
}

imax.new = function(params, params.activity.1, params.activity.2, params.scale.1, params.scale.2, lesd, lops, n){
  a0 = params[1] # See Wirtz, 2013
  v = params[2]
  
  #lops = a1 * lops
  a = n * a0 * (lesd -(-1)**n * lops)
  
  a.lim=2
  a[which(a < 0.0)] = 0
  a[which(a > a.lim)] = a.lim
  
  a.scale = activity(lesd, params.scale.1) * activity(lesd, params.scale.2)
  act  = activity(lesd, params.activity.1) * activity(lesd, params.activity.2)
  vdig = v * act
  
  # Eq.4 in Wirtz JPR 2013 with additional factor e^8: 172h-1*24*e^-8 = 1.4 d-1
  Imax = 6 * vdig * exp( (a + (a.lim - a) * lops + (a - a.lim-1) * lesd) * a.scale ) *24
  
  return( list(Imax, a) )
}

resp.imax = function(params, params.activity.1, params.activity.2, params.scale.1, params.scale.2, lesd){
  gamma = params[1]
  r0    = params[2]
  a_shift = params[3]

  a.scale = activity(lesd, params.scale.1) * activity(lesd, params.scale.2)
  a = activity(lesd, params.activity.1) * activity(lesd, params.activity.2)
  resp.mod = r0 * a * exp(gamma*lesd * a.scale) *24
  return(resp.mod)
}

## Fit the OPS, Imax, and Resp with one framework
fit.all = function(lesd.ops, y.ops, ops.group, lesd.imax, y.imax, imax.group, lesd.resp, y.resp, lb, ub, suggestpar, n=1){
  
  RMSE = function(params, lesd.ops, y.ops, lesd.imax, y.imax, lesd.resp, y.resp, n0){
    params = setNames(params, names.par)

    # Fix some paremeter values?
    # params[c('ac_1', 'ac_2')] = c(5.5, 6.5)
    
    params.ops  = params[c('s', 'm')]
    params.imax = params[c('b0', 'vdig')]
    params.resp = params[c('gamma', 'r0', 'af_scale')]
    
    #ops.mod  = ops.specialization(lesd.ops, params.ops) # No bifurcation
    params.a.ops.1 = params[c('ka', 'af_ops', 'ac_1')]
    params.a.ops.2 = params[c('ka', 'af_ops', 'ac_2')]
    
    # Special rule for the shift
    params.a.ops.1['af_ops'] = params.a.ops.1['af_ops'] * params.a.ops.1['ac_1']
    params.a.ops.2['af_ops'] = params.a.ops.2['af_ops'] * params.a.ops.2['ac_2']
    
    # ind.ops.low  = which(lesd.ops < params['ac_1'] | ops.group==1) # All points before the shift of in the 'low' group
    # ind.ops.mid  = which(lesd.ops < params['ac_1'] | ops.group==2) 
    # ind.ops.high = which(lesd.ops < params['ac_2'] | ops.group==3) 
    ind.ops.low  = which(ops.group==1) # All points before the shift of in the 'low' group
    ind.ops.mid  = which(ops.group==2) 
    # ind.ops.high = which(ops.group==4) 
    
    ops.mod.low  = ops.specialization(lesd.ops[ind.ops.low], params.ops, linear=T)
    ops.mod.mid  = ops.specialization(lesd.ops[ind.ops.mid], params.ops, params.a.ops.1)
    # ops.mod.high = ops.specialization(lesd.ops[ind.ops.high], params.ops, params.a.ops.1, params.a.ops.2)
    
    
    params.a.metab.1 = params[c('ka', 'af_metab', 'ac_1')]
    params.a.metab.2 = params[c('ka', 'af_metab', 'ac_2')]
    
    params.a.scale.1 = params[c('ka', 'af_scale', 'ac_1')]
    params.a.scale.2 = params[c('ka', 'af_scale', 'ac_2')]
    
    # imax.mod = imax.new(params.imax, lesd.imax, 
    #                     ops.specialization(lesd.imax, params.ops), n0)[[1]]
    # ind.imax.low  = which(lesd.imax < params['ac_1'] | imax.group==1) # All points before the shift of in the 'low' group
    # ind.imax.mid  = which(lesd.imax < params['ac_1'] | imax.group==2) 
    # ind.imax.high = which(lesd.imax < params['ac_2'] | imax.group==3) 
    ind.imax.low  = which(imax.group==1) # All points before the shift of in the 'low' group
    ind.imax.mid  = which(imax.group==2) 
    # ind.imax.high = which(imax.group==4) 
    
    imax.mod.low = imax.new(params.imax, params.a.metab.1, params.a.metab.2,
                            params.a.scale.1, params.a.scale.2, lesd.imax[ind.imax.low], 
                            ops.specialization(lesd.imax[ind.imax.low], params.ops, linear = T), n0)[[1]]
    imax.mod.mid = imax.new(params.imax, params.a.metab.1, params.a.metab.2, 
                            params.a.scale.1, params.a.scale.2, lesd.imax[ind.imax.mid], 
                             ops.specialization(lesd.imax[ind.imax.mid], params.ops, params.a.ops.1), n0)[[1]]
    # imax.mod.high = imax.new(params.imax, params.a.metab.1, params.a.metab.2, 
    #                          params.a.scale.1, params.a.scale.2, lesd.imax[ind.imax.high], 
    #                          ops.specialization(lesd.imax[ind.imax.high], params.ops, params.a.ops.1, params.a.ops.2), n0)[[1]]
    
    resp.mod = resp.imax(params.resp, params.a.metab.1, 
                         params.a.metab.2, params.a.scale.1, params.a.scale.2, lesd.resp)
    
    rmse.calculate = function(x, y, normalize=T, log.var=T){
      #if(sd.use) sd.y = sd(y, na.rm=T) else sd.y = 1
      if(log.var){ x=log(x); y=log(y) }
      if(normalize) sd.x = sd(x, na.rm=T) else sd.x = 1
      
      #sqrt( sum( (y - x)**2, na.rm=T ) / sd.y) / length( na.omit(x) )
      #sum( abs(y - x), na.rm=T ) / sd.y / length( na.omit(x) ) 
      sum( abs(y - x), na.rm=T ) / length( na.omit(x) ) / sd.x   # NMAE
      #sqrt( sum( (y - x)**2, na.rm=T ) / length( na.omit(x) ) ) # RMSD
    }

    #rmse.ops = rmse.calculate(y.ops, ops.mod) # No bifurcation
    
    # Choose the closest branch, can be sensitive to fit
    # rmse.ops.high = rmse.calculate(y.ops, ops.mod.high)
    # rmse.ops.low = rmse.calculate(y.ops, ops.mod.low)
    # rmse.ops = apply(rbind(rmse.ops.high, rmse.ops.low), 2, 'min', na.rm=T)
    
    # Based on group, data partitioning more robust
    # rmse.ops = rmse.calculate(y.ops[c(ind.ops.high, ind.ops.mid, ind.ops.low)], 
    #                           c(ops.mod.high, ops.mod.mid, ops.mod.low),
    #                           log.var=F)
    rmse.ops = rmse.calculate(y.ops[c(ind.ops.mid, ind.ops.low)], 
                              c(ops.mod.mid, ops.mod.low),
                              log.var=F)

    # rmse.imax = rmse.calculate(y.imax[c(ind.imax.high, ind.imax.mid, ind.imax.low)],
    #                            c(imax.mod.high, imax.mod.mid, imax.mod.low)
    #                            )
    rmse.imax = rmse.calculate(y.imax[c(ind.imax.mid, ind.imax.low)],
                               c(imax.mod.mid, imax.mod.low)
                              )
    
    rmse.resp = rmse.calculate(y.resp, resp.mod)
    rmse = ( rmse.ops + rmse.imax + rmse.resp ) /3
    # rmse = (rmse.resp + rmse.imax)/2
    # rmse = rmse.imax
    
    return( rmse )
  } 
  
  # NLOPT
  mod = isres(x0 = suggestpar, fn = RMSE,
              lower = lb, upper = ub, maxeval = 5e5L,
              lesd.ops=lesd.ops, y.ops=y.ops,
              lesd.imax=lesd.imax, y.imax=y.imax,
              lesd.resp=lesd.resp, y.resp=y.resp, n0=n)

  # suggestpar = mod$par
  # lb = suggestpar * 0.9 * ( (1-sign(mod$par))*0.1 + 1 )
  # ub = suggestpar * 1.1 * ( -(1-sign(mod$par))*0.1 + 1 )
  # 
  # mod = lbfgs(x0 = suggestpar, fn = RMSE,
  #             lower = lb, upper = ub,
  #             lesd.ops=lesd.ops, y.ops=y.ops,
  #             lesd.imax=lesd.imax, y.imax=y.imax,
  #             lesd.resp=lesd.resp, y.resp=y.resp, n0=n)
  
  print(mod$value)
  mod.par = mod$par
  
  return( mod.par )
}

# Remove the max extrema
# pref.max = pref.max[ which(pref.max$imax<=0.25*24), ]

# Parameters for the non-linear model fit
names.par = c('af_ops', 'af_metab', 'af_scale',
              'ka', 'ac_1', 'ac_2',
              's', 'm',
              'b0', 'vdig',
              'gamma', 'r0')

lb = setNames( c(0, 0, 0.2,
                 0, 5.5, 6.1,
                 
                 0, 0, 
                 0, 0,
                 -1.5, 0),
                 names.par )

ub = setNames( c(5, 10, 5,
                 60, 5.5, 6.1,
                 
                 10, 10,
                 0.5, 0.4,
                 -0.5, 1), 
                 names.par )

suggestpar = setNames( c(0.2, 0.2, 0.2,
                         2, 5.5, 6.1,
               
                         1, 2,
                         0.2, 0.15,
                         -0.75, 0.3),
                         names.par )

## Visualizing the distribution of observations - skewed to large copepods
fit.dens.plot = function(ops.dt=modb, imax.dt=pref.max, resp.dt=resp.data){
 
  # Refined rule from Hardle et al (2004)
  thumbs.rule = function( esd.data ){
    logesd = rep( esd.data, max( 1, 3-length(esd.data) ) )
    # Avoid error if only one point (if more than 1 point, does not change anything)
    
    R = quantile(logesd, probs=c(0.25, 0.75)) 
    h = 1.06 * min( sd(logesd) , (R[2] - R[1])/1.34 ) * length(logesd) ** (-1/5)
    return(h) } 
  
  # Or per dataset?
  x.obs = c(ops.dt$lesd, imax.dt$lesd, resp.dt$lesd[which(!is.na(resp.dt$lesd))])
  bw.bin = thumbs.rule( x.obs )
  print( paste('bw =', round(bw.bin, 2)) )
  
  x.lim = c(4, 8.5)
  
  # Get the colors
  col.low.ops  = unique(ops.dt$colo[ops.dt$gp.k==1])
  col.mid.ops  = unique(ops.dt$colo[ops.dt$gp.k==2])
  col.uncl.ops = unique(resp.dt$colo[is.na(resp.dt$gp.k)])
  # col.high.ops = unique(ops.dt$colo[ops.dt$gp.k==4])
  
  # Model fit
  par(fig=c(0,1,0.66,1), mgp=c(3, 0.3, 0), cex=1.5, mar=c(2.5, 4, 0.5, 0.1), tck=0.04, xpd=F)
  
  # Density of observations
  plot.dens = function(dat, per.gp=T){
    if(!per.gp){dat$gp.k=1; dat$colo='gray50'}
    
    for(gi in unique(dat$gp.k)){
      ind.gi = which(dat$gp.k==gi & !is.na(dat$lesd))
      
      d = density(dat$lesd[ind.gi], bw=thumbs.rule(dat$lesd[ind.gi]), 
                  from = min(x.lim), to = max(x.lim), n = 200)
      polygon(c(d$x, rev(d$x)), c(rep(0, length(d$x)), rev(d$y)), border=NA, 
              col=adjustcolor(unique(dat$colo[ind.gi]), alpha.f=0.7))
    }
  }
  
  # OPS
  plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
  plot.dens(ops.dt)
  axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), tck=0.04, las=1)
  axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=F, tck=0.04)
  mtext(side=3, adj=0.04, 'A - OPS observations', line=-1.4, cex=1.5, font=2)
  
  text(x=log(100), y=0.5,  'low OPS', adj=0., col=col.low.ops)
  text(x=log(900), y=0.5, 'high OPS', adj=0., col=col.mid.ops)
  
  # Imax
  par(fig=c(0,1,0.33,0.66), new=T)
  
  plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
  plot.dens(imax.dt)
  axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), tck=0.04, las=1)
  axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=F, tck=0.04)
  mtext(side=3, adj=0.04, 'B - Imax observations', line=-1.4, cex=1.5, font=2)
  mtext(side=2, 'distribution frequency', line=2, cex=1.5)
  
  # Respiration
  par(fig=c(0,1,0,0.33), new=T)
  
  plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, xaxt='n', yaxt='n', xaxs='i', yaxs='i')
  plot.dens(resp.dt, per.gp=F)
  axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), tck=0.04, las=1)
  axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=c(100, 200, 500, 1e3, 2e3), tck=0.04)
  mtext(side=3, adj=0.04, 'C - Respiration observations', line=-1.4, cex=1.5, font=2)
  mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.5, cex=1.5)
  
  text(x=log(400), y=0.5,  'unclear', adj=0., col=col.uncl.ops)
}

x11(height=12, width=7)

fit.dens.plot()

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/observations_densities.pdf')
dev.off()

# Calculate aggregated datasets
bins.aggregate = function(x){
  h = dpih(x[!is.na(x)])
  bins = seq(min(x, na.rm=T)-h, max(x, na.rm=T)+h, by=h)
  return(bins)
} # See Wand (1995)

ops.agg = data.frame() ; imax.agg = data.frame() ; resp.agg = data.frame()
ops.bin  = bins.aggregate(modb$lesd)
imax.bin = bins.aggregate(pref.max$lesd)
resp.bin = bins.aggregate(resp.data$lesd)

resp.data$lr = log(resp.data$r)
for( ix in 1:(length(resp.bin)-1) ){
  xi = resp.bin[ix] ; xf = resp.bin[ix+1]
  
  # No groups in resp
  ind.resp = which(resp.data$lesd >= xi & resp.data$lesd < xf )
  
  resp.agg = rbind( resp.agg,
                    setNames(resp.data[ind.resp,c('lesd', 'lr')] %>% summarise_all(list('mean', 'sd'), na.rm=T),
                             c( c('lesd', 'lr'), paste(c('lesd', 'lr'), 'sd', sep='_') ) )
  )
}

pref.max$li = log(pref.max$imax)
for( ix in 1:(length(imax.bin)-1) ){
  xi = imax.bin[ix] ; xf = imax.bin[ix+1]
  
  for(gi in 1:max(pref.max$gp.k)){
    ind.imax = which(pref.max$lesd >= xi & pref.max$lesd < xf & pref.max$gp.k == gi)

    imax.add = setNames(pref.max[ind.imax,c('lesd', 'li')] %>% summarise_all(list('mean', 'sd'), na.rm=T),
                        c( c('lesd', 'li'), paste(c('lesd', 'li'), 'sd', sep='_') ) )
    imax.add$gp.k = gi
    imax.add$colo = unique(pref.max$colo[pref.max$gp.k==gi])
    
    imax.agg = rbind(imax.agg, imax.add)
  }
}

for( ix in 1:(length(ops.bin)-1) ){
  xi = ops.bin[ix] ; xf = ops.bin[ix+1]
  
  for(gi in 1:max(modb$gp.k)){
    ind.ops  = which(modb$lesd  >= xi & modb$lesd  < xf & modb$gp.k  == gi)
    
    ops.add = setNames(modb[ind.ops,c('lesd', 'lops')] %>% summarise_all(list('mean', 'sd'), na.rm=T),
                       c( c('lesd', 'lops'), paste(c('lesd', 'lops'), 'sd', sep='_') ) )
    ops.add$gp.k = gi
    ops.add$colo = unique(modb$colo[modb$gp.k==gi])
    
    ops.agg = rbind(ops.agg, ops.add)
  }
}

# for( ix in 1:(length(bin.seq)-1) ){
#   xi = bin.seq[ix] ; xf = bin.seq[ix+1]
#   
#   # No groups in resp
#   ind.resp = which(resp.dt$lesd >= xi & resp.dt$lesd < xf )
#   
#   resp.local = rbind( resp.local,
#                       setNames(resp.dt[ind.resp,c('lesd', 'lr')] %>% summarise_all(list('mean', 'sd'), na.rm=T),
#                                c( c('lesd', 'lr'), paste(c('lesd', 'lr'), 'sd', sep='_') ) )
#   )
#   
#   for(gi in 1:max(imax.dt$gp.k)){
#     ind.imax = which(imax.dt$lesd >= xi & imax.dt$lesd < xf & imax.dt$gp.k == gi)
#     ind.ops  = which(ops.dt$lesd  >= xi & ops.dt$lesd  < xf & ops.dt$gp.k  == gi)
#     
#     imax.add = setNames(imax.dt[ind.imax,c('lesd', 'li')] %>% summarise_all(list('mean', 'sd'), na.rm=T),
#                         c( c('lesd', 'li'), paste(c('lesd', 'li'), 'sd', sep='_') ) )
#     imax.add$gp.k = gi
#     imax.add$colo = unique(imax.dt$colo[imax.dt$gp.k==gi])
#     
#     ops.add = setNames(ops.dt[ind.ops,c('lesd', 'lops')] %>% summarise_all(list('mean', 'sd'), na.rm=T),
#                        c( c('lesd', 'lops'), paste(c('lesd', 'lops'), 'sd', sep='_') ) )
#     ops.add$gp.k = gi
#     ops.add$colo = unique(ops.dt$colo[ops.dt$gp.k==gi])
#     
#     imax.local = rbind(imax.local, imax.add)
#     ops.local  = rbind(ops.local,  ops.add )
#   }
# }

resp.agg$pch = 19
imax.agg$pch = 19
ops.agg$pch  = 19

## Fitting the non-linear model to aggregated data points
agg.par = fit.all(ops.agg$lesd,  ops.agg$lops,     ops.agg$gp.k,
                  imax.agg$lesd, exp(imax.agg$li), imax.agg$gp.k,
                  resp.agg$lesd, exp(resp.agg$lr),
                  lb, ub, suggestpar)
agg.par = setNames(agg.par, names.par)

## Fitting the non-linear model to all data points
all.par = fit.all(modb$lesd,      modb$lops,     modb$gp.k,
                  pref.max$lesd,  pref.max$imax, pref.max$gp.k,
                  resp.data$lesd, resp.data$r,
                  lb, ub, suggestpar)
all.par = setNames(all.par, names.par)

## Plotting function for the model fit
pref.max$li  = log(pref.max$imax)
resp.data$lr = log(resp.data$r)
fit.plot = function(vec.par, ops.dt, imax.dt, resp.dt, agg=F){
  
  # Model predictions
  params.ops  = vec.par[c('s', 'm')]
  params.imax = vec.par[c('b0', 'vdig')]
  params.resp = vec.par[c('gamma', 'r0', 'af_scale')]
  
  params.a.ops.1 = vec.par[c('ka', 'af_ops', 'ac_1')]
  params.a.ops.2 = vec.par[c('ka', 'af_ops', 'ac_2')]
  
  # Special rule for the OPS shift
  params.a.ops.1['af_ops'] = params.a.ops.1['af_ops'] * params.a.ops.1['ac_1']
  params.a.ops.2['af_ops'] = params.a.ops.2['af_ops'] * params.a.ops.2['ac_2']
  
  params.a.metab.1 = vec.par[c('ka', 'af_metab', 'ac_1')]
  params.a.metab.2 = vec.par[c('ka', 'af_metab', 'ac_2')]
  params.a.scale.1 = vec.par[c('ka', 'af_scale', 'ac_1')]
  params.a.scale.2 = vec.par[c('ka', 'af_scale', 'ac_2')]
  
  # x-axis
  x.seq             = seq(3, 9, length.out=150)
  x.lim             = c(4, 8.5)
  lesd.shift        = vec.par[c('ac_1', 'ac_2')]
  
  ops.seq.low  = ops.specialization(x.seq, params.ops, linear=T)
  ops.seq.mid  = ops.specialization(x.seq, params.ops, params.a.ops.1)
  # ops.seq.high = ops.specialization(x.seq, params.ops, params.a.ops.1, params.a.ops.2)
  
  imax.seq.low = imax.new(params.imax, params.a.metab.1, params.a.metab.2, 
                          params.a.scale.1, params.a.scale.2, x.seq, ops.seq.low, 1)[[1]]
  imax.seq.mid = imax.new(params.imax, params.a.metab.1, params.a.metab.2, 
                          params.a.scale.1, params.a.scale.2, x.seq, ops.seq.mid, 1)[[1]]
  # imax.seq.high = imax.new(params.imax, params.a.metab.1, params.a.metab.2, 
  #                          params.a.scale.1, params.a.scale.2, x.seq, ops.seq.high, 1)[[1]]
  
  resp.seq =  resp.imax(params.resp, params.a.metab.1, params.a.metab.2, 
                        params.a.scale.1, params.a.scale.2, x.seq)
  
  # Functions for a nice plot
  axis.ticks = function(side.ax, tck = 0.04, log = F, axis.lab = T, ...){
    ticks = (-100:100)  # Tick positions
    
    # Sub-ticks
    for( xi in 1:10 ){
      if(log){
        axis(side.ax, at = log(10^ticks*xi), labels = F, tck = tck * 0.5, ...)
      }else{
        axis(side.ax, at = 10^ticks*xi, labels = F, tck = tck * 0.5, ...)
      }
    }
    
    # Main ticks
    for( xi in c(5,10) ){
      if(log){
        axis(side.ax, at = log(10^ticks*xi), 
             labels = 10^ticks*xi, #parse(text = paste0(xi,"^", ticks)), 
             tck = tck, ...)
      }else{
        axis(side.ax, at = 10^ticks*xi, tck = tck, ...)
      }
    }
    
    # if(axis.lab){
    #   axis.lab = parse(text = paste0("10^", ticks))
    # }
    
  } # Plot log10 ticks
  
  # Error bars
  plot.error.bars = function(mean.x.dt, sd.x.dt, mean.y.dt, sd.y.dt, log=F, col.error=rgb(0.5, 0.5, 0.5, alpha = 0.5)){
    
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
  
  # Line with color gradient 
  line.col.gradient = function(x, y, x.shift, col1, col2, ...){
    ind.shift = which.min( abs(x-x.shift) )
    n = length(x)
    
    col.middle = colorRampPalette( c(col1, col2) )(3)[2] # Middle color between blue and red
    
    col.seq.1 = colorRampPalette( c(col1, col.middle) )(ind.shift)
    col.seq.2 = colorRampPalette( c(col.middle, col2) )(n- ind.shift)
    cols = c(col.seq.1, col.seq.2)
    
    segments(x[-n], y[-n], x[-1], y[-1], col = cols, ...)
  }
  
  # Get the colors
  col.low.ops  = unique(ops.dt$colo[ops.dt$gp.k==1])
  col.mid.ops  = unique(ops.dt$colo[ops.dt$gp.k==2])
  # col.uncl.ops = unique(ops.dt$colo[ops.dt$gp.k==3])
  # col.high.ops = unique(ops.dt$colo[ops.dt$gp.k==4])
  
  par(fig=c(0,1,0.65,0.95), mgp=c(3, 0.5, 0), cex=1.5, mar=c(1, 4, 0.1, 0.1), tck=0.04, pch=19, xpd=F)
  
  # OPS
  if(agg){ y.lim=c(1.5, 6.5) }else{ y.lim=c(1.3, 6.5) }
  plot(ops.dt$lesd, ops.dt$lops, type='n', ann=F, ylim=y.lim, xlim=x.lim, xaxt='n', yaxt='n')
  axis(2, at=log(c(5, 12, 30, 60, 120, 240)), labels=c(5, 12, 30, 60, 120, 240), tck=0.04, las=1)
  axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=F, tck=0.04)
  
  if(agg){
    plot.error.bars(ops.dt$lesd, ops.dt$lesd_sd, ops.dt$lops, ops.dt$lops_sd)
  }
  
  lines(x.seq, x.seq + log(0.1), col='gray50', lwd=4, lty=1)
  
  # lines(x.seq, ops.seq.high, col='white', lwd=6)
  # lines(x.seq, ops.seq.high, col=col.high.ops, lwd=6)
  lines(x.seq, ops.seq.mid,  col='white', lwd=6)
  # lines(x.seq, ops.seq.mid,  col=col.mid.ops, lwd=4)
  line.col.gradient(x.seq, ops.seq.mid, lesd.shift[1], col.low.ops, col.mid.ops, lwd=4)
  lines(x.seq, ops.seq.low,  col='white', lwd=6)
  lines(x.seq, ops.seq.low,  col=col.low.ops, lwd=4)
  
  points(ops.dt$lesd, ops.dt$lops, cex=1.4, col='white')
  points(ops.dt$lesd, ops.dt$lops, col=adjustcolor(ops.dt$colo, alpha.f=0.7))
  abline(v = lesd.shift[1], col='gray50', lty=2, lwd=4)
  
  mtext(side=2, expression('OPS ESD, ' *mu *'m'), line=2.5, cex=1.5)
  mtext(side=3, adj=0.02, 'A', line=-1.1, cex=1.5, font=2)

  # text(x=log(1200), y=log(550), 'High OPS',      col=col.high.ops)
  # text(x=log(2200), y=log(80),  'Group unclear', col=col.uncl.ops)
  text(x=log(700), y=log(500), 'high OPS', col=col.mid.ops, adj=0)
  text(x=log(700), y=log(8),   'low OPS',  col=col.low.ops, adj=0)
  text(x=log(800), y=log(50),  'linear scaling', col='gray50', adj=0)
  
  # Imax
  par(fig=c(0,1,0.35,0.65), new=T)
  
  if(agg){ y.lim=log(c(0.1, 4)) }else{ y.lim=c(-6, -0.85)+log(24) }
  plot(imax.dt$lesd, imax.dt$li, type='n', ann=F, 
       ylim=y.lim, xlim=x.lim, xaxt='n', yaxt='n')
  axis.ticks(2, log=T, las=1)
  axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=F, tck=0.04)
  
  if(agg){
    plot.error.bars(imax.dt$lesd, imax.dt$lesd_sd, imax.dt$li, imax.dt$li_sd)
  }
  
  saiz.model = log( exp(0.225)* ( (exp(x.seq)**3*0.523)/(8.3*1e6) )**(0.703-1) )
  lines(x.seq, saiz.model, col='gray50', lwd=4, lty=1)
  
  # lines(x.seq, log(imax.seq.high), col='white',      lwd=6)
  # lines(x.seq, log(imax.seq.high), col=col.high.ops, lwd=4)
  lines(x.seq, log(imax.seq.mid),  col='white',      lwd=6)
  # lines(x.seq, log(imax.seq.mid),  col=col.mid.ops,  lwd=4)
  line.col.gradient(x.seq, log(imax.seq.mid), lesd.shift[1], col.low.ops, col.mid.ops, lwd=4)
  lines(x.seq, log(imax.seq.low),  col='white',      lwd=6)
  lines(x.seq, log(imax.seq.low),  col=col.low.ops,  lwd=4)
  
  points(imax.dt$lesd, imax.dt$li, cex=1.4, col='white')
  points(imax.dt$lesd, imax.dt$li, col=adjustcolor(imax.dt$colo, alpha.f=0.7), lwd=2.5)
  
  abline(v = lesd.shift, col='gray50', lty=2, lwd=4)
  
  mtext(side=2, expression('Imax, d'^-1), line=2.5, cex=1.5)
  mtext(side=3, adj=0.02, 'B', line=-1.1, cex=1.5, font=2)
  
  # Respiration
  par(fig=c(0,1,0.05,0.35), new=T)
  
  if(agg){ y.lim=log(c(0.01, 0.2)) }else{ y.lim=c(-8, -4)+log(24) }
  plot(resp.dt$lesd, resp.dt$lr, type='n', ann=F, 
       ylim=y.lim, xlim=x.lim, xaxt='n', yaxt='n')
  axis.ticks(2, log=T, las=1)
  axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=c(100, 200, 500, 1e3, 2e3), tck=0.04)
  
  if(agg){
    plot.error.bars(resp.dt$lesd, resp.dt$lesd_sd, resp.dt$lr, resp.dt$lr_sd)
  }
  
  lm.resp = lm(lr + 0.75 * lesd ~ 1, data=resp.dt) # Specified slope
  abline(lm.resp$coefficients, -0.75, lty=1, col='gray50', lwd=4)
  
  lines(x.seq, log(resp.seq), lwd=6, col='white')
  lines(x.seq, log(resp.seq), lwd=4, col='black')
  
  points(resp.dt$lesd, resp.dt$lr, cex=1.4, col='white')
  points(resp.dt$lesd, resp.dt$lr, lwd=2.5, col='black')#adjustcolor('gray50', alpha.f=0.7))
  
  # text(x=log(60), y=log(0.03), 'unclear OPS', col='black', adj=0)
  
  abline(v = lesd.shift, col='gray50', lty=2, lwd=4)
  
  mtext(side=2, expression('respiration rate, d'^-1), line=2.5, cex=1.5)
  mtext(side=1, expression('copepod body size ESD, ' *mu *'m'), line=1.8, cex=1.5)
  
  mtext(side=3, adj=0.02, 'C', line=-1.1, cex=1.5, font=2)
  
  # Legend
  par(fig=c(0,1,0.95,1), mar=c(0, 4, 0, 0.1), new=T)
  
  plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, axes=F)
  
  lines(x = c(0, lesd.shift[1])*c(1.01,0.99), y = rep(0, 2), lty = 1, lwd=4, col='gray50')
  lines(x = c(lesd.shift[1], lesd.shift[2])*c(1.01,0.99), y = rep(0, 2), lty = 1, lwd=4, col='gray50')
  lines(x = c(lesd.shift[2], 10)*c(1.01,0.99), y = rep(0, 2), lty = 1, lwd=4, col='gray50')
  
  text(x=(lesd.shift[1] + x.lim[1])/2,      y=0.4, 'Passive',  adj=0.5)
  text(x=(lesd.shift[1] + lesd.shift[2])/2, y=0.4, 'Switcher', adj=0.5)
  text(x=(lesd.shift[2] + x.lim[2])/2,      y=0.4, 'Active',   adj=0.5)
  
  # par(fig=c(0,1,0,1), new=T)
  # 
  # legend('topleft', c('non-linear model', 'classic linear model'),#, 'metabolic activity shift'), 
  #        pch=NA, col=c('black', 'gray50'), lwd=4, lty=1,
  #        bty='n', x.intersp=0.4, inset=c(0., -0.01), horiz=T, xpd=T)
}

## Plot
x11(height=12, width=7)

fit.plot(agg.par, ops.agg, imax.agg, resp.agg, agg=T)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/physio_model_agg.pdf')
dev.off()

x11(height=12, width=7)

fit.plot(all.par, modb, pref.max, resp.data)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/physio_model_all.pdf')
dev.off()

## Save the parameters
write.table(as.data.frame(agg.par), 
            '~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/agg.csv')

write.table(as.data.frame(all.par), 
            '~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/all.csv')

## Plot the ecdf of OPS data
ecdf.low.ops  = ecdf(modb$lesd[modb$gp.k==1])
ecdf.mid.ops  = ecdf(modb$lesd[modb$gp.k==2])
# ecdf.high.ops = ecdf(modb$lesd[modb$gp.k==4])

x.seq = seq(1, 10, by=0.05)
lesd.shift    = agg.par[c('ac_1', 'ac_2')]
min.ops.shift = ecdf.mid.ops(lesd.shift[1])
#mean( c(ecdf.mid.ops(lesd.shift[1]), ecdf.high.ops(lesd.shift[2])) )

x11(height=7, width=7)
par(fig=c(0,1,0,0.9), mgp=c(3, 0.3, 0), cex=1.5, mar=c(4, 4, 0.1, 0.1), tck=0.04, xpd=F)

plot(ecdf.low.ops, ylim=c(0,1), xlim=range(modb$lesd), ann=F, verticals=T, do.points=F,
     xaxt='n', yaxt='n', lwd=3, col=unique(modb$colo[modb$gp.k==1]))
plot(ecdf.mid.ops, add=T, verticals=T, do.points=F, lwd=3,
     col=unique(modb$colo[modb$gp.k==2]))
# plot(ecdf.high.ops, add=T, verticals=T, do.points=F, lwd=3,
#      col=unique(modb$colo[modb$gp.k==4]))

axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), las=1)
axis(2, at=min.ops.shift, labels=round(min.ops.shift, 2), font=2, las=1)
axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=c(100, 200, 500, 1e3, 2e3))
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=2, cex=1.5)
mtext(side=2, 'cumulative distribution function', line=2.5, cex=1.5)

# abline(v = lesd.shift[1], col='gray50', lty=2, lwd=4)
abline(v = lesd.shift[1], col='gray50', lty=2, lwd=4)
abline(h = min.ops.shift, col='gray40', lty=3, lwd=2)

# legend('bottomright', c('low OPS', 'Medium OPS', 'High OPS'),
#        pch=NA, col=c(unique(modb$colo[modb$gp.k==1]), unique(modb$colo[modb$gp.k==2]), unique(modb$colo[modb$gp.k==4])),
#        lwd=4, bty='n', x.intersp=0.4, inset=c(0., 0.25), xpd=T)
legend('bottomright', c('low OPS', 'high OPS'),
       pch=NA, col=c(unique(modb$colo[modb$gp.k==1]), unique(modb$colo[modb$gp.k==2])),
       lwd=4, bty='n', x.intersp=0.4, inset=c(0., 0.25), xpd=T)

# Legend
par(fig=c(0,1,0.9,1), mar=c(0.1,4,0.1,0.1), new=T)

x.lim = range(modb$lesd)
plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, axes=F)

lines(x = c(0, lesd.shift[1])*c(1.01,0.99), y = rep(0, 2), lty = 1, lwd=4, col='gray50')
lines(x = c(lesd.shift[1], 10)*c(1.01,0.99), y = rep(0, 2), lty = 1, lwd=4, col='gray50')
# lines(x = c(lesd.shift[2], 10)*c(1.01,0.99), y = rep(0, 2), lty = 1, lwd=4, col='gray50')

text(x=(lesd.shift[1] + x.lim[1])/2, y=0.3, 'Passive', adj=0.5)
text(x=(lesd.shift[1] + x.lim[2])/2, y=0.3, 'Switcher + Active', adj=0.5)
# text(x=(lesd.shift[2] + x.lim[2])/2, y=0.3, 'Active', adj=0.5)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/ecdf_OPS.pdf')
dev.off()

## Compare the MAEs of the reference models VS metabolic model
pref.max$li  = log(pref.max$imax)
resp.data$lr = log(resp.data$r)
MAE.model = function(vec.par, ops.dt, imax.dt, resp.dt){
  mae = function(y.obs, y.model){
    sum( abs(y.obs-y.model), na.rm=T ) / length( na.omit(y.obs) )
    #sqrt( sum( (y.obs-y.model)**2, na.rm=T ) / length( na.omit(y.obs) ) )
  }
  
  # Model predictions
  params.ops  = vec.par[c('s', 'm')]
  params.imax = vec.par[c('b0', 'vdig')]
  params.resp = vec.par[c('gamma', 'r0', 'af_scale')]
  
  params.a.ops.1 = vec.par[c('ka', 'af_ops', 'ac_1')]
  params.a.ops.2 = vec.par[c('ka', 'af_ops', 'ac_2')]
  
  # Special rule for the OPS shift
  params.a.ops.1['af_ops'] = params.a.ops.1['af_ops'] * params.a.ops.1['ac_1']
  params.a.ops.2['af_ops'] = params.a.ops.2['af_ops'] * params.a.ops.2['ac_2']
  
  params.a.metab.1 = vec.par[c('ka', 'af_metab', 'ac_1')]
  params.a.metab.2 = vec.par[c('ka', 'af_metab', 'ac_2')]
  params.a.scale.1 = vec.par[c('ka', 'af_scale', 'ac_1')]
  params.a.scale.2 = vec.par[c('ka', 'af_scale', 'ac_2')]
  
  # x-axis
  x.seq             = seq(3, 9, length.out=150)
  x.lim             = c(4, 8.5)
  lesd.shift        = vec.par[c('ac_1', 'ac_2')]
  
  ## Calculate the model series
  
  # Respiration
  lm.resp = lm(lr + 0.75 * lesd ~ 1, data=resp.dt) # Specified slope
  resp.mod.l  = resp.dt$lesd*(-0.75) + lm.resp$coefficients
  
  resp.mod.nl = resp.imax(params.resp, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2, resp.dt$lesd)
   
  # OPS
  ops.mod.l  = ops.dt$lesd + log(0.1)
  
  # ops.seq.low  = ops.specialization(ops.dt$lesd[ops.dt$gp.k==1], params.ops, linear=T)
  # ops.seq.mid  = ops.specialization(ops.dt$lesd[ops.dt$gp.k==2], params.ops, params.a.ops.1)
  # # ops.seq.high = ops.specialization(ops.dt$lesd[ops.dt$gp.k==3], params.ops, params.a.ops.1, params.a.ops.2)
  # # ops.mod.nl   = c(ops.seq.low, ops.seq.mid)#, ops.seq.high)
  # ops.mod.nl   = c(ops.specialization(ops.dt$lesd, params.ops, linear=T)[ops.dt$gp.k==1], 
  #                  ops.specialization(ops.dt$lesd, params.ops, params.a.ops.1)[ops.dt$gp.k==2])#, ops.seq.high)
  ops.mod.nl = ops.dt$gp.k
  ops.mod.nl[ops.dt$gp.k==1] = ops.specialization(ops.dt$lesd, params.ops, linear=T)[ops.dt$gp.k==1]
  ops.mod.nl[ops.dt$gp.k==2] = ops.specialization(ops.dt$lesd, params.ops, params.a.ops.1)[ops.dt$gp.k==2]
  
  # Imax
  imax.mod.l = log( exp(0.225)* ( (exp(imax.dt$lesd)**3*0.523)/(8.3*1e6) )**(0.703-1) )
  
  # imax.mod.low = imax.new(params.imax, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2,
  #                         imax.dt$lesd[imax.dt$gp.k==1], 
  #                         ops.specialization(imax.dt$lesd[imax.dt$gp.k==1], params.ops, linear = T), 1)[[1]]
  # imax.mod.mid = imax.new(params.imax, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2,
  #                         imax.dt$lesd[imax.dt$gp.k==2], 
  #                         ops.specialization(imax.dt$lesd[imax.dt$gp.k==2], params.ops, params.a.ops.1), 1)[[1]]
  # # imax.mod.high = imax.new(params.imax, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2,
  # #                          imax.dt$lesd[imax.dt$gp.k==3], 
  # #                          ops.specialization(imax.dt$lesd[imax.dt$gp.k==3], params.ops, params.a.ops.1, params.a.ops.2), 1)[[1]]
  # imax.mod.nl = c(imax.mod.low, imax.mod.mid)#, imax.mod.high)
  
  imax.mod.nl = imax.dt$gp.k
  imax.mod.nl[imax.dt$gp.k==1] = imax.new(params.imax, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2,
                                          imax.dt$lesd[imax.dt$gp.k==1], 
                                          ops.specialization(imax.dt$lesd[imax.dt$gp.k==1], params.ops, linear = T), 1)[[1]]
  imax.mod.nl[imax.dt$gp.k==2] = imax.new(params.imax, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2,
                                          imax.dt$lesd[imax.dt$gp.k==2], 
                                          ops.specialization(imax.dt$lesd[imax.dt$gp.k==2], params.ops, params.a.ops.1), 1)[[1]]
    
  ## Calculate the MAEs - per size group (P/S/A) and all
  size.gp.lb = c(0,  0,                  lesd.shift['ac_1'], lesd.shift['ac_2'])
  size.gp.ub = c(10, lesd.shift['ac_1'], lesd.shift['ac_2'], 10)
  # Total / First group / Second group / Third group
  
  mae.dt   = matrix(NA, nrow=6, ncol=length(size.gp.lb))
  name.col = c()
  
  for(si in 1:length(size.gp.lb)){
    # MAE of respiration
    ind.r = which(resp.dt$lesd >= size.gp.lb[si] & resp.dt$lesd < size.gp.ub[si])
    mae.r.nl = mae(resp.dt$lr[ind.r], log(resp.mod.nl)[ind.r]) # Model
    mae.r.l  = mae(resp.dt$lr[ind.r], resp.mod.l[ind.r])    # Linear regression
    
    # MAE of OPS
    ind.ops = which(ops.dt$lesd >= size.gp.lb[si] & ops.dt$lesd < size.gp.ub[si])
    mae.ops.l  = mae(ops.dt$lops[ind.ops], ops.mod.l[ind.ops]) 
    
    ind.ops.gp = c( which(ops.dt$lesd >= size.gp.lb[si] & ops.dt$lesd < size.gp.ub[si] & ops.dt$gp.k == 1),
                    which(ops.dt$lesd >= size.gp.lb[si] & ops.dt$lesd < size.gp.ub[si] & ops.dt$gp.k == 2)#,
                    # which(ops.dt$lesd >= size.gp.lb[si] & ops.dt$lesd < size.gp.ub[si] & ops.dt$gp.k == 4)
                  )
    mae.ops.nl = mae(ops.dt$lops[ind.ops.gp], ops.mod.nl[ind.ops.gp])
    
    # MAE of Imax
    ind.i = which(imax.dt$lesd >= size.gp.lb[si] & imax.dt$lesd < size.gp.ub[si])
    mae.imax.l = mae(imax.dt$li[ind.i], imax.mod.l[ind.i])
    
    ind.i.gp = c( which(imax.dt$lesd >= size.gp.lb[si] & imax.dt$lesd < size.gp.ub[si] & imax.dt$gp.k == 1),
                  which(imax.dt$lesd >= size.gp.lb[si] & imax.dt$lesd < size.gp.ub[si] & imax.dt$gp.k == 2)#,
                  # which(imax.dt$lesd >= size.gp.lb[si] & imax.dt$lesd < size.gp.ub[si] & imax.dt$gp.k == 4)
                )
    mae.imax.nl = mae(imax.dt$li[ind.i.gp], log( imax.mod.nl[ind.i.gp] ))
    
    # Save
    name.col = c(name.col, paste(size.gp.lb[si], size.gp.ub[si], sep=' / ') )
    mae.dt[,si] = c(mae.r.l, mae.r.nl, mae.ops.l, mae.ops.nl, mae.imax.l, mae.imax.nl)
  }
  
  colnames(mae.dt) = name.col
  rownames(mae.dt) = c('mae.r.l', 'mae.r.nl', 'mae.ops.l', 'mae.ops.nl', 'mae.imax.l', 'mae.imax.nl')
  
  return(mae.dt)
}

mae.agg = MAE.model(agg.par, ops.agg, imax.agg, resp.agg)
mae.all = MAE.model(all.par, modb,    pref.max, resp.data)

write.table(round(mae.agg, 3), 
            '~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/MAE_agg.csv')

write.table(round(mae.all, 3), 
            '~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/MAE_all.csv')

## ANOVA for the slope of respiration rates (size limit at 800 mu m ESD)
resp.data$size_group = 's'
resp.data$size_group[resp.data$lesd >= agg.par['ac_2']] = 'l'

linear_model = lm(lr ~ lesd * size_group, data = resp.data)#[resp.data$ac=='A',])
myanova = anova(linear_model)
qf(1-0.05, myanova[3,1], myanova[4,1]) # Critical F value

coef(lm(lr ~ lesd, resp.data[resp.data$size_group=='l',]))
coef(lm(lr ~ lesd, resp.data[resp.data$size_group=='s',]))

## Compare the classical model formulations and the non-linear model
source("get_data/get_zoo_time_series.R") # Extract a copepod time series with size > cop.ts

# Get parameters
params.ops  = agg.par[c('s', 'm')]
params.imax = agg.par[c('b0', 'vdig')]
params.resp = agg.par[c('gamma', 'r0', 'af_scale')]

params.a.ops.1 = agg.par[c('ka', 'af_ops', 'ac_1')]
params.a.ops.2 = agg.par[c('ka', 'af_ops', 'ac_2')]

# Special rule for the OPS shift
params.a.ops.1['af_ops'] = params.a.ops.1['af_ops'] * params.a.ops.1['ac_1']
params.a.ops.2['af_ops'] = params.a.ops.2['af_ops'] * params.a.ops.2['ac_2']

params.a.metab.1 = agg.par[c('ka', 'af_metab', 'ac_1')]
params.a.metab.2 = agg.par[c('ka', 'af_metab', 'ac_2')]
params.a.scale.1 = agg.par[c('ka', 'af_scale', 'ac_1')]
params.a.scale.2 = agg.par[c('ka', 'af_scale', 'ac_2')]

# Assign foraging behavior to the copepods time series
# Differenciate cruisers (high OPS) and current feeders (low OPS)
cop.ts$species = cop.ts$Species
cop.b = add_feeding_trait(cop.ts, feeding.behavior)
# unique(cop.b$detail) # No cruiser detected
cop.b = cop.b[ which( !is.na(cop.b$ac) ), ] 

# Calculate the grazing preferences
# Respiration
lm.resp = lm(lr + 0.75 * lesd ~ 1, data=resp.agg) # Specified slope
cop.b$r.linear = exp( log(cop.b$mean.esd)*(-0.75) + lm.resp$coefficients )
cop.b$r.metab  = resp.imax(params.resp, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2, log(cop.b$mean.esd))

# OPS
cop.b$ops.linear = log(cop.b$mean.esd) + log(0.1)
cop.b$ops.metab  = ops.specialization(log(cop.b$mean.esd), params.ops, linear=T) # No cruisers so no high OPS

# Imax
cop.b$imax.linear = exp(0.225)* ( (cop.b$mean.esd**3*0.523)/(8.3*1e6) )**(0.703-1)
cop.b$imax.metab  = imax.new(params.imax, params.a.metab.1, params.a.metab.2, params.a.scale.1, params.a.scale.2,
                             log(cop.b$mean.esd), cop.b$ops.metab, 1)[[1]]

# Calculate the biomass-weighted-means of the rates
# Statistical analysis functions
time.average = function(station, date, biomass, size, ops, imax, resp){
  # Reconstruct the samples, as they are exploded into multiple lines (one per species observed)
  usi = paste( station, "_", date, sep="" )
  data = data.frame( usi = usi, station = station,
                     biomass = biomass, size = size,
                     ops = ops, imax = imax, resp = resp) # Data set that will be used to compute the averages
  
  data.sample = ddply( data, .(usi, date), summarize, 
                       biomass.tot = sum( biomass, na.rm=T), 
                       size.mean = sum( biomass * size, na.rm = T ),
                       ops.mean  = sum( biomass * ops, na.rm = T ),
                       imax.mean = sum( biomass * imax, na.rm = T ),
                       resp.mean = sum( biomass * resp, na.rm = T ) ) # Total biomass in a sample
  data.sample$size.mean = data.sample$size.mean / data.sample$biomass.tot # Mean size in a sample
  data.sample$ops.mean  = data.sample$ops.mean / data.sample$biomass.tot 
  data.sample$imax.mean = data.sample$imax.mean / data.sample$biomass.tot 
  data.sample$resp.mean = data.sample$resp.mean / data.sample$biomass.tot 
  
  # Daily averages, for the whole Wadden Sea
  data.day = ddply( data.sample, .(date), summarize, 
                    biomass = exp( mean( log(biomass.tot), na.rm=T ) ),
                    size = exp( mean( log(size.mean), na.rm = T ) ),
                    ops  = exp( mean( log(ops.mean), na.rm = T ) ),
                    imax = exp( mean( log(imax.mean), na.rm = T ) ),
                    resp = exp( mean( log(resp.mean), na.rm = T ) ),
                    ops.sd  = exp( sd( log(ops.mean), na.rm = T ) ),
                    imax.sd = exp( sd( log(imax.mean), na.rm = T ) ),
                    resp.sd = exp( sd( log(resp.mean), na.rm = T ) ) ) # Mean per day
  
  # return(data.day)
  
  # Monthly averages
  data.day$yr = year(data.day$date) ; data.day$mth = month(data.day$date) # Month and year information

  data.month = ddply(data.day, .(yr, mth), summarize,
                     biomass = mean( log(biomass), na.rm = T ),
                     size = mean( log(size), na.rm = T ),
                     ops  = mean( log(ops),  na.rm = T ),
                     imax = mean( log(imax), na.rm = T ),
                     resp = mean( log(resp), na.rm = T ),
                     nb = sum( mth ) )
  data.month$nb = data.month$nb / data.month$mth # Number of samples per month
  # date.mth = paste( data.month$yr, "-", data.month$mth, "-15", sep = "" )
  # data.month$date = ymd( date.mth )

  # Climatology, as in Barton et al, 2003
  weight.mth = ddply( data.month, .(mth), summarize, tot.nb = sum( nb ) ) # Weights are the number of samples in a month
  data.month = merge( data.month, weight.mth,  by = "mth")
  data.month$wi = data.month$nb / data.month$tot.nb

  # data.sample$biomass = log(data.month$biomass) # Return to log scale for averaging
  # data.sample$size = log(data.month$size)

  climat = ddply(data.month, .(mth), summarize,
                 bio.clim  = sum( biomass * wi, na.rm = T ),
                 size.clim = sum( size * wi, na.rm = T ),
                 ops.clim  = sum( ops  * wi, na.rm = T ),
                 imax.clim = sum( imax * wi, na.rm = T ),
                 resp.clim = sum( resp * wi, na.rm = T ) ) # Climatology

  # Standard deviations
  sd.month = merge(data.month, climat, by = "mth", all.x=T)
  sd.month$sd.bio  = sd.month$wi * ( sd.month$biomass - sd.month$bio.clim )**2
  sd.month$sd.size = sd.month$wi * ( sd.month$size -    sd.month$size.clim )**2
  sd.month$sd.ops  = sd.month$wi * ( sd.month$ops -     sd.month$ops.clim )**2
  sd.month$sd.imax = sd.month$wi * ( sd.month$imax -    sd.month$imax.clim )**2
  sd.month$sd.resp = sd.month$wi * ( sd.month$resp -    sd.month$resp.clim )**2

  sd.clim = ddply(sd.month, .(mth), summarize,
                  sd.bio  = sqrt( sum(sd.bio, na.rm = T) ),
                  sd.size = sqrt( sum(sd.size, na.rm = T) ),
                  sd.ops  = sqrt( sum(sd.ops, na.rm = T) ),
                  sd.imax = sqrt( sum(sd.imax, na.rm = T) ),
                  sd.resp = sqrt( sum(sd.resp, na.rm = T) ) )

  climat = merge(climat, sd.clim, by="mth")

  return(climat)
}# Stop at climatology
  
cop.linear = time.average(cop.b$Event, cop.b$date, cop.b$C..mg.m..3.,
                          cop.b$mean.esd, exp(cop.b$ops.linear), cop.b$imax.linear,
                          cop.b$r.linear)

cop.metab  = time.average(cop.b$Event, cop.b$date, cop.b$C..mg.m..3.,
                          cop.b$mean.esd, exp(cop.b$ops.metab), cop.b$imax.metab,
                          cop.b$r.metab)

plot(ymd( date(cop.linear$date) ), cop.linear$imax, type='l')
lines(ymd( date(cop.metab$date) ), cop.metab$imax, type='l', col='red')

plot(ymd( date(cop.linear$date) ), cop.linear$imax/cop.metab$imax, type='l')
abline(h=1)
abline(h=mean(cop.linear$imax/cop.metab$imax, na.rm=T))

plot(ymd( date(cop.linear$date) ), cop.linear$resp/cop.metab$resp, type='l')
abline(h=1)
abline(h=mean(cop.linear$resp/cop.metab$resp, na.rm=T))

plot(ymd( date(cop.linear$date) ), cop.linear$ops, type='l', log='y')
lines(ymd( date(cop.linear$date) ), cop.metab$ops, col='red')

plot(ymd( date(cop.linear$date) ), cop.linear$size, type='l', log='y')
lesd.shift = agg.par[c('ac_1', 'ac_2')]
abline(h=exp(lesd.shift))

ts.plot = function(ts.linear=cop.linear, ts.metab=cop.metab, lesd.shift=agg.par[c('ac_1', 'ac_2')]){
  
  # Functions for a nice plot
  timeticks = function(date.x, side = 1, tck = 0.02, all=F, labs=T, ...){
    at.ax = range( year( ymd(date.x) ) )
    x.ats = ymd( paste(at.ax[1]:at.ax[2], "-01-15", sep = "") )
    x.labs = paste(at.ax[1]:at.ax[2], "-01", sep = "")
    
    if(labs){
      axis(side = side, at = x.ats, labels = x.labs, tck = tck,  ...)
    }else{
      axis(side = side, at = x.ats, labels = labs, tck = tck,  ...)
    }
    
    if(all){
      for(i in 2:12){
        if(i<10){ii=paste("0", i, sep="")}
        else{ii=i}
        x.labsmid = ymd( paste(at.ax[1]:at.ax[2], "-", ii, "-15", sep = "") )
        axis(side = side, at = x.labsmid, labels = F, tck = tck*0.5, ...)}
    }else{
      x.labsmid = ymd( paste(at.ax[1]:at.ax[2], "-07-01", sep = "") )
      axis(side = side, at = x.labsmid, labels = F, tck = tck*0.5, ...)}
  }
  
  # Error bars
  plot.error.bars = function(mean.x.dt, sd.x.dt, mean.y.dt, sd.y.dt, log=F, col.error=rgb(0.5, 0.5, 0.5, alpha = 0.5)){
    
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
  
  # Line with color gradient 
  line.col.gradient = function(x, y, x.shift, col1, col2, ...){
    ind.shift = which.min( abs(x-x.shift) )
    n = length(x)
    
    col.middle = colorRampPalette( c(col1, col2) )(3)[2] # Middle color between blue and red
    
    col.seq.1 = colorRampPalette( c(col1, col.middle) )(ind.shift)
    col.seq.2 = colorRampPalette( c(col.middle, col2) )(n- ind.shift)
    cols = c(col.seq.1, col.seq.2)
    
    segments(x[-n], y[-n], x[-1], y[-1], col = cols, ...)
  }
  
  # Polygon of relative values 
  pol.relative = function(x, y, limit=1){
    x.seq = seq(min(x), max(x), length.out=100)#min(x):max(x)
    y.seq = approx(x, y, xout=x.seq, na.rm=T)$y
    x.pol = c(x.seq, rev(x.seq))
    
    # Above the limit
    y.high = y.seq
    y.high[y.high<=limit] = limit
    
    # Below the limit
    y.low = y.seq
    y.low[y.low>=limit] = limit
    
    # Low polygon
    polygon(x.pol, c(rep(limit, length(x.seq)), rev(y.low)),
            border=NA, col=adjustcolor('royalblue', alpha=0.5))
    
    # High polygon
    polygon(x.pol, c(rep(limit, length(x.seq)), rev(y.high)),
            border=NA, col=adjustcolor('firebrick2', alpha=0.5))
  }
  
  # Calculate smooth lines
  smooth.lines = function(points, data, span=0.5){
    predni = seq( min(points), max(points), length.out = 1000) # Prediction points
    
    model  = loess( data ~ points, span = span )
    smooth = predict( model, newdata = predni )
    
    return( data.frame(x=predni, y=smooth))
  } 
  
  # Get the colors
  col.line = 'gray40'
  
  # x-axis
  x.lim = range(ts.linear$mth)
  month.labels = data.frame( mth = 1:12,
                             name = c( "Jan", "Feb", "Mar", "April", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec" ),
                             full.name = c( "January", "February", "March", "April", "May", "June", "July", "August",
                                            "September", "October", "November", "December" ))
  
  # Revert to non-log space
  vars.to.exp = c('ops.clim', 'imax.clim', 'resp.clim')
  ts.linear[vars.to.exp] = exp( ts.linear[vars.to.exp] )
  ts.metab[vars.to.exp]  = exp( ts.metab[vars.to.exp] )
  
  par(mfrow=c(4,1), mgp=c(3, 0.1, 0), cex=1.5, mar=c(1, 4, 0.1, 0.1), tck=0.04, pch=19, xpd=F)

  # OPS
  sm.vec = smooth.lines(ts.linear$mth, ts.linear$ops.clim/ts.metab$ops.clim)
  plot(sm.vec$x, sm.vec$y, type='l', col=col.line,
       ann=F, xlim=x.lim, xaxt='n', yaxt='n', lwd=2, xaxs='i', 
       ylim=c( min(sm.vec$y), 2.3) )
  pol.relative(sm.vec$x, sm.vec$y)
  abline(h=1, col=col.line, lwd=2)
  axis(2, at=seq(-2, 2, 0.5), labels=paste( round((seq(-2, 2, 0.5)-1)*100, 0), '%'), las=1 )
  # Mean diff or red area
  mean.diff = mean(ts.linear$ops.clim/ts.metab$ops.clim)-1
  text(x=abs(diff(par('usr')[1:2]))*0.3 + par('usr')[1],
       y=abs(diff(par('usr')[3:4]))*0.2 + par('usr')[3],
       paste('mean difference = ', round( mean.diff*100, 0), '%'), col='white', font=2)
  axis(1, 1:12, labels=F)
  
  mtext(side=3, adj=0.04, 'A - optimal prey size', line=-1.1, cex=1.5, font=2)
  
  # Imax
  sm.vec = smooth.lines(ts.linear$mth, ts.linear$imax.clim/ts.metab$imax.clim)
  plot(sm.vec$x, sm.vec$y, type='l', col=col.line,
       ann=F, xlim=x.lim, xaxt='n', yaxt='n', lwd=2, xaxs='i')
  pol.relative(sm.vec$x, sm.vec$y,)
  abline(h=1, col=col.line, lwd=2)
  axis(2, at=seq(-2, 2, 0.15), labels=paste( round((seq(-2, 2, 0.15)-1)*100, 0), '%'), las=1 )
  # Mean of red area
  rel.diff = ts.linear$imax.clim/ts.metab$imax.clim-1
  mean.diff = mean(rel.diff[rel.diff>0])
  text(x=abs(diff(par('usr')[1:2]))*0.9 + par('usr')[1],
       y=abs(diff(par('usr')[3:4]))*0.6 + par('usr')[3],
       paste(round( (mean.diff)*100, 0), '%'), col='white', font=2)
  # Mean of blue area
  mean.diff = mean(rel.diff[rel.diff<0])
  text(x=abs(diff(par('usr')[1:2]))*0.15 + par('usr')[1],
       y=abs(diff(par('usr')[3:4]))*0.15 + par('usr')[3],
       paste(round( mean.diff*100, 0), '%'), col='white', font=2)
  axis(1, 1:12, labels=F)
  
  mtext(side=2, 'difference of linear vs non-linear scaling', line=2.8, cex=1.5)
  mtext(side=3, adj=0.04, 'B - Imax', line=-1.3, cex=1.5, font=2)

  # Respiration
  sm.vec = smooth.lines(ts.linear$mth, ts.linear$resp.clim/ts.metab$resp.clim)
  plot(sm.vec$x, sm.vec$y, type='l', col=col.line,
       ann=F, xlim=x.lim, xaxt='n', yaxt='n', lwd=2, xaxs='i')
  pol.relative(sm.vec$x, sm.vec$y)
  abline(h=1, col=col.line, lwd=2)
  axis(2, at=seq(1, 2, 0.1), labels=paste( round((seq(1, 2, 0.1)-1)*100, 0), '%'), las=1 )
  # Mean of red area
  rel.diff = ts.linear$resp.clim/ts.metab$resp.clim-1
  mean.diff = mean(rel.diff[rel.diff>0])
  text(x=abs(diff(par('usr')[1:2]))*0.85 + par('usr')[1],
       y=abs(diff(par('usr')[3:4]))*0.6 + par('usr')[3],
       paste(round( (mean.diff)*100, 0), '%'), col='white', font=2)
  # Mean of blue area
  mean.diff = mean(rel.diff[rel.diff<0])
  text(x=abs(diff(par('usr')[1:2]))*0.15 + par('usr')[1],
       y=abs(diff(par('usr')[3:4]))*0.2 + par('usr')[3],
       paste(round( mean.diff*100, 0), '%'), col='white', font=2)
  axis(1, 1:12, labels=F)
  
  mtext(side=3, adj=0.04, 'C - respiration', line=-1.3, cex=1.5, font=2)

  # Mean size ESD
  sm.vec = smooth.lines(ts.linear$mth, ts.linear$size.clim)
  sd.vec = smooth.lines(ts.linear$mth, ts.linear$sd.size)
  plot(sm.vec$x, sm.vec$y, type='n', ann=F, xlim=x.lim, 
       xaxt='n', yaxt='n', xaxs='i', 
       ylim=range(c(sm.vec$y - sd.vec$y, 
                    sm.vec$y + sd.vec$y)) )
  polygon(c(sm.vec$x, rev(sm.vec$x)),
          c(sm.vec$y-sd.vec$y, rev(sm.vec$y+sd.vec$y)),
          border=NA, col=adjustcolor(col.line, alpha=0.25))
  
  lines(sm.vec$x, sm.vec$y, col=col.line, lwd=2)
  # points(ts.linear$mth, ts.linear$size.clim, pch=19) # Check
  # points(ts.linear$mth, ts.linear$size.clim-ts.linear$sd.size, pch=19) # Check
  # points(ts.linear$mth, ts.linear$size.clim+ts.linear$sd.size, pch=19) # Check
  
  axis(2, at=log(c(250, 300, 350, 400)), labels=c(250, 300, 350, 400), las=1 )
  # axis(1, 1:12, labels=month.labels$name)
  axis(1, 3:9, labels=month.labels$name[3:9])
  abline(h=lesd.shift, lty=2, lwd=4, col='gray50')

  text(x=3.3, y=c(lesd.shift[1]*0.99, lesd.shift[1]*1.01),
       labels=c('passive', 'switcher'), col=c(col.ambush, 'darkorange'), adj=0)
  
  mtext(side=3, adj=0.04, 'D - body size', line=-1.1, cex=1.5, font=2)
  mtext(side=2, expression('copepod ESD, ' *mu*'m'), line=2.8, cex=1.5)
}

x11(height=12, width=7)

ts.plot()

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/climatology_difference.pdf')
dev.off()

## Compute the bootstrapping CI for the parameters
bs.parameters = function(modb, pref.max, dat, nDraws=1000){
  
  bs.resample = function(dt){ # Bootstrap resample with replacement
    dt   = na.omit(dt)
    
    # Re-draw the complete dataset, independent on the feeding behavior
    # draw = sample(1:nrow(dt), nrow(dt), replace=T)
    # dt   = dt[draw,]
    
    # Re-draw passive anc active feeders independently
    ind.passive = which( dt$ac=='P' )
    dt.p        = dt[ind.passive,]
    dt.np       = dt[-ind.passive,]

    draw.p      = sample(1:nrow(dt.p),  nrow(dt.p),  replace=T)
    draw.np     = sample(1:nrow(dt.np), nrow(dt.np), replace=T)

    rbind(dt.p[draw.p,], dt.np[draw.np,])
    # Note: this assures that there are passive species in the fit, otherwise the baseline asumptions
    # of the model are broken
  }
  
  names.par           = c('a_f', 'a_shift', 'k_a', 's', 'm', 'a0', 'vdig', 'g', 'r0')
  params.bs           = matrix(NA, ncol=length(names.par), nrow=nDraws)
  colnames(params.bs) = names.par
  
  fit.bs = function(){
  # for(i in 1:nDraws){
  
    modi  = bs.resample( modb[c('lesd', 'lops', 'group', 'ac')] )
    prefi = bs.resample( pref.max[c('lesd', 'Imax.at.15.degreeC..mugC.mugC.1.h.1.', 'group', 'ac')] )
    rdi   = bs.resample( resp.data[c('lesd', 'r', 'ac')] )
    
    all.par = fit.all(modi$lesd,  modi$lops, modi$group,
                      prefi$lesd, prefi$Imax.at.15.degreeC..mugC.mugC.1.h.1., prefi$group,
                      rdi$lesd,   rdi$r,
                      lb = c(1, 1, 1,
                             0, 0, 
                             0, 0,
                             -2, 0),
                      ub = c(5, 10, 100,
                             10, 10,
                             0.3, 1,
                             0, 10),
                      suggestpar = c(4, 5, 20,
                                     0, 5,
                                     0.2, 0.67,
                                     -0.75, 0.1) )
    
    # params.bs[i,] = all.par
    return(all.par)
  }
  
  n.cores = max(1, detectCores() - 4)
  params.bs = mclapply(seq_len(nDraws), function(i) {
    fit.bs()
  }, mc.cores = n.cores)
  params.bs = do.call(rbind, params.bs) # rbind the lists
  
  colnames(params.bs) = names.par
  
  params.mean = apply(params.bs, 2, mean,     na.rm=T)
  params.sd   = apply(params.bs, 2, sd,       na.rm=T)
  params.qt   = apply(params.bs, 2, quantile, probs=c(0.025, 0.5, 0.975))
  
  return( list('mean'=params.mean, 'sd'=params.sd, 'qt'=params.qt) )
}

# @TODO rework so that the resample takes place per size class, not on the FM

bs.sensitivity = bs.parameters(modb, pref.max, dat, nDraws=1000)
for( ni in names(bs.sensitivity) ){
  file.name = paste('~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/bootstrap/bs_params_', ni, 
                    '.csv', sep='')
  write.csv(bs.sensitivity[[ni]], file=file.name, row.names = T)
}

## Conceptual figure of the metabolic speed of copepods
lesd.shift = agg.par[c('ac_1', 'ac_2')]

# Take only the adults
ind.adults.pref = intersect(is.na(pref.max$stage),
                            which(pref.max$stage %in% c('', 'F', 'M', 'CV + F', 'CV + A',
                                                        'A', 'CVI M', 'CVI F')) )

ind.adults.resp = which(resp.data$stage %in% c('?', 'F', 'M', 'A'))

# Only addults
fm.dt = na.omit( rbind( pref.max[ind.adults.pref,   c('species', 'pred.esd', 'ac')],
                        setNames(resp.data[ind.adults.resp, c('species', 'esd', 'ac')], c('species', 'pred.esd', 'ac')) )
)

# All copepods
# fm.dt = na.omit( rbind( pref.max[c('species', 'pred.esd', 'ac')],
#                         setNames(resp.data[c('species', 'esd', 'ac')], c('species', 'pred.esd', 'ac')) )
# )

fm.dt$lesd = log(fm.dt$pred.esd)

## Continuous distribution
# df has columns: size, group
fm.groups = split(fm.dt$lesd, fm.dt$ac)

# common lESD-grid and bandwidth so all groups are smoothed the same way
lesd.seq = seq(min(fm.dt$lesd), max(fm.dt$lesd), length.out = 512)
bw       = bw.nrd0(fm.dt$lesd)

# density for each group, weighted by sample size so group prevalence matters
dens_list = lapply(fm.groups, function(x) {
  d = density(x, bw = bw, from = min(lesd.seq), to = max(lesd.seq), n = length(lesd.seq))
  d$y * length(x)   # weight by n so bigger groups contribute proportionally more
})

dens_fm = do.call(cbind, dens_list)      # Normlized densities of feeding groups
prob_fm = dens_fm / rowSums(dens_fm)     # normalize across groups -> probability at each size
cum_mat  = t(apply(prob_fm, 1, cumsum))      # cumulative sums for stacking

# Plot
x11(height=7, width=14)

x.seq.theoric = seq(lesd.shift[1]-1, lesd.shift[2]+1, length.out=200)

params.a.metab.1 = agg.par[c('ka', 'af_metab', 'ac_1')]
params.a.metab.2 = agg.par[c('ka', 'af_metab', 'ac_2')]
ac = activity(x.seq.theoric, params.a.metab.1) * activity(x.seq.theoric, params.a.metab.2)

# Add limits where feeding mode becomes < 50%
esd.limits = c(lesd.seq[ which(prob_fm[,'P']<=0.5)[1] ], 
               lesd.seq[ which(prob_fm[,'S']<=0.5 & lesd.seq>=log(400))[1] ] )

# Size distribution of foraging modes
par(fig=c(0,0.5, 0,0.95), mgp=c(3, 0.3, 0), mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04, cex=1.6)

# Foraging modes distribution
cols = c('P'=col.ambush, 'S'='darkorange', 'A'=col.filter)
plot(NULL, xlim = range(x.seq.theoric), ylim = c(0, 1),
     ann = F, xaxt='n', yaxt='n',
     xaxs='i', yaxs='i')

# Stripped area, no adult data
polygon(c( c(0,min(lesd.seq)), rev(c(0,min(lesd.seq)))), 
        c( c(0,0), c(10,10)),
        density=4, angle=45, lwd=10,
        col = 'gray90', border = NA)

polygon(c( range(lesd.seq), rev(range(lesd.seq))), 
        c( c(0,0), c(10,10)), col = 'white', border = NA) # Clean background

# Polygons of behaviors
prev = rep(0, length(lesd.seq))
for( i in seq_len(ncol(prob_fm)) ){
  polygon(c(lesd.seq, rev(lesd.seq)), c(prev, rev(cum_mat[, i])),
          col = adjustcolor(cols[ colnames(prob_fm)[i] ], alpha.f=0.7), 
          border = NA)
  prev = cum_mat[, i]
}

abline(v=esd.limits,lty=2, lwd=4, col='white')
abline(h=0.5,lty=1, lwd=4, col='gray40')

axis(2, at=c(0, 0.5, 1), labels=c(0, 0.5, 1), tck=0.04, las=1)
axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=c(100, 200, 500, 1e3, 2e3), tck=0.04)

mtext(side=3, 'A', line=0.3, cex=1.5, font=2, 
      at=par('usr')[1] - abs(diff(par('usr')[1:2]))*0.1, xpd=T)
mtext(side=2, 'frequency of adults observations', line=2., cex=1.6)
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.6, cex=1.6)

# Legend - banner
par(fig=c(0,0.5,0.95,1), mar=c(0, 3.5, 0, 0.1), new=T)

# x.lim = range(lesd.seq)
x.lim = range(x.seq.theoric)
plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, axes=F, xaxs='i') # xaxs='i' for cancelling shift in text

text(x=(esd.limits[1] + min(lesd.seq))/2, y=0.5, 'Passive', col=col.ambush, adj=0.5)
text(x=(esd.limits[1] + esd.limits[2])/2, y=0.5, 'Switcher', col='darkorange', adj=0.5)
text(x=(esd.limits[2] + max(x.lim))/2, y=0.5, 'Active foraging', col=col.filter, adj=0.5)

# Theoretical figure
par(fig=c(0.5,1, 0.,0.95), mgp=c(3, 0.5, 0), mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04, new=T)

plot(x.seq.theoric, ac, lwd=5, 
     col='gray60', ann=F, type='l', yaxt='n', xaxt='n', xaxs='i')

axis(1, at=log(c(100, 200, 500, 1e3, 2e3)), labels=c(100, 200, 500, 1e3, 2e3), tck=0.02)
axis(2, at=c(min(ac), max(ac)),
     labels=c(1, expression(alpha['max'])), las=1 )

abline(v=esd.limits, lty=2, lwd=4, col='gray50')

# polygon(c(esd.limits[1]*c(0.99, 1.01), rev(esd.limits[1]*c(0.99, 1.01))),
#         c(rep(par('usr')[3] + abs(diff(par('usr')[3:4]))*0.7-1.2, 2),
#           rep(par('usr')[3] + abs(diff(par('usr')[3:4]))*0.7+0.2, 2) ),
#         col='white', border=NA) # Remove ablines behind the text
# text(x=rep(esd.limits[1], 3)*c(0.9, 0.92, 0.85), 
#      y=par('usr')[3] + abs(diff(par('usr')[3:4]))*0.7-c(0,1,2),
#      labels=c('+ metabolic rates', '+ swimming', '+ predation risk'), adj=0)
text(x=rep(esd.limits[1], 3)*1.01, 
     y=par('usr')[3] + abs(diff(par('usr')[3:4]))*0.5-c(0,1.2,2.4),
     labels=c('+ metabolic rates', '+ swimming', '+ predation risk'), adj=0, srt=55)

mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.5, cex=1.5)
mtext(side=2, expression('metabolic activity ' *alpha), line=1.2, cex=1.5)
mtext(side=3, 'B', line=0.3, cex=1.5, font=2, 
      at=par('usr')[1] - abs(diff(par('usr')[1:2]))*0.1, xpd=T)

# Legend - banner
par(fig=c(0.5,1,0.95,1), mar=c(0, 3.5, 0, 0.1), new=T)

x.lim = range(x.seq.theoric)
plot(0:1, 0:1, type='n', ann=F, ylim=c(0, 1), xlim=x.lim, axes=F) # xaxs='i' for cancelling shift in text

lines(x = c(0, esd.limits[1])*c(1,0.99), y = rep(0, 2), lty = 1, lwd=4, col=col.ambush)
lines(x = c(esd.limits[1], esd.limits[2]), y = rep(0, 2), lty = 1, lwd=4, col='darkorange')
lines(x = c(esd.limits[2], 10)*c(1.01,1), y = rep(0, 2), lty = 1, lwd=4, col=col.filter)

text(x=(esd.limits[1] + x.lim[1])/2, y=0.5, 'Passive', col=col.ambush, adj=0.5)
text(x=(esd.limits[1] + esd.limits[2])/2, y=0.5, 'Switcher', col='darkorange', adj=0.5)
text(x=(esd.limits[2] + x.lim[2])/2, y=0.5, 'Active', col=col.filter, adj=0.5)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/activity_ESD.pdf')
dev.off()

## Sensitivity of the metabolic activity to the k_a (slope), see bootstrap results
x11(height=4, width=7)
par(mgp=c(3, 0.5, 0), cex=1.5, mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04)

x.seq.theoric = seq(lesd.shift-3, lesd.shift+3, length.out=200)
params.activity.sens = params.activity

# k_a = 10
params.activity.sens[3] = 10
ac = activity(x.seq.theoric, params.activity.sens)
plot(x.seq.theoric, ac, lwd=5, 
     col='gray60', ann=F, type='l', yaxt='n', xaxt='n', xaxs='i')

# k_a = 20
params.activity.sens[3] = 20
ac = activity(x.seq.theoric, params.activity.sens)
lines(x.seq.theoric, ac, lwd=4, col='gray30', lty=2)

# k_a = 50
params.activity.sens[3] = 50
ac = activity(x.seq.theoric, params.activity.sens)
lines(x.seq.theoric, ac, lwd=3, col='black', lty=3)

log10ticks(1, tck=0.02, log=T)
axis( 2, at=c(min(ac), max(ac)),
      labels=c(expression(alpha['i']), expression(alpha['max'])), las=1 )
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.5, cex=1.5)
mtext(side=2, 'metabolic activity', line=1.2, cex=1.5)

legend('topleft', 
       legend=c(expression('k'[alpha]*'=10'),
                expression('k'[alpha]*'=20'),
                expression('k'[alpha]*'=50')),
       bty='n', pch=NA, col=c('gray60','gray30','black'), lwd=5, lty=1:3,
       inset=c(0.1,0.1))

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/activity_ka_sensitivity.pdf')
dev.off()

## Figure showing the error in Imax after the Q10 correction
hist(log(pref.max$Imax.sd / pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.), breaks=20,
     main = expression('Histogram of relative deviation to I'['max']), 
                       las=1, xaxt='n', xlim=c(-10, 2), xlab = '')
log10ticks(1, tck=-0.04, log=T, line=-0.45)
mtext( expression( sigma['measurement'] *' / I'['max']), side=1, line=1.8)
#mtext('log relative error of Imax measurements', side=1, line=1.5)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/Imax_relative_error_measurements.pdf')
dev.off()

## Figure showing the experiment temperature, before correction to 15 C
hist( as.numeric(pref.max$temperature),
     main = expression('Histogram of experiments temperature'), 
     las=1, xaxt='n', xlab = '', col=rgb(0.9, 0.4, 0.4))
axis(1, at=seq(0, 30, 2), tck=-0.02, line=-0.45, labels=F)
axis(1, at=seq(0, 30, 4), tck=-0.04, line=-0.45)

mtext('temperature, °C', side=1, line=1.8)
abline(v=15, lty=3, lwd=2)
mtext(side=3, 'reference temperature, 15°C', line=0)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/experiments_temperature.pdf')
dev.off()

# Impute the Imax error according to the distribution of relative errors
pref.max$Imax.sd[which(is.na(pref.max$Imax.sd))] =
  exp( median(log(pref.max$Imax.sd / pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.), na.rm=T) ) *
  pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.[which(is.na(pref.max$Imax.sd))]

# Calculating the final error in Imax after the Q10 correction
sd.Q10 = sqrt( (3.8-1.9)**2/12 ) # If Q10 follows a uniform law

# Error propagation
pref.max$imax.uncorrected = pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1. / 
  ( 2.8 ** ( 0.1*(as.numeric(pref.max$temperature) - 15)) )

df.dq10  = pref.max$imax.uncorrected * (as.numeric(pref.max$temperature) - 15)/10 * 2.8**( 0.1*(as.numeric(pref.max$temperature) - 15) -1)
df.dimax = 2.8**( 0.1*(as.numeric(pref.max$temperature) - 15))
pref.max$Imax.error.correction = sqrt( (df.dq10 * sd.Q10)**2 + (df.dimax * pref.max$Imax.sd)**2 )

# Plot
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

x11(height=5, width=8)
par(mgp=c(3, 0.5, 0), cex=1.5, mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04)

plot(pref.max$lesd, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., 
     type='n', lwd=2.5, ann=F, xlim=x.lim, xaxt='n', yaxt='n', log='y')
plot.error.bars(pref.max$lesd, 0, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., pref.max$Imax.error.correction,
                col.error = 'gray30')
points(pref.max$lesd, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., pch=19, col='gray60')
log10ticks(1, tck=0.04, log=T)
log10ticks(2, tck=0.04, log=F, las=1)

mtext(side=2, expression('Imax, h'^-1), line=2.1, cex=1.4)
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.6, cex=1.5)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/Imax_error_propagation.pdf')
dev.off()

## Boxplots of copepod sizes in Passive and Active
adult_groups = c("F", "M", "NA", "CVI", "A")
ind.adults = unique( c(grep( paste(adult_groups, collapse='|'), pref.max$stage ),
                       which(pref.max$stage=='' | is.na(pref.max$stage))) )
pref.max$stage[ind.adults] = "A"

pref.max$stage.group=NA
indN=grep("N", pref.max$stage) ; indC=grep("C", pref.max$stage)
pref.max$stage.group[ indN ] = "N"
pref.max$stage.group[ indC ] = "C"
pref.max$stage.group[ intersect(indN, indC) ] = "NC"
pref.max$stage.group[ pref.max$stage == 'A' | is.na(pref.max$stage) ] = "A"

x11(height=5, width=8)
par(mgp=c(3, 0.5, 0), cex=2, mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04)

bp = boxplot(pref.max$pred.esd[pref.max$ac == 'P' & pref.max$stage.group %in% c('C','A')],
             pref.max$pred.esd[pref.max$ac == 'S' & pref.max$stage.group %in% c('C','A')],
             pref.max$pred.esd[pref.max$ac == 'A' & pref.max$stage.group %in% c('C','A')],
             pch=1, log='y',
             col=adjustcolor('darkorange', alpha.f=0.7), # c(col.ambush, col.filter), 
             ann=F, cex.axis=1, las=1, yaxt='n',
             names=c('Passive', 'Switcher', 'Active'))
log10ticks(2, log=F, cex.axis=1, las=1, tck=0.05)
mtext(expression('copepod body size, ESD ' *mu*'m'), side=2, line=2.3, cex=2, adj=1)

for (i in seq_along(bp$names)) {
  # 'i' is the x-position of each box
  text(x = i,
       y = bp$stats[5,i]*0.85, 
       labels = paste0("n=", bp$n[i]), col='black', font=2,
       pos = 3,   # above the point
       cex = 0.8)
}

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/FM_body_size.pdf')
dev.off()


## Plots of OPS and metabolic rates following the cyclopoid / calanoid traits

cyclopoida <- c(
  "Corycaeus anglicus",
  "Diacyclops thomasi",
  "Oithona davisae",
  "Oithona nana",
  "Oithona similis",
  "Oithona spinirostris",
  "Oncaea mediterranea"
)

ind.imax = which(pref.max$species %in% cyclopoida)

plot(pref.max$lesd[-ind.imax], log(pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.[-ind.imax]), 
     pch=19, col=col.filter, lwd=2.5,
     ann=F, xaxt='n', yaxt='n')
points(pref.max$lesd[ind.imax], log(pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.[ind.imax]), 
       pch=1, col=col.ambush, lwd=2.5)
log10ticks(1, tck=0.04, log=T)
log10ticks(2, tck=0.04, log=T, las=1)
mtext(side=2, expression('Imax, h'^-1), line=2.4, cex=1.4)
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=2.5, cex=1.5)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/Imax_error_propagation.pdf')
dev.off()

# Not much difference. Only a few cyclopoida are active feeders in the Imax dataset
# > Oncaea mediterranea and Diacyclops thomasi, otherwise no change 
