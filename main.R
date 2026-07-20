#### Script to adapt the specialization model and the 2013 Imax formula to the copepods dataset Imax/OPS

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Change path to R script directory

library(plyr)
library(dplyr)
library(nloptr)
library(mixtools)
library(parallel)
library(gdata)
library(readxl)

## Load other source files
source("get_other_compilations.R")    # Download datasets from Brun et al. (2016) and Pata and Hunt (2023)
source("get_respiration_data.R")      # Extract the respiration rates from Brun et al. (2016)
source("get_feeding_behavior_data.R") # Extract the feeding behavior information

source("~/PhD/Library_scripts/Rfunctions_plot_TS.R") # Functions for plotting

## Function to add the feeding behavior trait to the metabolic rates datasets
add_feeding_trait = function(data_species, fm.dataset){

  # Checking the species that are not in the dataset
  spec.meas = unique( data_species$species )
  spec.dat  = unique( fm.dataset$species )
  notreco   = spec.meas[ which( !(spec.meas %in% spec.dat) ) ] # Names not recognized 
  
  # Merging datasets
  strait = merge( data_species, 
                  trait[ c("species", "type", "detail", "ac") ], 
                  by = "species", all.x = T )
  #strait = strait[ - which( duplicated(strait) ), ] # Some rows are duplicated
  
  # Filling for the non recognized phylum (due to formatting)
  for( ni in notreco ){
    indni = which( strait$species == ni )
    type=NA
    detail=NA
    ac=NA
    
    if( length( grep("/", ni) ) > 0 ){ # Two genera given
      n0 = strsplit(ni, "/")[[1]] # In case it is a species name
      
      n1 = n0[1]
      n1 = strsplit(n1, " ")[[1]][1] # In case it is a species name
      #n1 = n1[ -which(n1=="") ]   # Do not take white spaces into account
      
      n2 = n0[2]
      n2 = strsplit(n2, " ")[[1]]
      indn2 = -which(n2=="")
      if( length( indn2 ) > 0 ){
        n2 = n2[ indn2 ][1] # Do not take white spaces into account
      }else{
        n2 = n2[1]
      }
      ind = c(which(n1==fm.dataset$species), which(n2==fm.dataset$species)) #c( grep( n1, trait$species ), grep( n2, trait$species ) )
      
    }else{ # Only one genera is given
      ni = strsplit(ni, " ")[[1]][1] # In case it is a species name
      ind = which(ni==fm.dataset$species) #grep( ni, trait$species )
    }
    
    if( length(ind) > 0 ){
      type   = unique(fm.dataset$type[ind])   # If one is detected, then it is marked
      detail = unique(fm.dataset$detail[ind]) # If one is detected, then it is marked
      ac     = unique(fm.dataset$ac[ind])
      
      if( is.character(type) ){strait$type[ indni ] = type}
      if( is.character(detail) ){strait$detail[ indni ] = detail}
      if( is.character(ac) ){strait$ac[ indni ] = ac}
      
    }
  }
  
  return(strait)
}

# Checking - Fine to use the feeding.behavior from the script! ----
trait = read.csv('~/PhD/Work/Copepods project/feeding_modes.csv')

modb = read.csv(file = "~/PhD/Work/Copepods project/Latex-feeding-modes/copepod_OPS_database_feeding_behavior.csv", header = T)

modb = modb[c('species')]

modb.t = add_feeding_trait(modb, trait)
modb.f = add_feeding_trait(modb, feeding.behavior)

modb.t$ac == modb.f$ac
which( !(modb.t$ac == modb.f$ac) )
ind = which(is.na(modb.t$ac == modb.f$ac))
modb.t[ind,]
modb.f[ind,] 

# >> No big change

pref = read.csv(file = "~/PhD/Work/Copepods project/NEW_TAKE_SCRIPTS/copepod_fmax_imax.csv")

# Remove any complicated datasets
datasets.not.to = c("Uye and Kasahara (1983)", "Storms (1974)", "Rao and Kumar (2002)", "Vogt et al (2013)")
# Probably not Imax
pref = pref[ -which(pref$primary.reference %in% datasets.not.to), ]
pref = pref[ -which(pref$prey.size.unit != "ESD"), ] # Remove Sommer and Sommer (2006)

## Select the max Ingestion rate per study and species/stage; a bit less evolved than the OPS detection
pref$ind.row = row(pref)[,1]
pref.max = ddply(pref[which(!is.na(pref$Imax.at.15.degreeC..mugC.mugC.1.h.1.)),], 
                 .(species, stage, primary.reference), summarize,
                 i.max = ind.row[ which.max(Imax.at.15.degreeC..mugC.mugC.1.h.1.) ] )

pref.max = pref[pref.max$i.max,]

pref.mt = add_feeding_trait(pref.max, trait)
pref.mf = add_feeding_trait(pref.max, feeding.behavior)

pref.mt$ac == pref.mf$ac
which( !(pref.mt$ac == pref.mf$ac) )
ind = which(is.na(pref.mt$ac == pref.mf$ac))
pref.mt[ind, c('species', 'ac', 'primary.reference')]
pref.mf[ind, c('species', 'ac', 'primary.reference')] 

# >> No big change ----

## Colors per feeding mode - @TODO change names
col.filter = rgb(0.2, 0.3, 0.9, alpha = 0.9)
col.ambush = rgb(0.9, 0.1, 0.2, alpha = 0.9)
col.switcher = 'gray50'

## Read the OPS dataset
modb = read.csv(file = "~/PhD/Work/Copepods project/Latex-feeding-modes/copepod_OPS_database_feeding_behavior.csv", header = T)
modb$lops = log(modb$ops)
modb$lesd = log(modb$pred)

modb$s = modb$lops - mean(modb$lops, na.rm=T)
plot(modb$lesd, modb$s)

## Split the OPS values into groups
# Fit a gaussian mixture model with the AIC criteria to find the specialization groups, using the s
wmult = makemultdata(modb$lops, cuts = seq(1, 7, by=0.5)) # Prepared subsets of data to fit mixed gaussians on
multmixmodel.sel(wmult, comps=1:length(seq(1, 7, by=0.5)), epsilon=1e-3) # AIC selects 2 kernels
w1 = normalmixEM(modb$lops, lambda = 0.5, mu = c(3,6), sigma=1)

# Choose the group by using the posterior probabilities of observation
g=rep(NA, nrow(modb))
g[w1$posterior[,1] < w1$posterior[,2]] = 2
g[w1$posterior[,1] >= w1$posterior[,2]] = 1

modb$gp = g
plot(modb$lesd, modb$lops, col = modb$gp, pch=19)

# Give a 'Low' or 'High' criteria for the groups
gp.mops = ddply(modb, .(gp), summarize, mops = mean(lops, na.rm=T))
gp.mops$group = c('low', 'high')[c( which.min(gp.mops$mops), which.max(gp.mops$mops) )]

modb = merge(modb, gp.mops[c('gp', 'group')], by='gp')

## Read the Imax data

# Plot the raw Imax dataset
pref = read.csv(file = "~/PhD/Work/Copepods project/NEW_TAKE_SCRIPTS/copepod_fmax_imax.csv")

# Remove any complicated datasets
datasets.not.to = c("Uye and Kasahara (1983)", "Storms (1974)", "Rao and Kumar (2002)", "Vogt et al (2013)")
# Probably not Imax
pref = pref[ -which(pref$primary.reference %in% datasets.not.to), ]
pref = pref[ -which(pref$prey.size.unit != "ESD"), ] # Remove Sommer and Sommer (2006)

# Removing the Imax calculated with the linear formula, uncertainty is much too large
# ind.no.imax = which(pref$Imax.method=='Linear')
# pref$Imax.at.15.degreeC..mugC.mugC.1.h.1.[ind.no.imax] = NA

pref = add_feeding_trait(pref)

## Select the max Ingestion rate per study and species/stage; a bit less evolved than the OPS detection
pref$ind.row = row(pref)[,1]
pref.max = ddply(pref[which(!is.na(pref$Imax.at.15.degreeC..mugC.mugC.1.h.1.)),], 
                 .(species, stage, primary.reference), summarize,
                 i.max = ind.row[ which.max(Imax.at.15.degreeC..mugC.mugC.1.h.1.) ] )

pref.max = pref[pref.max$i.max,]

## Give an OPS group to Imax

# Based on distance of prey_esd to the mean OPS
pref.max$group = NA
dist.high = abs( log(pref.max$prey.size) - gp.mops$mops[gp.mops$group=='high'] )
dist.low  = abs( log(pref.max$prey.size) - gp.mops$mops[gp.mops$group=='low'] )

pref.max$group[dist.high < dist.low]  = 'high'
pref.max$group[dist.high >= dist.low] = 'low'

pref.max$pch = 19
pref.max$pch[pref.max$group=='high']=1

## Adding the feeding behavior traits to the respiration dataset
resp.data = add_feeding_trait(resp.data)

## Function to add the plot coding (pch, color) to the dataset
add_graphics_trait = function(dt, is.group.ops=F){
  
  if(is.group.ops){
    dt$colo                   = col.filter
    dt$colo[dt$group=='high'] = col.ambush
  }
  
  dt$pch = 18 # Other
  dt$pch[dt$ac=='P'] = 1
  dt$pch[dt$ac=='A'] = 19

  return(dt)
}

modb      = add_graphics_trait(modb,     is.group.ops=T)
pref.max  = add_graphics_trait(pref.max, is.group.ops=T)
resp.data = add_graphics_trait(resp.data)

## Functions to fit a new specialization model
activity = function(lesd, params){
  a_f=params[1]; a_shift=params[2]; k_a=params[3]
  
  #a = ( d_multiplicator + exp(-(lesd**2-a_shift)) ) / ( 1 + exp(-(lesd**2-a_shift)) )
  a = ( a_f + exp(-k_a*(lesd-a_shift)) ) / ( 1 + exp(-k_a*(lesd-a_shift)) )
  
  return( a )
}

ops.specialization = function(lesd, params, params.activity, linear=F){
  specialization=params[1]; m=params[2]
  #lops.estim = exp(-specialization**2) * lesd + m1 # Original
  #lops.estim = exp(-specialization**2) * lesd / (1 + exp(-specialization**2) * lesd) *m1 +m2
  #lops.estim = exp(-specialization**2)*lesd**2 / (m2 + exp(-specialization**2)*lesd**2) *m1 #+ m2
  #lops.estim = exp(-specialization**2) * lesd + m2/(1 + exp(-f*lesd)) + m1
  #lops.estim = exp(-specialization**2) * lesd + m2/(1 + exp(-(lesd**2-f))) + m1

  if(!linear){ # Bifurcation
    lops.estim = exp(-specialization**2) * lesd + m * activity(lesd, params.activity)
  }else{
    lops.estim = exp(-specialization**2) * lesd + m 
  }
  
  return(lops.estim)
}

imax.new = function(params, params.activity, lesd, lops, n){
  a0 = params[1] # See Wirtz, 2013
  v = params[2]
  
  #lops = a1 * lops
  a = n * a0 * (lesd -(-1)**n * lops)
  
  a.lim=2
  a[which(a < 0.0)] = 0
  a[which(a > a.lim)] = a.lim
  
  vdig = v * activity(lesd, params.activity)
  
  # Eq.4 in Wirtz JPR 2013 with additional factor e^8: 172h-1*24*e^-8 = 1.4 d-1
  Imax = 6 * vdig * exp( (a + (a.lim - a) * lops + (a - a.lim-1) * lesd) )
  
  return( list(Imax, a) )
}

resp.imax = function(params, params.activity, lesd){
  gamma = params[1]
  r0    = params[2]
  
  resp.mod = r0 * activity(lesd, params.activity) * exp(gamma*lesd)
  return(resp.mod)
}

## Fit the OPS, Imax, and Resp with one framework
fit.all = function(lesd.ops, y.ops, ops.group, lesd.imax, y.imax, imax.group, lesd.resp, y.resp, lb, ub, suggestpar, n=1){
  
  RMSE = function(params, lesd.ops, y.ops, lesd.imax, y.imax, lesd.resp, y.resp, n0){
    params.activity = params[1:3]
    params.ops  = params[4:5]
    params.imax = params[6:7]
    params.resp = params[8:9]
    
    #ops.mod  = ops.specialization(lesd.ops, params.ops) # No bifurcation
    ind.ops.low  = which(lesd.ops < params.activity[2] | ops.group=='low') # All points before the shift of in the 'low' group
    ind.ops.high = which(lesd.ops < params.activity[2] | ops.group=='high') 
    
    ops.mod.high = ops.specialization(lesd.ops, params.ops, params.activity)
    ops.mod.low  = ops.specialization(lesd.ops, params.ops, params.activity, linear=T)
    
    # imax.mod = imax.new(params.imax, lesd.imax, 
    #                     ops.specialization(lesd.imax, params.ops), n0)[[1]]
    ind.imax.low  = which(lesd.imax < params.activity[2] | imax.group=='low') # All points before the shift of in the 'low' group
    ind.imax.high = which(lesd.imax < params.activity[2] | imax.group=='high') 
    
    imax.mod.high = imax.new(params.imax, params.activity, lesd.imax, 
                             ops.specialization(lesd.imax, params.ops, params.activity), n0)[[1]]
    imax.mod.low = imax.new(params.imax, params.activity, lesd.imax, 
                            ops.specialization(lesd.imax, params.ops, params.activity, linear = T), n0)[[1]]
    
    resp.mod = resp.imax(params.resp, params.activity, lesd.resp)
    
    rmse.calculate = function(x, y, normalize=T){
      #if(sd.use) sd.y = sd(y, na.rm=T) else sd.y = 1
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
    rmse.ops.high = rmse.calculate(y.ops[ind.ops.high], ops.mod.high[ind.ops.high])
    rmse.ops.low  = rmse.calculate(y.ops[ind.ops.low],  ops.mod.low[ind.ops.low])
    rmse.ops = rmse.ops.high + rmse.ops.low 
    
    #rmse.imax = rmse.calculate(y.imax, imax.mod)
    rmse.imax.high = rmse.calculate(y.imax[ind.imax.high], imax.mod.high[ind.imax.high])
    rmse.imax.low  = rmse.calculate(y.imax[ind.imax.low],  imax.mod.low[ind.imax.low])
    rmse.imax = rmse.imax.high + rmse.imax.low 
    
    rmse.resp = rmse.calculate(y.resp, resp.mod)
    rmse = ( rmse.ops + rmse.imax + rmse.resp ) /3
    
    return( rmse )
  } 
  
  mod = isres(x0 = suggestpar, fn = RMSE,
              lower = lb, upper = ub, maxeval = 5e5L,
              lesd.ops=lesd.ops, y.ops=y.ops,
              lesd.imax=lesd.imax, y.imax=y.imax, 
              lesd.resp=lesd.resp, y.resp=y.resp, n0=n)
  
  suggestpar = mod$par
  lb = suggestpar * 0.9 * ( (1-sign(mod$par))*0.1 + 1 ) 
  ub = suggestpar * 1.1 * ( -(1-sign(mod$par))*0.1 + 1 )
  
  mod = lbfgs(x0 = suggestpar, fn = RMSE, 
              lower = lb, upper = ub, 
              lesd.ops=lesd.ops, y.ops=y.ops,
              lesd.imax=lesd.imax, y.imax=y.imax, 
              lesd.resp=lesd.resp, y.resp=y.resp, n0=n)
  
  return( mod$par )
}

pref.max$lesd  = log(pref.max$pred.esd)
resp.data$lesd = log(resp.data$esd)

all.par = fit.all(modb$lesd, modb$lops, modb$group,
                  pref.max$lesd, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., pref.max$group,
                  resp.data$lesd, resp.data$r,
                  lb = c(1, 1, 10,
                         0, 0, 
                         0, 0,
                         -2, 0),
                  ub = c(5, 10, 100,
                         5, 10,
                         0.3, 1,
                         0, 10),
                  suggestpar = c(4, 5, 20,
                                 0, 5,
                                 0.2, 0.67,
                                 -0.75, 0.1) )

## Plot
params.activity = all.par[1:3]
x.seq           = seq(3, 9, length.out=150)
x.lim           = c(4, 8.5)
lesd.shift      = params.activity[2]

#ops.seq  = ops.specialization(x.seq, all.par[1:5]) # No bifurcation
ops.seq.high = ops.specialization(x.seq, all.par[4:5], params.activity)
ops.seq.low  = ops.specialization(x.seq, all.par[4:5], params.activity, linear=T)

imax.seq.low = imax.new(all.par[6:7], params.activity, x.seq, ops.seq.low, 1)[[1]]
imax.seq.high = imax.new(all.par[6:7], params.activity, x.seq, ops.seq.high, 1)[[1]]

resp.seq = resp.imax(all.par[8:9], params.activity, x.seq)

# Gradient color lines
line.col.gradient = function(x, y, x.shift, ...){
  ind.shift = which.min( abs(x-x.shift) )
  n = length(x)
  
  col.middle = colorRampPalette( c(col.filter, col.ambush) )(3)[2] # Middle color between blue and red
  
  col.seq.1 = colorRampPalette( c(col.filter, col.middle) )(ind.shift)
  col.seq.2 = colorRampPalette( c(col.middle, col.ambush) )(n- ind.shift)
  cols = c(col.seq.1, col.seq.2)
  
  segments(x[-n], y[-n], x[-1], y[-1], col = cols, ...)
}

x11(height=12, width=7)
par(mfrow=c(3,1), mgp=c(3, 0.5, 0), cex=1.5, mar=c(2.5, 4, 0.1, 0.1), tck=0.04)

# OPS
plot(modb$lesd, modb$lops, col=adjustcolor(modb$colo, alpha.f=0.7), pch=modb$pch, lwd=2.5,
     ann=F, ylim=c(1, 7), xlim=x.lim, xaxt='n', yaxt='n')
log10ticks(1, tck=0.04, log=T, axis.lab=F)
log10ticks(2, tck=0.04, log=T, las=1)

abline(v = lesd.shift, col='gray50', lty=3, lwd=3)
lines(x.seq, ops.seq.high, col='white', lwd=6)
line.col.gradient(x.seq, ops.seq.high, lesd.shift, lwd=4)
lines(x.seq, ops.seq.low, col='white', lwd=6)
lines(x.seq, ops.seq.low, col=col.filter, lwd=4)
mtext(side=2, expression('OPS, ESD' *mu *'m'), line=2.5, cex=1.4)

lines(x.seq, x.seq + log(0.1), col='gray20', lwd=4, lty=2)

# Legend
legend('topleft', c('passive', 'active', 'other/NA'),
       pch=c(1, 19, 18), col='black', lwd=3, lty=NA,
       bty='n', x.intersp=-0.4, cex=1, inset=c(0, -0.05))

text(x=log(1700), y=log(800), 'High OPS', col=col.ambush, srt=5)
text(x=log(2000), y=log(45),  'Low OPS',  col=col.filter, srt=5)

mtext(side=3, adj=0.99, 'A', line=-1, cex=1.5, font=2)

# Imax
plot(pref.max$lesd, log(pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.), 
     pch=pref.max$pch, col=adjustcolor(pref.max$colo, alpha.f=0.7), lwd=2.5,
     ann=F, xlim=x.lim, xaxt='n', yaxt='n')
log10ticks(1, tck=0.04, log=T, axis.lab=F)
log10ticks(2, tck=0.04, log=T, las=1)

abline(v = lesd.shift, col='gray50', lty=3, lwd=3)
lines(x.seq, log(imax.seq.high), col='white', lwd=6)
line.col.gradient(x.seq, log(imax.seq.high), lesd.shift, lwd=4)
lines(x.seq, log(imax.seq.low), col='white', lwd=6)
lines(x.seq, log(imax.seq.low), col=col.filter, lwd=4)
mtext(side=2, expression('Imax, h'^-1), line=2.5, cex=1.4)

#lm.imax = lm(log(Imax.at.15.degreeC..mugC.mugC.1.h.1.) ~ lesd, data=pref.max)
#abline(lm.imax, lty=2, col='gray10', lwd=3)
saiz.model = log( exp(0.225)* ( (exp(x.seq)**3*0.523)/(8.3*1e6) )**(0.703-1)/24 )
lines(x.seq, saiz.model, col='gray20', lwd=4, lty=2)

# Legend
legend('bottomright', 'linear model', 
       pch=NA, col='gray20', lwd=4, lty=2,
       bty='n', x.intersp=0.5, inset=c(0.15, 0.1))

mtext(side=3, adj=0.99, 'B', line=-1, cex=1.5, font=2)

# Respiration
plot(resp.data$lesd, log(resp.data$r), pch=resp.data$pch, col='gray70', lwd=2.5,
     ann=F, xlim=x.lim, xaxt='n', yaxt='n')
log10ticks(1, tck=0.04, log=T)
log10ticks(2, tck=0.04, log=T, las=1)

abline(v = lesd.shift, col='gray50', lty=3, lwd=3)
lines(x.seq, log(resp.seq), lwd=4, col='black')
mtext(side=2, expression('respiration rate, h'^-1), line=2.5, cex=1.4)
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.5, cex=1.4)

#lm.resp = lm(log(value) ~ lesd, data=dat)
#abline(lm.resp, lty=2, col='gray20', lwd=4)

lm.resp = lm(log(r) + 0.75 * lesd ~ 1, data=resp.data) # Specified slope
abline(lm.resp$coefficients, -0.75, lty=2, col='gray20', lwd=4)

mtext(side=3, adj=0.99, 'C', line=-1, cex=1.5, font=2)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/physio_model_copepod.pdf')
dev.off()

## ANOVA for the slope of respiration rates (size limit at 800 mu m ESD)
dat$size_group = 's'
dat$size_group[dat$esd >= 600] = 'l'

linear_model = lm(log(value) ~ log(esd) * size_group, data = dat[dat$ac=='A',])
anova(linear_model)

## MAEs of Ingestion models
mae = function(y.obs, y.model){
  sum( abs(y.obs-y.model), na.rm=T ) / length( na.omit(y.obs) )
  #sqrt( sum( (y.obs-y.model)**2, na.rm=T ) / length( na.omit(y.obs) ) )
}

# Saiz and Calbet (2007) model
saiz.pred = log( exp(0.225)* ( (pref.max$pred.esd**3*0.523)/(8.3*1e6) )**(0.703-1)/24 )

# Non-linear model
ops.seq.high = ops.specialization(pref.max$lesd[pref.max$group=='high'], all.par[4:5], params.activity)
ops.seq.low  = ops.specialization(pref.max$lesd[pref.max$group=='low'], all.par[4:5], params.activity, linear=T)

imax.seq.low  = imax.new(all.par[6:7], params.activity, pref.max$lesd[pref.max$group=='low'], ops.seq.low, 1)[[1]]
imax.seq.high = imax.new(all.par[6:7], params.activity, pref.max$lesd[pref.max$group=='high'], ops.seq.high, 1)[[1]]

nlm.pred = log( c(imax.seq.low, imax.seq.high) )
nlm.y    = log( c(pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.[pref.max$group=='low'],
                  pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.[pref.max$group=='high']) )

# MAEs
mae.saiz = mae( log(pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.), saiz.pred )
mae.nlm  = mae( nlm.y, nlm.pred )

## Compute the bootstrapping CI for the parameters
bs.parameters = function(modb, pref.max, dat, nDraws=1000){
  
  bs.resample = function(dt){ # Bootstrap resample with replacement
    dt   = na.omit(dt)
    draw = sample(1:nrow(dt), nrow(dt), replace=T)
    dt   = dt[draw,]
  }
  
  names.par           = c('a_f', 'a_shift', 'k_a', 's', 'm', 'a0', 'vdig', 'g', 'r0')
  params.bs           = matrix(NA, ncol=length(names.par), nrow=nDraws)
  colnames(params.bs) = names.par
  
  fit.bs = function(){
  # for(i in 1:nDraws){
  
    modi  = bs.resample( modb[c('lesd', 'lops', 'group')] )
    prefi = bs.resample( pref.max[c('lesd', 'Imax.at.15.degreeC..mugC.mugC.1.h.1.', 'group')] )
    rdi   = bs.resample( resp.data[c('lesd', 'r')] )
    
    all.par = fit.all(modi$lesd,  modi$lops, modi$group,
                      prefi$lesd, prefi$Imax.at.15.degreeC..mugC.mugC.1.h.1., prefi$group,
                      rdi$lesd,   rdi$r,
                      lb = c(1, 1, 1,
                             0, 0, 
                             0, 0,
                             -2, 0),
                      ub = c(5, 10, 100,
                             5, 10,
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

bs.sensitivity = bs.parameters(modb, pref.max, dat, nDraws=1000)
for( ni in names(bs.sensitivity) ){
  file.name = paste('~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/bootstrap/bs_params_', ni, 
                    '.csv', sep='')
  write.csv(bs.sensitivity[[ni]], file=file.name, row.names = T)
}

## Conceptual figure of the metabolic speed of copepods
x11(height=5, width=9)
par(mgp=c(3, 0.5, 0), cex=1.5, mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04)

x.seq.theoric = seq(lesd.shift-3, lesd.shift+3, length.out=200)
ac = activity(x.seq.theoric, params.activity)
plot(x.seq.theoric, ac, lwd=5, 
     col='gray60', ann=F, type='l', yaxt='n', xaxt='n', xaxs='i')

log10ticks(1, tck=0.02, log=T)
axis( 2, at=c(min(ac), max(ac)),
      labels=c(expression(alpha['i']), expression(alpha['max'])), las=1 )
mtext(side=1, expression('copepod body size, ESD ' *mu *'m'), line=1.5, cex=1.5)
mtext(side=2, 'metabolic activity', line=1.2, cex=1.5)

text(x = c(min(x.seq.theoric), max(x.seq.theoric)) * c(1.4, 0.9), 
     y = c(max(ac), min(ac)) * c(0.97, 1.05), 
     c('Passive', 'Active'), cex = 1.2)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/activity_ESD.pdf')
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
x11(height=5, width=8)
par(mgp=c(3, 0.5, 0), cex=2, mar=c(2.5, 3.5, 0.1, 0.1), tck=0.04)

bp = boxplot(pref.max$pred.esd[pref.max$ac == 'P'], pref.max$pred.esd[pref.max$ac == 'A'],
             pch=1, log='y',
             col=adjustcolor('darkorange', alpha.f=0.7), # c(col.ambush, col.filter), 
             ann=F, cex.axis=1, las=1, yaxt='n',
             names=c('Passive', 'Active'))
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


### Beyond is subject to deletion

## Showing the change in metabolic speed of copepods
dat$value.detrended = dat$value / exp(-0.75*log(dat$esd))
pref.max$imax.detrended = pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1. / exp(-0.75*log(pref.max$pred.esd))

plot(dat$lesd, dat$value.detrended, pch=19, col=dat$fm.col, ann=F, log='y')
plot(pref.max$lesd, pref.max$imax.detrended, pch=19, col=pref.max$fm.col, ann=F, log='y')


x11(height=12, width=7)
par(mfrow=c(2,1), mgp=c(3, 0.5, 0), cex=1.5, las=1, mar=c(2.5, 5, 0.5, 0.5), tck=0.04)

plot(dat$lesd, dat$value.detrended/max(dat$value.detrended, na.rm=T), 
     pch=19, col=dat$fm.col, ann=F, log='y')
title(ylab='Detrended respiration rate, scaled', xlab='', line=3, cex=1.4)

plot(pref.max$lesd, pref.max$imax.detrended/max(pref.max$imax.detrended, na.rm=T), 
     pch=19, col=pref.max$fm.col, ann=F, log='y')
title(ylab='Detrended ingestion rate, scaled', xlab='', line=3, cex=1.4)
mtext(side=1, 'log ESD', line=1.5, cex=1.4)

dev.copy2pdf(file='~/PhD/Work/Copepods project/Latex-feeding-modes/Imax_model_adaptation/detrended_metabolism.pdf')
dev.off()

## Bootstrap of the OPS parameters
bs.sample = function(data, nDraws=1000, nParams=3, namesParams=NULL, funcBS){
  xseq = seq(2, 10, by=0.02)

  params.list <- matrix(ncol=nParams, nrow=nDraws)
  result.list <- matrix(ncol=length(xseq), nrow=nDraws)
  
  # Bootstrap
  for(ii in 1:nDraws){
    draw = sample(1:nrow(data), nrow(data), replace=T)
    
    params.bootstrap = funcBS(xseq, data[draw,])
    
    # Collect the results
    result.list[ii,] = params.bootstrap[[2]]
    params.list[ii,] = params.bootstrap[[1]]
  }
  
  # Naming and formatting of the bootstrap datasets
  params.list = as.data.frame(params.list)
  result.list = as.data.frame(result.list)

  names(result.list) = xseq
  if (is.null(namesParams)) namesParams=1:nParams
  names(params.list) = namesParams

  # Calculate the interquartile range and the sd, as well as the parameters sensitivity
  model.uncertainty = data.frame(lesd = xseq)
  names.f = 'lesd'

  res.mean = apply(result.list, 2, 'mean', na.rm=T)
  res.sd   = apply(result.list, 2, 'sd', na.rm=T)
  res.qt   = apply(result.list, 2, 'quantile', probs=c(0.05, 0.5, 0.95), na.rm=T)
  
  model.uncertainty$m     = res.mean
  model.uncertainty$sd    = res.sd
  model.uncertainty$qtmin = res.qt[1,]
  model.uncertainty$qtmed = res.qt[2,]
  model.uncertainty$qtmax = res.qt[3,]
    
  # Parameters
  par.mean = apply(params.list, 2, 'mean', na.rm=T)
  par.sd   = apply(params.list, 2, 'sd', na.rm=T)
  #par.qt   = apply(params.list, 2, 'quantile', probs=c(0.05, 0.5, 0.95), na.rm=T)
  
  params.uncertainty = cbind(par.mean, par.sd)
  names(params.uncertainty) = c( paste(namesParams, '.m', sep=''),
                                 paste(namesParams, '.sd', sep='') )
                                 #paste(namesParams, '.m', sep='') )
  
  return( list(model.uncertainty, params.uncertainty) ) 
}

bs.ops.sample = function(xseq, data){
  params.bootstrap = fit.specialization(data$lesd, data$lops)
  ops.bootstrap    = ops.specialization(xseq, params.bootstrap)

  return( list(params.bootstrap, ops.bootstrap) )
}

ind = which(!is.na(modb$imax))

bs.ops = bs.sample(modb, nDraws=100, nParams=4, namesParams=c('s', 'm1', 'm2', 'f'), bs.ops.sample)

plot(modb$lesd, modb$lops, type='n', lwd=2, xlab='', ylab='', las=1, ylim=c(1, 8))

polygon(x=c( bs.ops[[1]][,1], rev(bs.ops[[1]][,1]) ),
        y=c( bs.ops[[1]][,2]-bs.ops[[1]][,3], rev(bs.ops[[1]][,2]+bs.ops[[1]][,3]) ),
        col='gray80', border=NA)
# polygon(x=c( bs.ops[[1]][,1], rev(bs.ops[[1]][,1]) ),
#         y=c( bs.ops[[1]][,4], rev(bs.ops[[1]][,6]) ),
#         col='gray80', border=NA)

points(modb$lesd, modb$lops, pch=19, col='gray50')
lines(bs.ops[[1]][,1], bs.ops[[1]][,2], lwd=2)

## Using exactly the same data
pref.max$lops = log(pref.max$prey.size)
bs.ops = bs.sample(pref.max, nDraws=100, nParams=4, namesParams=c('s', 'm1', 'm2', 'f'), bs.ops.sample)

plot(pref.max$lesd, pref.max$lops, type='n', lwd=2, xlab='', ylab='', las=1, ylim=c(1, 8))

polygon(x=c( bs.ops[[1]][,1], rev(bs.ops[[1]][,1]) ),
        y=c( bs.ops[[1]][,2]-bs.ops[[1]][,3], rev(bs.ops[[1]][,2]+bs.ops[[1]][,3]) ),
        col='gray80', border=NA)

points(pref.max$lesd, pref.max$lops, pch=19, col='gray50')
lines(bs.ops[[1]][,1], bs.ops[[1]][,2], lwd=2)

## Change the fitting process slightly, so that the OPS parameters are also adapted
fit.ops.imax = function(lesd.imax, y.imax, lesd.ops, y.ops, lb, ub, suggestpar, n=1){
  
  RMSE = function(params, lesd.imax, y.imax, lesd.ops, y.ops, n0){
    params.ops  = params[1:4]
    params.imax = params[5:7]
    
    ops   = ops.specialization(lesd.ops, params.ops)
    y.mod = imax.new(params.imax, lesd.imax, ops, n0)[[1]]
    
    #rmse = sum( abs( log(y.mod) - log(y.obs) ), na.rm=T ) / length(lesd) # MAE
    rmse = ( sqrt( sum( ( log(y.mod) - log(y.imax) )**2, na.rm=T ) / sd(log(y.mod), na.rm=T) ) / length(lesd.imax) +
             sqrt( sum( (ops - y.ops)**2, na.rm=T ) / sd(ops, na.rm=T) ) / length(lesd.ops) ) / 2 # RMSD
    # rmse = ( sqrt( sum( ( log(y.mod) - log(y.imax) )**2, na.rm=T ) ) / length(lesd.imax) +
    #          sqrt( sum( (ops - y.ops)**2, na.rm=T ) ) / length(lesd.ops) ) / 2 # RMSD
    
    # rmse = ( sum( abs( log(y.mod) - log(y.imax) ), na.rm=T ) / sd(log(y.mod), na.rm=T) / length(lesd.imax) +
    #          sum( abs(ops - y.ops), na.rm=T ) / sd(ops, na.rm=T) / length(lesd.ops) ) / 2 # RMSD
    # rmse = ( sum( abs( log(y.mod) - log(y.imax) ), na.rm=T ) / length(lesd.imax) +
    #          sum( abs(ops - y.ops), na.rm=T ) / length(lesd.ops) ) / 2 # RMSD
    
    #rmse = sqrt( sum( ( log(y.mod) - log(y.imax) )**2, na.rm=T ) ) / length(lesd.imax) 
    
    return( rmse )
  } 
  
  mod = isres(x0 = suggestpar, fn = RMSE,
              lower = lb, upper = ub, maxeval = 5e5L,
              lesd.imax=lesd.imax, y.imax=y.imax, 
              lesd.ops=lesd.ops, y.ops=y.ops, n0=n)
  
  suggestpar = mod$par
  lb = suggestpar * 0.9 * ( (1-sign(mod$par))*0.1 + 1 ) 
  ub = suggestpar * 1.1 * ( -(1-sign(mod$par))*0.1 + 1 )
  
  mod = lbfgs(x0 = suggestpar, fn = RMSE, 
              lower = lb, upper = ub, 
              lesd.imax=lesd.imax, y.imax=y.imax, 
              lesd.ops=lesd.ops, y.ops=y.ops, n0=n)
  
  return( mod$par )
}

imax.par = fit.ops.imax(pref.max$lesd,
                        pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1.,
                        modb$lesd,
                        modb$lops,
                        lb = c(0.3, 0.,  1, 32, 0,   0,   0),
                        ub = c(1.7, 0.8, 5, 48, 0.4, 10, 2),
                        suggestpar = c(1, 0.05, 3, 35, 0.1, 0.67, 0.2) )

ops.seq = ops.specialization(x.seq, imax.par[1:4])

plot(modb$lesd, modb$lops, pch=19, col='gray60', ann=F, log='y')
lines(x.seq, ops.seq, lwd=2)

imax.seq = imax.new(imax.par[5:7], x.seq, ops.specialization(x.seq, imax.par[1:4]), 1)

plot(pref.max$lesd, pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1., 
     pch=19, col='gray60', ann=F, log='y')
lines(x.seq, imax.seq[[1]], lwd=2)
title(ylab='Imax', xlab='ESD', line=2, cex=1.4)


## Fixing the Imax with the OPS: maybe the ESD != OPS - does not fix!
ind.adult = unique( c(grep('A', pref.max$stage), which(pref.max$stage %in% c('', 'F', 'M') )) )
ind.adult = unique( c(ind.adult, which( is.na(pref.max$stage) )) )
pref.max$stage[ ind.adult ] = 'A'

pref.max$lops = NA # Loop to affect the OPS
for(i in 1:nrow(pref.max)){
  indi = which(pref.max$species[i] == modb$species & pref.max$stage[i] == modb$stage)
  
  if(length(indi) > 0){
    ind.ops = indi[which.min( (modb$lesd[indi] - log(pref.max$prey.size[i]))**2 )] # Keeps the index traceable
    pref.max$lops[i] = modb$lops[ind.ops]
  }
}

pref.max$imax.fix = pref.max$Imax.at.15.degreeC..mugC.mugC.1.h.1. * exp(3/2 * (pref.max$lops - log(pref.max$prey.size))**2)
plot(pref.max$pred.esd, pref.max$imax.fix, log='xy', pch=19, ylim=c(1e-3, 1))
