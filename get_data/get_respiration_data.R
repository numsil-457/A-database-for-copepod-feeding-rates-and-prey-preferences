#### Script that downloads and reads the Respiration data from https://doi.org/10.1594/PANGAEA.862968
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Change path to R script directory

# Convert body mass to carbon
resp.data$`Body mass (mg)`[ resp.data$`Body mass type` == "body mass dry" & !is.na(resp.data$`Body mass (mg)`) ] = 
  resp.data$`Body mass (mg)`[ resp.data$`Body mass type` == "body mass dry" & !is.na(resp.data$`Body mass (mg)`) ] * 0.48
# Conversion factor from Kiorboe (2013) https://doi.org/10.4319/lo.2013.58.5.1843

# Convert respiration rates from muL02 to mugC (as in Serra-Pompei, 2020; https://doi.org/10.1016/j.pocean.2020.102473)
resp.data$`Specific respiration at 15 °C` = resp.data$`Specific respiration at 15 °C` * 0.36 * 1000 / 1e6
# respiration in mgC.mgC-1.h-1

resp.data$mugC = resp.data$`Body mass (mg)`*1000  # In micro gC

convert.mugC.to.esd = function(mugC){
  (mugC / (0.12*10**(-6)) * 3/4/pi)**(1/3) * 2 
}
# As in Hansen et al, (1994) https://doi.org/10.4319/lo.1994.39.2.0395

resp.data$esd = convert.mugC.to.esd(resp.data$mugC) 

# Homogenize the dataset names with the others
resp.data = resp.data[ c('Taxon', 'Life form', 'esd', 'Specific respiration at 15 °C', 'Reference key') ]
names(resp.data) =     c("species", "stage",   "esd", "r",                             "reference")

plot(resp.data$esd, resp.data$r, log = "xy", pch=19, col='gray50',
     ylab='Specific respiration at 15 °C, h-1', 
     xlab=expression('copepod body size, ESD ' *mu*'m') )
