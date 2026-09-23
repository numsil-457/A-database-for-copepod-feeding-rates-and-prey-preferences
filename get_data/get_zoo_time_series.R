#### Downloads the zooplankton time series in the East Frisian Wadden Sea from Imbert et al (2026; https://doi.org/10.3389/fmars.2026.1858560)
library(worrms)
library(lubridate)

## Download datasets
if( !file.exists('WS_2015_2024_zooplankton.tab') ){
  download.file("https://doi.pangaea.de/10.1594/PANGAEA.995417?format=textfile", 
                destfile="WS_2015_2024_zooplankton.tab")
}

# Read the dataset
zoo.ts = read.csv('WS_2015_2024_zooplankton.tab', sep = '\t',header = TRUE, skip=43)
aphia_id = as.numeric(sub(".*:taxname:", "", zoo.ts$Species.UID..Semantic.URI.)) # IDs for searching in WORMS

# Convert date to readable format
zoo.ts$date = ymd( date(zoo.ts$Date.Time) )

# Calculate the mean ESD
zoo.ts$mean.volume = zoo.ts$Biovol.v..mm..3.m..3. / zoo.ts$Abund.v....m..3. # In mm3
zoo.ts$mean.esd    = (3/4/pi * zoo.ts$mean.volume*1e9) ** (1/3) * 2 # In micro m

zoo.ts$class = NA
for( ai in unique(aphia_id) ){
  class = wm_classification(ai)
  class = class[ which(class$rank=='Class'), 'scientificname' ]
  
  if(dim(class)[1]<1){
    class=NA
  }
  
  zoo.ts$class[ which(ai == aphia_id ) ] = class
}

zoo.ts$class = unlist(zoo.ts$class)

# Copepod time series
cop.ts = zoo.ts[which(zoo.ts$class=='Copepoda'),]

# Fill in the carbon weight, when missing (using the factors of Kiorboe, 2013; https://aslopubs.onlinelibrary.wiley.com/doi/10.4319/lo.2013.58.5.1843)
cop.ts$C..mg.m..3.[is.na(cop.ts$C..mg.m..3.)] = 0.48 * cop.ts$Dry.m..mg.m..3.[is.na(cop.ts$C..mg.m..3.)]

# Gelatinous zooplankton time series
gel.ts = zoo.ts[which(zoo.ts$class %in% c('Scyphozoa', 'Hydrozoa', 'Gymnolaemata', 'Tentaculata',
                                          'Nuda')),]

# Split gelatinous zooplankton into subgroups to calculate the carbon content
gel.ts$gel.group = NA
gel.ts$gel.group[which(gel.ts$class %in% c('Scyphozoa', 'Hydrozoa'))] = 'cnidaria'
gel.ts$gel.group[which(gel.ts$class %in% c('Gymnolaemata', 'Tentaculata', 'Nuda'))] = 'ctenophora'

gel.ts$C..mg.m..3.[ which(is.na(gel.ts$C..mg.m..3.) & gel.ts$gel.group == "cnidaria") ] =
  0.132 * gel.ts$Dry.m..mg.m..3.[ which(is.na(gel.ts$C..mg.m..3.) & gel.ts$gel.group == "cnidaria") ]

gel.ts$C..mg.m..3.[ which(is.na(gel.ts$C..mg.m..3.) & gel.ts$gel.group == "ctenophora") ] =
  0.051 * gel.ts$Dry.m..mg.m..3.[ which(is.na(gel.ts$C..mg.m..3.) & gel.ts$gel.group == "ctenophora") ]
