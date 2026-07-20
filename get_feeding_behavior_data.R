#### Script that reads the feeding behavior information from the datasets of Brun et al (2016) and Pata and Hunt (2023)

## Read the dataset by Pata and Hunt (2023) and give a table of activity / feeding mode per species
read.feeding_modes.pata = function(dat){
  dat = dat[ which( dat$class == 'Copepoda' & dat$traitName == 'feedingMode' ), ]
  
  trait.db = data.frame()
  for( si in unique(dat$scientificName) ){
    dati = dat[ which(dat$scientificName == si),]
    
    type=NA
    detail=NA
    fm.val = unique(dati$traitValue)
    
    ## Fixed feeding mode
    if( length(grep("active", fm.val))>0 ){type='active'}
    if( length(grep("passive", fm.val))>0 ){type='passive'}
    if( length(grep("active", fm.val))>0 & length(grep("passive", fm.val))>0 ){ type='mixed' }
    
    if( length(grep("ambush", fm.val))>0 ){detail='ambusher'}
    if( length(grep("cruise", fm.val))>0 ){detail='cruiser'}
    if( length(grep("current", fm.val))>0 ){detail='filter'}
    if( length(grep("current", fm.val))>0 & length(grep("cruise", fm.val))>0 ){detail='filter-cruiser'}
    if( length(grep("ambush", fm.val))>0 & length(grep("cruise", fm.val))>0 ){detail='ambusher-cruiser'}
    if( length(grep("current", fm.val))>0 & length(grep("ambush", fm.val))>0 ){detail = 'mixed'}
    
    if( length(grep("particle feeder", fm.val))>0 ){detail='particle feeder'}
    if( length(grep("parasite", fm.val))>0 ){detail='parasite'}
    
    dlref = paste(c( unique(dati$secondaryReference), 'Pata2023'), collapse='; ')
    
    df.si = data.frame(species = si, type=type, detail=detail, 
                       download_ref = dlref, 
                       primary_ref = dati$primaryReference)
    trait.db = rbind(trait.db, df.si)
  }
  
  # Clean up the parenthesis in the species names
  trait.db$species = gsub("\\s*\\([^\\)]+\\)", "", trait.db$species)
  
  return(trait.db)
}

read.feeding_modes.brun = function(trait){
  names(trait)[1] = "species"
  
  trait$type = NA
  trait$detail = NA
  
  for(si in 1:length(trait$species)){
    feedings = trait[si, c("Passive", "Mixed", "Active") ]
    detail = trait[si, c("Feeding current", "Cruise feeder", "Ambush feeder") ]
    detail[ is.na(detail) ] = 0
    
    type = find_activity_fm_brun(feedings, detail)
    detail = type[2] ; type=type[1]
    
    if( is.character(type) ){trait$type[si] = type}
    if( is.character(detail) ){trait$detail[si] = detail}
  }
  
  # Clean up the parenthesis in the species names
  trait$species = gsub("\\s*\\([^\\)]+\\)", "", trait$species)
  
  trait$download_ref = 'Brun2016'
  trait$primary_ref = trait$`Reference key`
  
  return( trait[c('species', 'type', 'detail', 'download_ref', 'primary_ref')] )
}

# Merge the feeding mode from the database of Brun et al (2016)
find_activity_fm_brun = function(feedings, detail){
  type=NA
  fm = NA
  
  # Activity
  if( feedings[2] == 1 ){
    fm = "mixed"
    type = "mixed"
  }else if( feedings[1] == 1 ){ type="passive"
  }else if( feedings[3] == 1 ){ type="active" }
  
  # Feeding mode
  if(feedings[2] != 1){
    if(detail[1]==1 & detail[2]!=1 & detail[3]!=1){ fm = "filter"
    }else if(detail[1]!=1 & detail[2]==1 & detail[3]!=1){ fm = "cruiser"
    }else if(detail[1]!=1 & detail[2]!=1 & detail[3]==1){ fm = "ambusher"
    }else if( detail[1]==1 & detail[2]==1 & detail[3]!=1 ){fm = "filter-cruiser"}
  }
  
  return( c(type, fm) )
}

fb.pata = read.feeding_modes.pata( fb.pata ) # Reading Pata and Hunt (2023) data
fb.brun = read.feeding_modes.brun( fb.brun ) # Reading Brun et al (2016) data

# Remove the repeated species names in Pata2023
d = strsplit(fb.pata$species, split=" ")
for(i in 1:length(d)){
  d[[i]] = paste(unique(d[[i]]), collapse = ' ')
}
fb.pata$species=unlist(d)

## Check which species of Brun (2016) are not in the Pata and Hunt (2023) dataset
notreco = unique(fb.brun$species)[ which( !(unique(fb.brun$species) %in% unique(fb.pata$species)) ) ] 

# Further check if the name did not change formatting
not.at.all = c()
for( si in notreco ){
  ind=grep( si, unique(fb.pata$species) )
  if(length(ind)==0){
    not.at.all = c( not.at.all, si )
  }
}

# Add them the missing species the dataset
ind = which(fb.brun$species %in% not.at.all)

fb.dat = rbind(fb.pata, fb.brun[ind,])
fb.dat = fb.dat[ - which( duplicated(fb.dat) ), ] # Some species are duplicated

# Check if some species are repeated, but just with different references
feeding.behavior = fb.dat |>
                   group_by(species, type, detail, download_ref) |>
                   summarise(primary_ref = paste(primary_ref, collapse = "; "), .groups = "drop")

## Add the feeding behavior classification as 'Passive' or 'Active'
ind.passive = which( (feeding.behavior$type == 'passive' | feeding.behavior$detail %in% c('ambusher')) & 
                      !(feeding.behavior$detail %in% c('mixed', 'filter')) & !(feeding.behavior$type %in% c('active')) )

ind.active = which( (feeding.behavior$type == 'active' | feeding.behavior$detail %in% c('filter-cruiser', 'cruiser', 'filter')) &
                     !(feeding.behavior$detail %in% c('mixed', 'ambush')) )

feeding.behavior$ac=NA
feeding.behavior$ac[feeding.behavior$type == 'mixed' | feeding.behavior$detail == 'mixed'] = 'S'
feeding.behavior$ac[ind.passive] = 'P'
feeding.behavior$ac[ind.active] = 'A'

## Remove the feeding modes datasets to save memory
rm(fb.pata, fb.brun, fb.dat)