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

fb.pata = read.feeding_modes.pata( fb.pata ) # Reading Pata and Hunt (2025) data
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

# Check if some species are repeated, but with different references
# feeding.behavior = fb.dat |>
#                    group_by(species, type, detail, download_ref) |>
#                    summarise(primary_ref = paste(primary_ref, collapse = "; "), .groups = "drop")
feeding.behavior = fb.dat[!duplicated(fb.dat[c('species', 'type', 'detail', 'download_ref')]), ]

## Add the feeding behavior classification as 'Passive' or 'Active'
ind.passive = which( (feeding.behavior$type == 'passive' | feeding.behavior$detail %in% c('ambusher')) & 
                      !(feeding.behavior$detail %in% c('filter-cruiser', 'cruiser', 'filter', 'mixed')) & !(feeding.behavior$type %in% c('active')) )

ind.active = which( (feeding.behavior$type == 'active' | feeding.behavior$detail %in% c('filter-cruiser', 'cruiser', 'filter')) &
                     !(feeding.behavior$detail %in% c('ambush', 'mixed')) )

feeding.behavior$ac=NA
# feeding.behavior$ac[feeding.behavior$type == 'mixed' | (is.na(feeding.behavior$type) & feeding.behavior$detail == 'mixed')] = 'S'
feeding.behavior$ac[feeding.behavior$type == 'mixed' | feeding.behavior$detail == 'mixed'] = 'S'
feeding.behavior$ac[ind.passive] = 'P'
feeding.behavior$ac[ind.active] = 'A'

## Remove the feeding modes datasets to save memory
rm(fb.pata, fb.brun, fb.dat)

## Function to add the feeding behavior trait to the copepod datasets
add_feeding_trait = function(data_species, fm.dataset){
  
  # Reset the columns of foraging behavior, if any
  data_species = data_species[, !(names(data_species) %in% c('type', 'detail', 'ac')) ]
  
  # Checking the species that are not in the dataset
  spec.meas = unique( data_species$species )
  spec.dat  = unique( fm.dataset$species )
  notreco   = spec.meas[ which( !(spec.meas %in% spec.dat) ) ] # Names not recognized 
  
  # Merging datasets
  strait = merge( data_species, 
                  fm.dataset[ c("species", "type", "detail", "ac") ], 
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

## Get the foraging behavior of Acartia and Centropages species
# switcher.behavior = read.csv('switcher_behavior.csv')
# 
# add_foraging_switchers = function(dat, sw.dt = switcher.behavior){
#   
#   dat.stage = dat$stage
#   
#   # Homogenize the stage variable to compare
#   adult_groups = c("F", "M", "NA", "CVI", "A")
#   ind.adults = unique( c(grep( paste(adult_groups, collapse='|'), dat$stage ),
#                          which(dat$stage=='' | is.na(dat$stage))) )
#   dat.stage[ind.adults] = "A"
#   
#   # Check the namings of columns, as it can change between datasets
#   ref.name  = c('ref', 'primary.reference')[which(c('ref', 'primary.reference') %in% names(dat))]
#   prey.name = c('prey', 'prey.species')[which(c('prey', 'prey.species') %in% names(dat))]
#   
#   # Find the foraging behavior of switchers
#   for(il in 1:nrow(sw.dt)){
#     swi = sw.dt[il,]
#     
#     ind = intersect( which(dat$species == swi$species & dat.stage == swi$stage ),
#                      grep(swi$primary.ref, dat[[ref.name]], fixed=T) )
#     
#     if(length(ind)>0){
#       f.behav = c('P', 'A')[which( c('passive', 'active') == swi$behavior )]
#       if( !is.na(swi$behavior) ){ # Don't override if nothing is found
#         dat$ac[ind] = f.behav }
#     }
#   }
#   return(dat)
# }
