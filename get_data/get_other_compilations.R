#### Downloads the respiration rates and feeding behavior datasets from Brun at al. (2016) and Pata and Hunt (2023) 

## Download datasets
if( !file.exists('Brun2016.xlsx') ){
  download.file("https://store.pangaea.de/Publications/BrunP-etal_2016/Brun-etal_2016_Copepode_trait.xlsx", destfile="Brun2016.xlsx")
}

if( !file.exists('Pata2023.csv') ){ # See https://doi.org/10.5281/zenodo.8377969
  file.name = "https://raw.github.com/Pelagic-Ecosystems/Zooplankton_trait_database/master/data_input/Trait_dataset_level1/trait_dataset_level1-2023-08-15.csv"
  download.file(file.name,
                destfile="Pata2023.csv", method='wget', mode='w')
}

## Getting the respiration rates
resp.data = read_excel("Brun2016.xlsx", sheet = 'Respiration rates')

## Getting the feeding behavior datasets
fb.brun = read_excel("Brun2016.xlsx", sheet = 'Feeding mode')

fb.pata = read.csv(file='Pata2023.csv')