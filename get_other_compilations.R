#### Downloads the respiration rates and feeding behavior datasets from Brun at al. (2016) and Pata and Hunt (2023) 
library(gdata)
library(readxl)

## Getting the respiration rates
download.file("https://store.pangaea.de/Publications/BrunP-etal_2016/Brun-etal_2016_Copepode_trait.xlsx", destfile="Brun2016.xlsx")

resp.data = read_excel("Brun2016.xlsx", sheet = 'Respiration rates')

## Getting the feeding behavior datasets
fb.brun = read_excel("Brun2016.xlsx", sheet = 'Feeding mode')

# See https://doi.org/10.5281/zenodo.8377969
file.name = "https://raw.github.com/Pelagic-Ecosystems/Zooplankton_trait_database/master/data_input/Trait_dataset_level1/trait_dataset_level1-2023-08-15.csv"
download.file(file.name,
              destfile="Pata2023.csv", method='wget', mode='w')
fb.pata = read.csv(file='Pata2023.csv')