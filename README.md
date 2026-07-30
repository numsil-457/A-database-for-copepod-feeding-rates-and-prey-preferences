# Overall description

This repository contains the code necessary to reproduce the figures and analysis of: Imbert T., Garcia-Oliva O., and Wirtz K. (Year). Title. Journal.

The grazing rates were digitized from previous publications in the literature. The copepod and prey sizes were extracted from the respective studies. See the file 'metadata.csv' for more details on the compilation of grazing experiments. 

The present analysis also uses open access data published in the literature.

# Data availability

- The compilation of grazing rates is detailed in 'metadata.csv'.
- The feeding behavior information of copepods was extracted from the studies by Brun et al. (2016), available on PANGAEA (https://doi.pangaea.de/10.1594/PANGAEA.862968), and by Pata and Hunt (2025; https://aslopubs.onlinelibrary.wiley.com/doi/10.1002/lno.12478).
- The specific respiration rates of copepods were also extracted from the study of Brun et al. (2016)

# Script files

The required software is R. This repository contains X scripts.

- 'main.R' creates the figures and fits the size-based model to the prey preferences, and grazing and respiration rates of copepods.

- 'get_other_compilations.R' downloads the online databases of Brun et al. (2016) and Pata and Hunt (2025).

- 'get_feeding_behavior_data.R' extracts the feeding behavior information from the downloaded datasets.

- 'get_respiration_data.R' extracts the specific respiration rates from the downloaded datasets.
