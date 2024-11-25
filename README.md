# Ecological harshness mediates reproductive trade-offs in a wild bird population

This repository hosts data, R code, and Stan code for Bliard L, Martin JS, Childs DZ, Cole EF, Sheldon BC, Paniw M, Ozgul A. Ecological harshness mediates reproductive trade-offs in a great tit population.

Preprint version: 

The repository is archived on Zenodo at 

## GENERAL INFORMATION

1. Title: Data and scripts from "Ecological harshness mediates reproductive trade-offs in a great tit population".

2. Author Information:
	
        A.  Name: Louis Bliard
		Institution: University of Zurich
		Address: Winterthurerstrasse 190, 8057 Zurich, Switzerland
		Email: bliard.louis@gmail.com
	
        B.  Name: Jordan S Martin
		Institution: Eawag - Swiss Federal Institute of Aquatic Science and Technology
		
        C.  Name: Dylan Z Childs
		Institution: University of Sheffield

        D.  Name: Ella F Cole
		Institution: University of Oxford

        E.  Name: Ben C Sheldon
		Institution: University of Oxford

        F.  Name: Maria Paniw
		Institution: Estación Biológica de Doñana
		
        G.  Name: Arpat Ozgul
		Institution: University of Zurich
		Email: arpat.ozgul@uzh.ch
		
		
4. Date of data collection: 1961-2018

5. Geographic location of data collection: 51°77′N, 1°32′W


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC-BY 4.0

2. Links to publications that cite or use the data: NA

3. Links to other publicly accessible locations of the data: SPI-Birds https://nioo.knaw.nl/en/projects/spi-birds/study/wytham-woods 

4. Was data derived from another source? NA

5. Recommended citation for this dataset: Bliard L, Martin JS, Childs DZ, Cole EF, Sheldon BC, Paniw M, Ozgul A (XXXX). Ecological harshness mediates reproductive trade-offs in a great tit population [Data set].



## DATA & FILE OVERVIEW

1. File List: 
- `data_wytham_tits.RData`
- `df.txt`
- `df_fecundity.txt`

- `analysis_wytham_tits.R`
- `model_wytham_tits.stan`

- `posterior_draws.txt`


2. Relationship between files, if important: 

The RData file `data_wytham_tits.RData` contains the datasets `df.txt` and `df_fecundity.txt`, as well as 20 alternative versions of these datasets with imputed missing data `mi.df` and `mi_df_fecundity`. These datasets are needed to run the Stan model `model_wytham_tits.stan` using the R script `analysis_wytham_tits.R`. The non-imputed version of the dataset are also provided as txt file.

The text file `posterior_draws.txt` contains the output of the model. More precisely, it contains the polled posterior draws (60000 posterior samples per parameter) for all model paramters of interest.


## METHODOLOGICAL INFORMATION
 
1. Methods for processing the data: R. Only the formatted data for analysis is provided. Note that individuals were anonymised (they will not match the individual IDs from other great tit publications). Please contact the data custodian in charge of the great tit data (Ben Sheldon) if you wish to access the great tit and blue tit database for your own analyses.

2. Instrument- or software-specific information needed to interpret the data: 
- R v.4.3.2 https://www.r-project.org/
- CmdStanR v.0.8.1 https://mc-stan.org/cmdstanr/

3. People involved with sample collection: GREAT TIT PEOPLE

4. People involved with data formatting and analysis: Bliard L.

5. for more general informations regarding the methods, see https://doi.org/10.32942/X2D89H and the related code https://github.com/Jordan-Scott-Martin/covariance-reaction-norms and https://doi.org/10.1111/1365-2656.14173 and the related code https://github.com/lbiard/detecting_tradeoffs_crn_models 


### DATA-SPECIFIC INFORMATION FOR: `data_wytham_tits.RData`

1. Number of objects: 4

2. Objects List: 
- "df.txt" = data frame containing data for the offspring mass model.
- "df_fecundity.txt" = data frame containing data for the fecundity and recruitment models.
- "mi_df" = list of 20 data frames, with same format as "df.txt" but with imputed missing data.
- "mi_df_fecundity" =  list of 20 data frames, with same format as "df_fecndity.txt" but with imputed missing data.


### DATA-SPECIFIC INFORMATION FOR: `df.txt`

1. Number of variables: 19

2. Number of cases/rows: 53753
   Each row correspond to one offspring.

3. Variable List: 
- "BreedingSeason" = year of the observation.
- "FemaleID" = unique identifier of the mother.
- "BroodIDLaid" = unique identifier of the brood.
- "IndivID" = unique identifier of the offpsring.
- "CapturePlot" = location within Wytham woods.
- "Mass.x" = mass of the offpsring (in grams).
- "ChickAge" = age of the offspring at the time of mass measurement.
- "Plot" = location within Wytham woods.
- "LocationID" = unique identifier of the nextbox.
- "LayDate_observed" = lay date of the brood in Julian days starting on April 1st.
- "BroodSize_observed" = size of the brood.
- "Mass.y" = mass of the mother (in grams).
- "breeding_age" = age of the mother (0 = 1st year bird; 1 = 2+ year bird).
- "spring_temperature" = daily mean temperature from March 1st to May 9th (in Celsius).
- "spring_precipitation" = sum of precipitation from April 1st to May 31st (in mm).
- "population_density" = number of females hatching at least one egg in the given breeding season.
- "Half_fall_date" = day of peak abundance of winter moth larvae.
- "mast.score" = beech mast index, scored as an ordinal variable of increasing beech mast from 0 to 2 
- "synchrony" = difference (in days) between the half fall date and the mother's laying date.

4. Missing data codes: NA


### DATA-SPECIFIC INFORMATION FOR: `df_fecundity.txt`

1. Number of variables: 16

2. Number of cases/rows: 7287
   Each row correspond to one brood.

3. Variable List:
- "BroodID" = unique identifier of the brood.
- "BreedingSeason" = year of the observation.
- "FemaleID" = unique identifier of the mother.
- "Plot" = location within Wytham woods.
- "LocationID" = unique identifier of the nextbox.
- "LayDate_observed" = lay date of the brood in Julian days starting on April 1st.
- "BroodSize_observed" = size of the brood.
- "Mass" = mass of the mother (in grams).
- "breeding_age" = age of the mother (0 = 1st year bird; 1 = 2+ year bird).
- "spring_temperature" = daily mean temperature from March 1st to May 9th (in Celsius).
- "spring_precipitation" = sum of precipitation from April 1st to May 31st (in mm).
- "population_density" = number of females hatching at least one egg in the given breeding season.
- "n_recruits" = number of offpsring from the brood that have recruited into the population in following years.
- "Half_fall_date" = day of peak abundance of winter moth larvae.
- "mast.score" = beech mast index, scored as an ordinal variable of increasing beech mast from 0 to 2 
- "synchrony" = difference (in days) between the half fall date and the mother's laying date.

4. Missing data codes: NA


### DATA-SPECIFIC INFORMATION FOR: `posterior_draws.txt`

1. Number of variables: 62

2. Number of cases/rows: 60000
   Each row correspond to one posterior sample.
