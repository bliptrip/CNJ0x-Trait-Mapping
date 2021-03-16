# CNJ0x Trait QTL Mapping Project

## Author: *Andrew Maule*
## Objective
Map quantitative trait loci (QTL) of cranberry fruiting upright traits in Rutger's elite breeding populations CNJ02 and CNJ04.
Uprights were sampled from bogs at the New Jersey Agricultural Experiment Station (NJAES) of Rutgers University
located in Chatsworth, NJ, for the years 2011-2014 (non-continuous).

## Project Layout
### Data
The `Data/` folder is where the original and cleaned phenotypic and genotypic data are stored.

* Data
    * `genetic_data/` - Genotype marker data and linkage maps.
        * `DerivedData/` - This is an empty placeholder at the moment, as the raw genetic map data provided has not been modified.
        * `RawData/`
            * CNJ02_AllASMapData.csv - The genetic map for the CNJ02 population.
            * CNJ04_AllASMapData.csv - The genetic map for the CNJ04 population.
            * Gryg_AllASMapData.csv  - The genetic map for the GRYG population.
            * consensusMapAll2.csv   - A consensus genetic map of the CNJ02, CNJ04, and GRYG populations.
            * consensusMapAll2.csv.rds
            * consensusMapAll2_withGenes.xlsx - A consensus genetic map with gene annotations.
            * ParentalMaps_HKs.xlsx - Parental genetic maps.
            * SegregationDistortion.xlsx - Segregation distortion of markers.
    * `phenotypic data/` - The phenotypic data.
        * `DerivedData/`
            * `Data-combined-collated.xlsx` - A tidy-version of the data collected for all years and both populations.  This is generated from the 
            `Scripts/Python/collate_all_accessions.sh` script, invoked with `Makefile` target `collate`.
        * `RawData/` - The most important part of each of these folders is the recorded upright data, one per sheet.  These are the raw, unedited
        data collected.
            * `2011 Data/`
            * `2012 Data/`
            * `2013 Data/`
            * `2014 data/`
    * `imgs/` - Contains pictures of NJAES bog 5 lower and upper, where the samples were taken from.
    * `historical_data/` - Chill model and GDD data for NJAES for years 2011-2014.  Potentially useful when modeling year-by-year covariances.

### Scripts
The `Scripts/` folder contains the scripts to preprocess the phenotypic data (predominately in `Python/` subfolder) and fit models, find
significant QTLs, and generate plots and summary tables from these generated datasets (predominately in `R/` subfolder).  Each subfolder has a `Makefile`, with defined targets and user-specified variables that allow user to sequentially and reproducibly execute scripts.

### Workflows
The `Workflows/` folder contains specific configurations and datasets that were run when executing QTL mapping against this population.  Each subfolder represents a different configuration (sets of traits) 
and the associated data derived from calculating QTLs, plots, etc.

* Workflows
    * templates - This folder should be copied to a new unique workflow folder that represents a new configuration
        * configs This is the input configuration folder that scripts use to dictate types of analysis, traits to run permutations and calculate qtls on, etc.  
           * model-traits.cfg.csv - Semicolon-separated file. Mixed model configuration file.  Each row represents a mixed model analysis to do.  Includes options for multitrait analysis, along with configuration parameters to the mixed model function (method to estimate model parameters, etc.)
           * model.cfg - A configuration file used as input when calculating the mixed model BLUPS and parameters when running QTL analysis.  It is meant to be dually interpreted from either python or R.
           * perms.mk - A makefile-formatted file with variables that define how many permutations to run per trait, the number of parallel clusters (CHTC) to use when calculating permutations, the starting seed, and the version of R to use.
        * plots - Any plot template files
           * circos - Folder of circos-related template files
              * drawcircos.html - An interactive html template for displaying a circos graph.  NOTE: This file currently depends on having a localhost http server.
              * drawcircos.js   - The supporting javascript file for drawcircos.html.  This file depends on 
    * \<workflow folder\> - contains all of the same files as templates/ folder, plus the following
        


## Summary of Traits Collected
These upright traits include:

* number of pedicels
* number of pedicels w/o fruit
* number of berries
* number of aborted flowers
* total berry weight
* upright length
* secondary growth
* dry weight of leaves
* largest berry length
* largest berry width
* largest berry weight
* calyx diameter
* calyx lobe form
* calyx lobe size
* calyx end shape
* stem end shape
* berry shape
* number of seeds of largest fruit
* rebud

Because a lot of these traits were gathered over multiple years, using different undergraduate student helpers and sometimes ambiguous classifiers, only some of these traits are truly trustworthy.  Of these the most reliable include total berry length, largest berry length, largest berry width, and largest berry weight, and number of seeds.

