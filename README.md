#CNJ0x Trait QTL Mapping Project

##Author: *Andrew Maule*
##Objective
Map traits taken from uprights sent to the Zalapa lab in the years 2011-2014 from Vorsa's CNJ02, CNJ04, and their parental populations.

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

Because a lot of these traits were gathered over multiple years, using different undergraduate student helpers and sometimes ambiguous classifiers, only some of these traits are truly trustworth.
Of these the most trustworthy include total berry length, largest berry length, largest berry width, and largest berry weight, and number of seeds.

##Project Organization
###Data
The data folder is where the original and cleaned phenotypic and genotypic data are stored.

* Data
    * genetic_data - The genotypic data.
        * DerivedData - This is an empty placeholder at the moment, as the raw genetic map data provided has not been modified.
        * RawData
            * CNJ02_AllASMapData.csv - The genetic map for the CNJ02 population.
            * CNJ04_AllASMapData.csv - The genetic map for the CNJ04 population.
            * Gryg_AllASMapData.csv  - The genetic map for the GRYG population.
            * consensusMapAll2.csv   - A consensus genetic map of the CNJ02, CNJ04, and GRYG populations.
            * consensusMapAll2.csv.rds
            * consensusMapAll2_withGenes.xlsx - A consensus genetic map with gene annotations.
            * ParentalMaps_HKs.xlsx - Parental genetic maps.
            * SegregationDistortion.xlsx - Segregation distortion of markers.
    * phenotypic_data - The phenotypic data.
        * DerivedData - The 
        * RawData - The most important part of each of these folders is the recorded upright data, one per sheet
            * 2011 Data
            * 2012 Data
            * 2013 Data
            * 2014 data

###Workflows
The Workflows/ folder contains specific configurations and datasets that were run when executing QTL mapping against this population.  Each subfolder represents a different configuration (sets of traits) 
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
        

