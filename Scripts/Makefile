#Makefile for running workflows/scripts to do a particular function

.PHONY: compile collate fill_accession subset 

DATA_DIR=../Data/phenotypic\ data/NJ\ Data

compile: collate fill_accession

fill_accession: $(DATA_DIR)/Compiled/AccessionMap.xlsx
		./fill_accessions.py -i $(DATA_DIR)/Compiled/AccessionMap.xlsx --input_ws "Map" -o $(DATA_DIR)/Compiled/Data-combined-collated.xlsx --output_ws "Combined CNJ Upright Data"

collate: collate_datasets.txt $(DATA_DIR)/Compiled/AccessionMap.xlsx
		./collate_all_accessions.sh $(DATA_DIR) Data-combined-template.xlsx Data-combined-collated.xlsx

subset: $(DATA_DIR)/Compiled/Data-combined-collated.xlsx 
		./random_subset.py -o $(DATA_DIR)/Compiled/Data-compiled-collated-subset.xlsx --output_ws "Combined CNJ Upgright Data" -i $(DATA_DIR)/Compiled/Data-combined-collated.xlsx --input_ws "Combined CNJ Upright Data" -p 1	
