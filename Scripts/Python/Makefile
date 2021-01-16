#Makefile for running workflows/scripts to do a particular function

.PHONY: compile collate fill_accession subset find_raw_differences find_differences find_differences_2011 find_differences_2012 copy_fixes_2011 copy_fixes_2012 diff_spreadsheets_2011 diff_spreadsheets_2012

DEBUGP= 
ifneq ($(DEBUG),)
		DEBUGP=python3 -m pdb
endif

DATA_DIR=../../Data/phenotypic\ data
RAW_DIR=RawData
DERIVED_DIR=DerivedData
BERRY_TEMPLATE_IFOLDER ?= $(DATA_DIR)/$(RAW_DIR)/berry_templates
BERRY_TEMPLATE_OFOLDER ?= $(DATA_DIR)/$(DERIVED_DIR)/berry_templates

compile: collate fill_accession

fill_accession: $(DATA_DIR)/$(DERIVED_DIR)/AccessionMap.xlsx
		./fill_accessions.py -i $(DATA_DIR)/$(DERIVED_DIR)/AccessionMap.xlsx --input_ws "Map" -o $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated.xlsx --output_ws "Combined CNJ Upright Data"

collate: collate_datasets.txt $(DATA_DIR)/$(DERIVED_DIR)/AccessionMap.xlsx
		./collate_all_accessions.sh $(DATA_DIR) $(RAW_DIR) $(DERIVED_DIR) Data-combined-template.xlsx Data-combined-collated.xlsx

subset: $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated.xlsx 
		./random_subset.py -o $(DATA_DIR)/$(DERIVED_DIR)/Data-compiled-collated-subset.xlsx --output_ws "Combined CNJ Upgright Data" -i $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated.xlsx --input_ws "Combined CNJ Upright Data" -p 1	

find_differences:
		./compare_raw_collate_to_precurate.sh $(DATA_DIR) $(DERIVED_DIR) Data-combined-collated.xlsx Data-combined.xlsx Data-combined-template.xlsx Data-diffs.xlsx 5 --convert_year

find_differences_2011:
		./compare_raw_collate_to_precurate.sh $(DATA_DIR) $(DERIVED_DIR) Data-combined-collated.xlsx Data-combined-collated_2011.xlsx Data-combined-template.xlsx Data-diffs-2011.xlsx 3

find_differences_2012:
		./compare_raw_collate_to_precurate.sh $(DATA_DIR) $(DERIVED_DIR) Data-combined-collated.xlsx Data-combined-collated-2012.xlsx Data-combined-template.xlsx Data-diffs-2012.xlsx 3

copy_fixes_2011:
		./copy_curates_to_original_files.py -c $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated_2011.xlsx --curate_ws "Combined CNJ Upright Data" --data_prefix=$(DATA_DIR)/$(RAW_DIR) 2> errors

copy_fixes_2012:
		./copy_curates_to_original_files.py -c $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated-2012.xlsx --curate_ws "Combined CNJ Upright Data" --data_prefix=$(DATA_DIR)/$(RAW_DIR) 2> errors

diff_spreadsheets_2011:
		./diff_spreadsheets.py -1 $(DATA_DIR)/$(RAW_DIR)/2011\ Data/BogLower5-R3.xlsx --file_1_ws "Sheet1" -2 $(DATA_DIR)/$(RAW_DIR)/2011\ Data/BogLower5-R3-curated.xlsx --file_2_ws "Sheet1" -o $(DATA_DIR)/$(RAW_DIR)/2011\ Data/BogLower5-R3-diffs.xlsx --output_ws "Sheet1" 2> errors
		./diff_spreadsheets.py -1 $(DATA_DIR)/$(RAW_DIR)/2011\ Data/BogLower5-R4.xlsx --file_1_ws "Sheet1" -2 $(DATA_DIR)/$(RAW_DIR)/2011\ Data/BogLower5-R4-curated.xlsx --file_2_ws "Sheet1" -o $(DATA_DIR)/$(RAW_DIR)/2011\ Data/BogLower5-R4-diffs.xlsx --output_ws "Sheet1" 2> errors

diff_spreadsheets_2012:
		./diff_spreadsheets.py -1 $(DATA_DIR)/$(RAW_DIR)/2012\ Data/BogLower5-R3-1.xlsx --file_1_ws "Sheet1" -2 $(DATA_DIR)/$(RAW_DIR)/2012\ Data/BogLower5-R3-1-curated.xlsx --file_2_ws "Sheet1" -o $(DATA_DIR)/$(RAW_DIR)/2012\ Data/BogLower5-R3-1-diffs.xlsx --output_ws "Sheet1" 2> errors
		./diff_spreadsheets.py -1 $(DATA_DIR)/$(RAW_DIR)/2012\ Data/BogLower5-R4.xlsx --file_1_ws "Sheet1" -2 $(DATA_DIR)/$(RAW_DIR)/2012\ Data/BogLower5-R4-curated.xlsx --file_2_ws "Sheet1" -o $(DATA_DIR)/$(RAW_DIR)/2012\ Data/BogLower5-R4-diffs.xlsx --output_ws "Sheet1" 2> errors

find_raw_differences:
		./diff_all_spreadsheets.sh ../../Data/phenotypic\ data RawData orig

berry_normalize:
		mkdir -p $(BERRY_TEMPLATE_OFOLDER)
		$(DEBUGP) ./berry_templates_normalize.py -i $(BERRY_TEMPLATE_IFOLDER)/fruit_template_oblong_binary.png -i $(BERRY_TEMPLATE_IFOLDER)/fruit_template_oval_binary.png -i $(BERRY_TEMPLATE_IFOLDER)/fruit_template_pyriform_binary.png -i $(BERRY_TEMPLATE_IFOLDER)/fruit_template_round_binary.png -i $(BERRY_TEMPLATE_IFOLDER)/fruit_template_spindle_binary.png -o $(BERRY_TEMPLATE_OFOLDER)
