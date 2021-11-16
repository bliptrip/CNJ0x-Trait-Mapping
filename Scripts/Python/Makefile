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

berry_chimera_downsample:
		$(DEBUGP) ./berry_chimera_downsample.py -d 4 -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.png -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.png
		$(DEBUGP) ./berry_chimera_downsample.py -d 4 -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.png -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.png
		$(DEBUGP) ./berry_chimera_downsample.py -d 4 -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.png -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.png
		$(DEBUGP) ./berry_chimera_downsample.py -d 4 -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.png -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.png
		$(DEBUGP) ./berry_chimera_downsample.py -d 4 -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.png -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.png

berry_chimera_templates:
		$(DEBUGP) ./berry_chimera.py -p $(BERRY_TEMPLATE_OFOLDER) -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.npy "oblong"
		$(DEBUGP) ./berry_chimera.py -p $(BERRY_TEMPLATE_OFOLDER) -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.npy "round"
		$(DEBUGP) ./berry_chimera.py -p $(BERRY_TEMPLATE_OFOLDER) -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.npy "oval"
		$(DEBUGP) ./berry_chimera.py -p $(BERRY_TEMPLATE_OFOLDER) -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.npy "spindle"
		$(DEBUGP) ./berry_chimera.py -p $(BERRY_TEMPLATE_OFOLDER) -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.npy "pyriform"

berry_chimeras:
		$(DEBUGP) ./berry_chimera.py -d $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated.xlsx -o $(DATA_DIR)/$(DERIVED_DIR)/Data-combined-collated.all.csv

berry_umcc:
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.npy			-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_5050.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_5050.small.umcc 	  
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_7525.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_7525.small.umcc 	  
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_2575.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_2575.small.umcc 	  
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_5050.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_5050.small.umcc  
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_7525.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_7525.small.umcc  
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_2575.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_2575.small.umcc  
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_5050.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_5050.small.umcc 
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_7525.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_7525.small.umcc 
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_2575.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_2575.small.umcc 
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_5050.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_5050.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_7525.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_7525.small.umcc
		$(DEBUGP) ./berry_umcc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_2575.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_2575.small.umcc

berry_scc:
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.npy			-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_5050.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_5050.small.scc 	  
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_7525.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_7525.small.scc 	  
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_2575.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_2575.small.scc 	  
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_5050.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_5050.small.scc  
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_7525.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_7525.small.scc  
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_2575.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_2575.small.scc  
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_5050.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_5050.small.scc 
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_7525.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_7525.small.scc 
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_2575.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_2575.small.scc 
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_5050.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_5050.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_7525.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_7525.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_2575.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_2575.small.scc
		$(DEBUGP) ./berry_scc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_wavy_binary.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_wavy_binary.small.scc

berry_escc:
		$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oblong_binary.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_round_binary.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.npy			-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_oval_binary.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_spindle_binary.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.npy		-o $(BERRY_TEMPLATE_OFOLDER)/fruit_template_pyriform_binary.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_5050.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_5050.small.escc 	  
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_7525.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_7525.small.escc 	  
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_2575.small.npy 	   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oval_2575.small.escc 	  
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_5050.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_5050.small.escc  
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_7525.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_7525.small.escc  
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_2575.small.npy   -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_oblong_2575.small.escc  
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_5050.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_5050.small.escc 
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_7525.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_7525.small.escc 
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_2575.small.npy  -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_spindle_2575.small.escc 
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_5050.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_5050.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_7525.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_7525.small.escc
		#$(DEBUGP) ./berry_escc.py -i $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_2575.small.npy -o $(BERRY_TEMPLATE_OFOLDER)/fruit_chimera_round_pyriform_2575.small.escc


berry_umcc_plot:
		#$(DEBUGP) ./plot_umcc_diffs.py -i $(BERRY_TEMPLATE_OFOLDER) --names oblong round oval spindle pyriform --files fruit_template_oblong_binary.small.umcc.npy fruit_template_round_binary.small.umcc.npy fruit_template_oval_binary.small.umcc.npy fruit_template_spindle_binary.small.umcc.npy fruit_template_pyriform_binary.small.umcc.npy
		$(DEBUGP) ./plot_umcc_diffs.py -i $(BERRY_TEMPLATE_OFOLDER) --names oblong round oval rov5050 rov7525 rov2575 rob5050 rob7525 rob2575 --files fruit_template_oblong_binary.small.umcc.npy fruit_template_round_binary.small.umcc.npy fruit_template_oval_binary.small.umcc.npy fruit_chimera_round_oval_5050.small.umcc.npy fruit_chimera_round_oval_7525.small.umcc.npy fruit_chimera_round_oval_2575.small.umcc.npy fruit_chimera_round_oblong_5050.small.umcc.npy fruit_chimera_round_oblong_7525.small.umcc.npy fruit_chimera_round_oblong_2575.small.umcc.npy
