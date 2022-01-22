#Makefile for building/packaging data & scripts for HTCondor submission

.PHONY: install GenCross PermMerge Perm2DMerge Rbuild PermBuild StepwiseQTLBuild StepwiseQTL

HTCONDOR=HTCondorFiles
SANDBOX?=../../sandbox
CIRCOS_FOLDER=circos-cnj02

export WORKFLOW ?= ../../Workflows/1
export P1_Name ?= Mullica_Queen
export P2_Name ?= Crimson_Queen
export POP_Name ?= CNJ02
export POP_Circos = #Force user to define
export RDS_DIR ?= "../../Data/phenotypic\ data/DerivedData/cleanup_data.R.output/Robjs"
export COLLATE_EFFECTS ?= TRUE

#Set the following for configuring the BLUPs summary table in target GenQTLTable: NUM_TOP configures how many of the top-ranked QTLs to display per trait:model combination
QTL_METHOD ?= stepwiseqtl
NUM_TOP ?= 2

install:
		R CMD BATCH packagedeps.R

means:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_00_1_means.R 

plotcors:
		R CMD BATCH qtl_pipeline_00_2_plotcors.R 

genmodels:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_01_genmodels.R

genanalyses:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_01_analyzemodels.R

genblupplots:
		touch "$(WORKFLOW)/traits/blups.summary.table.png"
		touch "$(WORKFLOW)/traits/raw.summary.table.png"
		R CMD BATCH '--args workflow="$(WORKFLOW)" P1_Name="$(P1_Name)" P2_Name="$(P2_Name)" POP_Name="$(POP_Name)"' qtl_pipeline_01_03_genBLUPplots.R

genh2plot:
		R CMD BATCH '--args workflow="$(WORKFLOW)" POP_Name="$(POP_Name)"' qtl_pipeline_01_04_genh2plots.R

genbrcorrs:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_01_05_genBLUPCorrPlots.R

#Rbuild tars up the packages necessary for installing an R-build and then transferring to submit server.  Untar Rpackage.tgz, submit an interactive session using
#condor_submit -i Rpackage.sub
#Comment-out any packages you don't want to install if they're already in the tarball.
Rbuild:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

#Rdownload - Pulls the generaed RPackage file from the submit server.
Rdownload:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@
		
#Permbuild is meant for building a tarball of necessary scripts/etc. needed for running single-QTL (scanone) permutation test on cluster.
PermBuild:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

#RunPerms will extract the permutations folder and submit the condor jobs.
RunPerms:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

#PullPerms will extract the calculated permutations from the submit server
PullPerms:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

PermMerge:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_03_mergperms.R

ScanOne:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_04_scanone.R

#Perm2DBuild is meant for building a tarball of necessary scripts/etc. needed for running 2D-QTL (scantwo) permutation tests on cluster.
Perm2DBuild:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

#RunPerms2D will extract the 2D permutations folder and submit the condor jobs.
RunPerms2D:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

#PullPerms2D will extract the calculated 2D permutations from the submit server
PullPerms2D:
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

Perm2DMerge:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_06_mergeperms2D.R

#StepwiseQTLBuild is meant for bundling/building an HTCondor dataset for running stepwiseqtl()
StepwiseQTLBuild::
		make -C $(HTCONDOR) WORKFLOW=../$(WORKFLOW) $@

#StepwiseQTL is meant to run stepwiseqtl() locally (as sometimes the HTCondor system is slammed)
StepwiseQTL:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_07_stepwiseqtl.local.R
		
#FitQTL Makes the QTL and fits the QTL model derived from stepwiseqtl()
FitQTL:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_08_makefitqtls.R

#Collate all QTL info into one file
CollateQTL:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_09_collate.R

#Generate QTL table
GenQTLTable:
		R CMD BATCH '--args workflow="$(WORKFLOW)" qtl_scan_method="$(QTL_METHOD)" num_top_qtls="$(NUM_TOP)"' qtl_pipeline_10_01_genQTLtables.R

GenIntTable:
		R CMD BATCH '--args workflow="$(WORKFLOW)" num_top_qtls="$(NUM_TOP)"' qtl_pipeline_10_02_genInteractionTables.R

#Collate all QTL info into one file
GenConsensus:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_10_genconsensus.R

#Generate Circos Plot
GenCircos:
		R CMD BATCH '--args workflow="$(WORKFLOW)"' qtl_pipeline_11_gencircos.R

#Generate QTL table
GenEffects:
		R CMD BATCH '--args workflow="$(WORKFLOW)" num_top_qtls="$(NUM_TOP)"' qtl_pipeline_12_geneffects.R

GenIntEffects:
		R CMD BATCH '--args workflow="$(WORKFLOW)" collate_effects="$(COLLATE_EFFECTS)" num_top_qtls="$(NUM_TOP)"' qtl_pipeline_13_geninteffects.R

#Copy all circos files and their dependencies to a 'sandbox' folder where all the necessary files are available for rendering the circos plots from a web browser
#NOTE: Must specify SANDBOX on command-line to give destination directory to copy to.
GenCircosSandbox:
		mkdir -p $(SANDBOX)
		cp -r ../../Workflows/templates/traits/plots/circos/package*.json  $(SANDBOX)/
		cp -r ../../Workflows/templates/traits/plots/circos/drawcircos.*  $(SANDBOX)/
		cp -r ../../Workflows/templates/traits/plots/circos/*.css $(SANDBOX)/
		cp -r ../../Workflows/templates/traits/plots/circos/fonts $(SANDBOX)/
		cp -r $(WORKFLOW)/configs $(SANDBOX)/
	cd $(SANDBOX) && npm install

LaunchCircosSandbox:
		cd $(SANDBOX) && npm run -t start &
		sleep 5
		cd $(SANDBOX) && npm run -t launch
		
GenReactSandbox:
		rm -rf $(SANDBOX)
		cp -r ../../Workflows/templates/traits/plots/circos/circos-viewer  $(SANDBOX)
		cp -r $(WORKFLOW)/configs $(SANDBOX)/public/
		cp -r $(WORKFLOW)/traits/blups_collated.wide.json $(SANDBOX)/public/configs/
		cp -r $(WORKFLOW)/traits/plots/corrplot.$(POP_Circos).csv $(SANDBOX)/public/configs/corrplot.csv
		cp -r $(WORKFLOW)/traits/effects_collated.json $(SANDBOX)/public/configs/
		cp -r $(WORKFLOW)/traits/lod_profiles.json $(SANDBOX)/public/configs/
		cd $(SANDBOX) && npm install

GenReactSandboxConfigs:
		cp -r ../../Workflows/templates/traits/plots/circos/circos-viewer/public  $(SANDBOX)
		cp -r $(WORKFLOW)/configs $(SANDBOX)/public/
		cp -r $(WORKFLOW)/traits/blups_collated.wide.json $(SANDBOX)/public/configs/
		cp -r $(WORKFLOW)/traits/plots/corrplot.$(POP_Circos).csv $(SANDBOX)/public/configs/corrplot.csv
		cp -r $(WORKFLOW)/traits/effects_collated.json $(SANDBOX)/public/configs/
		cp -r $(WORKFLOW)/traits/lod_profiles.json $(SANDBOX)/public/configs/
		cd $(SANDBOX) && npm install

LaunchReactSandbox:
		cd $(SANDBOX) && npm start

BuildReactSandbox:
		cd $(SANDBOX) && npm run build
		cp $(SANDBOX)/build/index.html $(SANDBOX)/build/index.orig.html


CopyReactSandbox:
		cd $(SANDBOX)/build && sed -e 's/\/\(static\)/\1/g' -e 's/\/\(logo192.png\)/\1/g'  -e 's/\/\(favicon.ico\)/\1/g' index.orig.html > index.html
		#cd $(SANDBOX)/build && sed --in-place=.bak -e 's/\/(static\/js\/)/\0/g' -e 's/\/static\/css\//\/$(CIRCOS_FOLDER)\0/g' -e 's/logo192.png/$(CIRCOS_FOLDER)\/\0/g'  -e 's/favicon.ico/$(CIRCOS_FOLDER)\/\0/g' index.html
		ssh -X -Y -p 7262 maule@vaccinium.vcru.wisc.edu rm -fr /var/www/html/$(CIRCOS_FOLDER)
		scp -r -P 7262 $(SANDBOX)/build maule@vaccinium.vcru.wisc.edu:/var/www/html/$(CIRCOS_FOLDER)

clean:
		rm *.Rout
		make -C $(HTCONDOR) $@
