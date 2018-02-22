#Makefile for building/packaging data & scripts for HTCondor submission

.PHONY: install GenCross PermMerge Perm2DMerge Rbuild PermBuild StepwiseQTLBuild StepwiseQTL

HTCONDOR=HTCondorFiles

export WORKFLOW ?= ../../Workflows/1
export RDS_DIR ?= "../../Data/phenotypic\ data/DerivedData/cleanup_data.R.output/Robjs"
export SCAN_DIR ?= ../../Analysis
export SCANONE_MERGE ?= $(SCAN_DIR)/runperms/cnjpop.scanone.p2.rds
export SCANTWO_MERGE ?= $(SCAN_DIR)/runperms2D/cnjpop.scantwo.p2.rds
#Set the following to "scanone" or "stepwiseqtl" depending on what dataset you want to use for generating the circos plot.
export CIRCOS_QTL_TYPE ?= stepwiseqtl

install:
		R CMD BATCH packagedeps.R

means:
		R CMD BATCH qtl_pipeline_00_1_means.R 

plotcors:
		R CMD BATCH qtl_pipeline_00_2_plotcors.R 

genmodels:
		R CMD BATCH '--args workflow="../../Workflows/1"' qtl_pipeline_01_genmodels.R

#Rbuild tars up the packages necessary for installing an R-build and then transferring to submit server.  Untar Rpackage.tgz, submit an interactive session using
#condor_submit -i Rpackage.sub
#Comment-out any packages you don't want to install if they're already in the tarball.
Rbuild:
		make -C $(HTCONDOR) QTL_DIR=../$(QTL_DIR) RDS_DIR=../$(RDS_DIR) $@
		
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

#Collate all QTL info into one file
GenConsensus:
		R CMD BATCH '--args workflow="$(WORKFLOW)" qtl_type="$(CIRCOS_QTL_TYPE)"' qtl_pipeline_10_genconsensus.R

#Generate Circos Plot
GenCircos:
		R CMD BATCH '--args workflow="$(WORKFLOW)" qtl_type="$(CIRCOS_QTL_TYPE)"' qtl_pipeline_11_gencircos.R

#Merge H2
MergeH2:
		echo "model,mtraits,trait,h2" > $(WORKFLOW)/traits/h2.csv; \
		IFS=','; for t in $(WORKFLOW)/traits/*--*; do \
				ft=$$t; \
				t=$${t##*/}; \
				model=$${t%--*}; \
				mtraits=$${t#*--}; \
				mtraitsb=$$(echo $$mtraits | sed -e 's/__/,/g'); \
				eval "declare -a mtraitsa=($$(echo $$mtraitsb))"; \
				numtraits=$${#mtraitsa[*]}; \
				if [ -f $${ft}/h2.csv ]; then \
						tail -n +2 $${ft}/h2.csv | while read trait h2; do \
								if [ "$$h2" != "0" ] && [ "$$h2" != "NA" ]; then \
										if [ $$numtraits -eq 1 ]; then \
												echo "$$model,\"$$mtraitsb\",$${mtraitsa[0]},$$h2"; \
										else \
												echo "$$model,\"$$mtraitsb\",$${trait//\"},$$h2"; \
										fi \
								fi \
						done >> $(WORKFLOW)/traits/h2.csv; \
				fi \
		done

test:
				IFS=','; tail -n +2 vcov.csv | while read m1 m2 m3 m4 m5 m6 m7; do \
						echo $$m1; \
				done; \


clean:
		rm *.Rout
		make -C $(HTCONDOR) $@
