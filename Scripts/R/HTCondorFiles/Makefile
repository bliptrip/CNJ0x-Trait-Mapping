#Makefile for building/packaging data & scripts for HTCondor submission

.PHONY: Rbuild PermBuild

QTL_DIR=../../../Data/genetic_data/DerivedData/qtl
RDS_DIR=../../../Data/phenotypic\ data/DerivedData/cleanup_data.R.output/Robjs
TOTAL_MMERS=4 #4 year
TOTAL_PERMS_PER_TRAIT=10
NUM_CLUSTERS_PER_TRAIT=1
#For now I have to harcode this, as I can't figure out how to get the calculation to work in this makefile.  Just divide TOTAL_PERMS_PER_TRAIT/NUM_CLUSTERS_PER_TRAIT
NUM_PERMS_PER_CLUSTER=10
START_SEED=54955149
#Specify the following based on teh versionof R you want to build/use
R_VERSION=R-3.4.1
#Set the following to determine which analysis(es) and trait(s) to run permutations on.  If multiple are to be selected, then separate them with commas.
#Remember to wrap the following variables in quotes!
MMER="2011"
#If dealing with multiple traits, comma-separate them
TRAITS="total_berry_weight"
#Specify a subtrait if dealing with a multivariate trait analysis
SUBTRAIT=


#Rbuild tars up the packages necessary for installing an R-build and then transferring to submit server.  Untar Rpackage.tgz, submit an interactive session using
#condor_submit -i Rpackage.sub
#Comment-out any packages you don't want to install if they're already in the tarball.
Rbuild:
		sed -E -i.bak -e 's/^R_VERSION\s*=.*/R_VERSION=$(R_VERSION)/g' Rpackage.install.sh
		tar -czf Rpackage.tgz packagedeps.R Rpackage.install.sh Rpackage.sub $(R_VERSION).tar.gz
		scp Rpackage.tgz maule2@submit-5.chtc.wisc.edu:~/
		
#Permbuild is meant for building a tarball of necessary scripts/etc. needed for running single-QTL (scanone) permutation test on cluster.
PermBuild:
		#Untar the R-package and only extract binaries for transferring to each cluster so as to reduce size of package.
		if [[ ! -f $(R_VERSION)-chtc.tar.gz ]]; then \
				tar -xzf $(R_VERSION).tar.gz ; \
				tar -czf $(R_VERSION)-chtc.tar.gz $(R_VERSION)/lib64/R ; \
				rm -rf $(R_VERSION) ; \
		fi
		#Modify the submit file and the script file to reflect the number of permutations and the number of queued jobs needed
		sed -E -i.bak -e 's/^arguments = .*/arguments = input.seed $(NUM_PERMS_PER_CLUSTER) $(R_VERSION)-chtc.tar.gz $(MMER) $(TRAITS) $(SUBTRAIT)/g' qtl_pipeline_02_runperms.sub
		sed -E -i.bak -e 's/^queue [[:digit:]]+/queue $(NUM_CLUSTERS_PER_TRAIT)/g' qtl_pipeline_02_runperms.sub
		#Add in the appropriate version of R
		sed -E -i.bak -e 's/^transfer_input_files\s*=.*/transfer_input_files=$(R_VERSION)-chtc.tar.gz,$$(Process)\/input.seed,qtl_pipeline_02_runperms.R,qtl_pipeline_02_runperms.datafiles.tar.gz/g' qtl_pipeline_02_runperms.sub
		#Build the support datafiles archive
		tar -czf qtl_pipeline_02_runperms.datafiles.tar.gz $(QTL_DIR) $(RDS_DIR)/cnjpop.mmer2.p2.rds
		tar -cf qtl_pipeline_02_runperms.tar ../qtl_pipeline_02_runperms.R qtl_pipeline_02_runperms.sh qtl_pipeline_02_runperms.sub $(R_VERSION)-chtc.tar.gz qtl_pipeline_02_runperms.datafiles.tar.gz
		number=0; while [[ $$number -lt $(NUM_CLUSTERS_PER_TRAIT) ]]; do \
				mkdir $$number ; \
				((process_seed=$(START_SEED) + number)) ; \
				echo "$$process_seed" > $$number/input.seed ; \
				tar -rf qtl_pipeline_02_runperms.tar $$number/ ; \
				rm -rf $$number/ ; \
				((number=number+1)) ; \
		done;
		gzip -9 qtl_pipeline_02_runperms.tar 
		scp qtl_pipeline_02_runperms.tar.gz maule2@submit-5.chtc.wisc.edu:~/


clean:
		rm -rf Rpackage.tgz
		rm -rf $(R_VERSION)-chtc.tar.gz
		rm -f *.bak*
