#Makefile for building/packaging data & scripts for HTCondor submission

USERNAME ?= maule2
SUBMIT_SERVER ?= submit-1.chtc.wisc.edu
include $(WORKFLOW)/configs/perms.mk

.PHONY: Rbuild PermBuild

#Rbuild tars up the packages necessary for installing an R-build and then transferring to submit server.  Untar Rpackage.tgz, submit an interactive session using
#condor_submit -i Rpackage.sub
#Comment-out any packages you don't want to install if they're already in the tarball.
Rbuild:
		gsed -e 's/{R_VERSION}/$(R_VERSION)/g' Rpackage.install.sh.template > Rpackage.install.sh
		gsed -e 's/{R_VERSION}/$(R_VERSION)/g' Rpackage.sub.template > Rpackage.sub
		chmod u+x Rpackage.install.sh
		ssh $(USERNAME)@$(SUBMIT_SERVER) 'mkdir -p projects/'
		scp ../packagedeps.chtc.R Rpackage.install.sh Rpackage.sub $(USERNAME)@$(SUBMIT_SERVER):~/projects/
		ssh $(USERNAME)@$(SUBMIT_SERVER) 'cd ~/projects; condor_submit -batch-name "$(USERNAME)_Rpackage" Rpackage.sub'

#Run the following after checking that the Rbuild condor task completed successfully.
Rdownload:
		scp $(USERNAME)@$(SUBMIT_SERVER):~/projects/packages-R$(R_VERSION).tar.gz .
		
		
#Permbuild is meant for building a tarball of necessary scripts/etc. needed for running single-QTL (scanone) permutation test on cluster.
PermBuild:
		#Untar the R-package and only extract binaries for transferring to each cluster so as to reduce size of package.
		tar -cf qtl_pipeline_02_runperms.tar ../qtl_pipeline_02_runperms.R qtl_pipeline_02_runperms.sh qtl_pipeline_02_runperms.all.sh packages-R$(R_VERSION).tar.gz
		./PermBuildSubmitData.py --total_perms_per_trait $(TOTAL_PERMS_PER_TRAIT) --num_cluster_per_trait $(NUM_CLUSTERS_PER_TRAIT) -s $(START_SEED) -w "$(WORKFLOW)" -r $(R_VERSION) --template qtl_pipeline_02_runperms.template.sub --tarfile qtl_pipeline_02_runperms.tar 
		gzip -f -9 qtl_pipeline_02_runperms.tar 
		scp qtl_pipeline_02_runperms.tar.gz $(USERNAME)@$(SUBMIT_SERVER):~/projects/

#RunPerms will extract the permutations folder and submit the condor jobs.
RunPerms:
		ssh $(USERNAME)@$(SUBMIT_SERVER) 'mkdir -p projects/CNJ0x-Trait-Mapping-runperms; cd projects/CNJ0x-Trait-Mapping-runperms; rm -rf *; tar -xzf ../qtl_pipeline_02_runperms.tar.gz; ./qtl_pipeline_02_runperms.all.sh'

#PullPerms will extract the calculated permutations from the submit server
PullPerms:
		ssh $(USERNAME)@$(SUBMIT_SERVER) 'cd projects/; tar -C CNJ0x-Trait-Mapping-runperms -czf CNJ0x-Trait-Mapping-runperms.tar.gz ./'
		scp $(USERNAME)@$(SUBMIT_SERVER):~/projects/CNJ0x-Trait-Mapping-runperms.tar.gz $(WORKFLOW)/
		mkdir -p $(WORKFLOW)/perms
		tar -C $(WORKFLOW)/perms -xzf $(WORKFLOW)/CNJ0x-Trait-Mapping-runperms.tar.gz
		rm $(WORKFLOW)/CNJ0x-Trait-Mapping-runperms.tar.gz


#Perm2DBuild is meant for building a tarball of necessary scripts/etc. needed for running 2D-QTL (scantwo) permutation tests on cluster.
Perm2DBuild:
		#Untar the R-package and only extract binaries for transferring to each cluster so as to reduce size of package.
		tar -cf qtl_pipeline_05_runperms2D.tar ../qtl_pipeline_05_runperms2D.R qtl_pipeline_05_runperms2D.sh qtl_pipeline_05_runperms2D.all.sh packages-R$(R_VERSION).tar.gz
		./PermBuildSubmitData.py --total_perms_per_trait $(TOTAL_PERMS_PER_TRAIT) --num_cluster_per_trait $(NUM_CLUSTERS_PER_TRAIT) -s $(START_SEED) -w "$(WORKFLOW)" -r $(R_VERSION) --template qtl_pipeline_05_runperms2D.template.sub --tarfile qtl_pipeline_05_runperms2D.tar
		gzip -f -9 qtl_pipeline_05_runperms2D.tar 
		scp qtl_pipeline_05_runperms2D.tar.gz $(USERNAME)@$(SUBMIT_SERVER):~/

#PullPerms2D will extract the calculated 2D permutations from the submit server
PullPerms2D:
		ssh $(USERNAME)@$(SUBMIT_SERVER) 'cd projects/; tar -C CNJ0x-Trait-Mapping-runperms2D -czf CNJ0x-Trait-Mapping-runperms2D.tar.gz ./'
		scp $(USERNAME)@$(SUBMIT_SERVER):~/projects/CNJ0x-Trait-Mapping-runperms2D.tar.gz $(WORKFLOW)/
		mkdir -p $(WORKFLOW)/perms2D
		tar -C $(WORKFLOW)/perms2D -xzf $(WORKFLOW)/CNJ0x-Trait-Mapping-runperms2D.tar.gz
		rm $(WORKFLOW)/CNJ0x-Trait-Mapping-runperms2D.tar.gz

#RunPerms2D will extract the 2D permutations folder and submit the condor jobs.
RunPerms2D:
		ssh $(USERNAME)@$(SUBMIT_SERVER) 'mkdir -p projects/CNJ0x-Trait-Mapping-runperms2D; cd projects/CNJ0x-Trait-Mapping-runperms2D; rm -rf *; tar -xzf ../../qtl_pipeline_05_runperms2D.tar.gz; ./qtl_pipeline_05_runperms2D.all.sh'

StepwiseQTLBuild::
		#Untar the R-package and only extract binaries for transferring to each cluster so as to reduce size of package.
		if [[ ! -f $(R_VERSION)-chtc.tar.gz ]]; then \
				tar -xzf $(R_VERSION).tar.gz ; \
				tar -cz -C $(R_VERSION)/lib64/ -f $(R_VERSION)-chtc.tar.gz R ; \
				rm -rf $(R_VERSION) ; \
		fi
		#Build the support datafiles archive
		tar -cf qtl_pipeline_07_stepwiseqtl.tar ../qtl_pipeline_07_stepwiseqtl.R qtl_pipeline_07_stepwiseqtl.sh $(R_VERSION)-chtc.tar.gz
		./StepwiseQTLBuildSubmitData.py -w "$(WORKFLOW)" -r $(R_VERSION) --template qtl_pipeline_07_stepwiseqtl.template.sub --tarfile qtl_pipeline_07_stepwiseqtl.tar 
		gzip -f -9 qtl_pipeline_07_stepwiseqtl.tar 
		scp qtl_pipeline_07_stepwiseqtl.tar.gz $(USERNAME)@$(SUBMIT_SERVER):~/

clean:
		rm -rf qtl_pipeline_02_runperms.tar.gz 
		rm -rf Rpackage.tgz
		rm -rf $(R_VERSION)-chtc.tar.gz
		rm -f *.bak*
