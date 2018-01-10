#!/bin/sh

#This script simply runs all of the traits using separate condor_submit requests,
#but eventually will take a parameter that specifies that batch subsets can be run so as
#not to slam the system.
#
#NOTE: This must be run from folder that contains all mmer-trait_grouping-subtrait subfolders in order to
#work.

for file in *; do
   if [ -d $file ]; then
		cd $file; condor_submit -batch-name "AM_RUNPERMS2D" qtl_pipeline_05_runperms2D.sub; cd ..
   fi
done
