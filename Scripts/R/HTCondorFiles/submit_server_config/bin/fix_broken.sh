#!/bin/sh

~/bin/find_all_empty_perms2d.sh | while read size file; do 
	file_path=${file%/*}
	id=$(basename $file_path)
	cp qtl_pipeline_05_runperms2D.fixmissing.sub ${file_path}
	cd ${file_path}
	sed -i -e "s/BROKEN_ID=0/BROKEN_ID=${id}/" qtl_pipeline_05_runperms2D.fixmissing.sub
	condor_submit qtl_pipeline_05_runperms2D.fixmissing.sub
	cd -
done


