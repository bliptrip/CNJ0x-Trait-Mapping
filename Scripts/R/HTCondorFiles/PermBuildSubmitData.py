#!/usr/bin/env python

import csv
import os
import pandas as pd
import sys
import shutil
import subprocess
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Build submit files and dependent folders using model-traits.cfg.csv configuration file provided from the workflow.")
    parser.add_argument('-w', '--workflow', action='store', default='1', help="Location of workflow folder to access configuration file, etc.")
    parser.add_argument('--total_perms_per_trait', action='store', type=int, default=1000, help="Total number of permutations to run on a given trait.")
    parser.add_argument('--num_cluster_per_trait', action='store', type=int, default=10, help="Total number of clusters to queue up for each set of multi-traits analyzed.")
    parser.add_argument('-s', '--seed', action='store', type=int, required=True, help="Initialization seed.")
    parser.add_argument('-r', '--r_version', action='store', required=True, help="Version of R used.")
    parser.add_argument('--template', action='store', required=True, help="Template submit file to copy to subfolders.")
    parser.add_argument('--tarfile', action='store', required=True, help="Tar file to add generate submit data to.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

def generate_submit_data(num_clusters, num_perms, seed, workflow, model, trait, qtlmethod, r_version, template, tarfile):
    masterfold = model + '--' + trait 
    srcfold = workflow + "/traits/" + masterfold
    cross_df = pd.read_csv(srcfold + "/cross.csv")
    if cross_df[trait].any(skipna = True): #Make sure that there is at least one non-zero BLUP value -- otherwise, mmer failed to estimate BLUPs for this trait
        os.mkdir(masterfold)
        #Copy the template submit file to this folder
        subfile = re.sub(r'(.+)\.template\.sub', r'\1.sub', template)
        subfile = masterfold + '/' + subfile
        shutil.copy(template, subfile)
        #Edit the sub file to match configuration settings
        cmd = "/usr/local/bin/gsed -i -e 's/{R_VERSION}/" + r_version + "/g' " + subfile
        print(cmd)
        subprocess.call(cmd, shell=True)
        cmd = "/usr/local/bin/gsed -E -i -e 's/^(queue[[:space:]]*)[[:digit:]]+/\\1" + str(num_clusters) + "/g' " + subfile
        print(cmd)
        subprocess.call(cmd, shell=True)
        #Copy the appropriate cross.csv into the current folder.
        shutil.copy(srcfold + "/cross.csv", masterfold + "/cross.csv")
        shutil.copy(srcfold + "/cross.rds", masterfold + "/cross.rds")
        #Loop through the folder and generate subfolders for each process
        for i in range(0,num_clusters):
            subfold = masterfold + '/' + str(i)
            os.mkdir(subfold)
            seed_fd   = open(subfold + '/input.seed', 'w+')
            seed_fd.write(str(seed + i))
            seed_fd.close()
            nperms_fd = open(subfold + '/input.nperms', 'w+')
            nperms_fd.write(str(num_perms))
            nperms_fd.close()
            model_fd = open(subfold + '/input.model', 'w+')
            model_fd.write(model)
            model_fd.close()
            traitg_fd = open(subfold + '/input.trait', 'w+')
            traitg_fd.write(trait)
            traitg_fd.close()
            qtlmethod_fd = open(subfold + '/input.qtlmethod', 'w+')
            qtlmethod_fd.write(qtlmethod)
            qtlmethod_fd.close()
        subprocess.call("tar -rf %s %s" % (tarfile,masterfold), shell=True)
        #Now that the master folder has been added to the tarball, remove the folder from current path.
        shutil.rmtree(masterfold)


if __name__ == '__main__':
    parsed                 = parse_args()
    total_perms_per_trait  = parsed.total_perms_per_trait
    num_cluster_per_trait  = parsed.num_cluster_per_trait
    num_perms_per_trait    = total_perms_per_trait/num_cluster_per_trait
    initial_seed           = parsed.seed
    workflow               = parsed.workflow
    r_version              = parsed.r_version
    template               = parsed.template
    tarfile				   = parsed.tarfile

    #Read in configuration file.
    #These should be read in from workflow/configs/model.cfg.R
    exec(open(workflow+'/configs/model.cfg', 'r').read())

    #Loop through entries in the config file and generate a submit file for each entry
    #Don't forget to look for masked items and ignore them.
    #Generate a submit file for each model-trait under consideration.
    with open(workflow+'/configs/model-traits.cfg.csv', 'rU') as csvfile:
        config_reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        i = 0
        for model_trait in config_reader:
            model = model_trait["model"]
            trait = model_trait["trait"]
            mask = model_trait["mask"]
            if( (not mask) or (mask != "TRUE") ):
                generate_submit_data(num_cluster_per_trait,num_perms_per_trait,initial_seed+(i),workflow,model,trait,qtl_method,r_version,template,tarfile)
            i += 1
        csvfile.close()
