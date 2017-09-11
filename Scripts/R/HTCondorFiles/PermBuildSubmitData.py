#!/opt/local/bin/python

import os
import sys
import shutil
import subprocess
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Build submit files and dependent folders using the mmers/traits/subtraits specified on the comamnd-line.")
    parser.add_argument('--total_perms_per_trait', action='store', type=int, default=1000, help="Total number of permutations to run on a given trait.")
    parser.add_argument('--num_cluster_per_trait', action='store', type=int, default=10, help="Total number of clusters to queue up for each subtrait analyzed.")
    parser.add_argument('-s', '--seed', action='store', type=int, required=True, help="Initialization seed.")
    parser.add_argument('-m', '--mmers', action='store', required=True, help="Specific mixed models to use (years, all-years)")
    parser.add_argument('-t', '--traits', action='store', required=True, help="Traits (groupings) to run on -- this represents if multivariate analysis was done on a group of traits vs. not.")
    parser.add_argument('-l', '--subtraits', action='store', required=True, help="Subtraits to iterate over within traits.")
    parser.add_argument('-r', '--r_version', action='store', required=True, help="Version of R used.")
    parser.add_argument('--template', action='store', required=True, help="Template submit file to copy to subfolders.")
    parser.add_argument('--tarfile', action='store', required=True, help="Tar file to add generate submit data to.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

def generate_submit_data(num_clusters, num_perms, seed, mmer, trait_grouping, subtrait, r_version, template, tarfile):
    #Generate a trait master folder with mmer__trait_grouping__subtrait designator
    masterfold = mmer + '--' + '__'.join(trait_grouping.split(',')) + '--' + subtrait
    os.mkdir(masterfold)
    #Copy the template submit file to this folder
    subfile = re.sub(r'(.+)\.template\.sub', r'\1.sub', template)
    subfile = masterfold + '/' + subfile
    shutil.copy(template, subfile)
    #Edit the sub file to match configuration settings
    cmd = "sed -E -i \"\" -e 's/^(arguments[[:space:]]*=.*)/\\1 " + r_version + "-chtc.tar.gz/g' " + subfile
    subprocess.call(cmd, shell=True)
    cmd = "sed -E -i \"\" -e 's/^(transfer_input_files[[:space:]]*=.*)/\\1,..\/..\/" + r_version + "-chtc.tar.gz/g\' " + subfile
    subprocess.call(cmd, shell=True)
    cmd = "sed -E -i \"\" -e 's/^(queue[[:space:]]*)[[:digit:]]+/\\1" + str(num_clusters) + "/g' " + subfile
    subprocess.call(cmd, shell=True)
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
        mmer_fd = open(subfold + '/input.mmer', 'w+')
        mmer_fd.write(mmer)
        mmer_fd.close()
        traitg_fd = open(subfold + '/input.traitg', 'w+')
        traitg_fd.write(trait_grouping)
        traitg_fd.close()
        subtrait_fd = open(subfold + '/input.subtrait', 'w+')
        subtrait_fd.write(subtrait)
        subtrait_fd.close()
    subprocess.call("tar -rf %s %s" % (tarfile,masterfold), shell=True)
    #Now that the master folder has been added to the tarball, remove the folder from current path.
    shutil.rmtree(masterfold)


if __name__ == '__main__':
    parsed                 = parse_args()
    total_perms_per_trait  = parsed.total_perms_per_trait
    num_cluster_per_trait  = parsed.num_cluster_per_trait
    num_perms_per_trait    = total_perms_per_trait/num_cluster_per_trait
    initial_seed           = parsed.seed
    mmers                  = parsed.mmers
    traits                 = parsed.traits
    subtraits              = parsed.subtraits
    r_version              = parsed.r_version
    template               = parsed.template
    tarfile				   = parsed.tarfile

    #Generate a submit file for each mmer-trait-subtrait under consideration.
    traits_array = traits.split(';')
    subtraits_array = subtraits.split(';')
    for mmer in mmers.split(';'):
        for i in range(0,len(traits_array)):
            trait_grouping     = traits_array[i]
            subtrait_grouping  = subtraits_array[i]
            for j,subtrait in enumerate(subtrait_grouping.split(',')):
                generate_submit_data(num_cluster_per_trait,num_perms_per_trait,initial_seed+(i*j),mmer,trait_grouping,subtrait,r_version,template,tarfile)
