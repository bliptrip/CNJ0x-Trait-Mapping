#!/opt/local/bin/python

import os
import sys
import shutil
import subprocess
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Build submit files and dependent folders (for stepwiseqtl) using the mmers/traits/subtraits specified on the comamnd-line.")
    parser.add_argument('-m', '--mmers', action='store', required=True, help="Specific mixed models to use (years, all-years)")
    parser.add_argument('-t', '--traits', action='store', required=True, help="Traits (groupings) to run on -- this represents if multivariate analysis was done on a group of traits vs. not.")
    parser.add_argument('-l', '--subtraits', action='store', required=True, help="Subtraits to iterate over within traits.")
    parser.add_argument('-r', '--r_version', action='store', required=True, help="Version of R used.")
    parser.add_argument('--template', action='store', required=True, help="Template submit file to copy to subfolders.")
    parser.add_argument('--tarfile', action='store', required=True, help="Tar file to add generated submit data to.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

def generate_submit_data(masterfold, cluster_number, mmer, trait_grouping, subtrait, r_version, template, tarfile):
    #Generate a trait master folder with mmer__trait_grouping__subtrait designator
    #Loop through the folder and generate subfolders for each process
    subfold = masterfold + '/' + str(cluster_number)
    os.mkdir(subfold)
    mmer_fd = open(subfold + '/input.mmer', 'w+')
    mmer_fd.write(mmer)
    mmer_fd.close()
    traitg_fd = open(subfold + '/input.traitg', 'w+')
    traitg_fd.write(trait_grouping)
    traitg_fd.close()
    subtrait_fd = open(subfold + '/input.subtrait', 'w+')
    subtrait_fd.write(subtrait)
    subtrait_fd.close()

if __name__ == '__main__':
    parsed                 = parse_args()
    mmers                  = parsed.mmers
    traits                 = parsed.traits
    subtraits              = parsed.subtraits
    r_version              = parsed.r_version
    template               = parsed.template
    tarfile				   = parsed.tarfile

    #Copy the template HTCondor submit file and edit
    masterfold = sys.argv[0] + "--sandbox"
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

    #Generate a submit file for each mmer-trait-subtrait under consideration.
    traits_array = traits.split(';')
    subtraits_array = subtraits.split(';')

    total_clusters_needed  = 0
    for mmer in mmers.split(';'):
        for i in range(0,len(traits_array)):
            trait_grouping     = traits_array[i]
            subtrait_grouping  = subtraits_array[i]
            for j,subtrait in enumerate(subtrait_grouping.split(',')):
                generate_submit_data(masterfold, total_clusters_needed, mmer,trait_grouping,subtrait,r_version,template,tarfile)
                total_clusters_needed = total_clusters_needed + 1

    #Finish editing the queue to reflect the number of clusters
    cmd = "sed -E -i \"\" -e 's/^(queue[[:space:]]*)[[:digit:]]+/\\1" + str(total_clusters_needed) + "/g' " + subfile
    subprocess.call(cmd, shell=True)

    #Now add the sandbox foler to the tar file and remove it
    subprocess.call("tar -r -C %s -f %s ./" % (masterfold, tarfile), shell=True)
    #Now that the master folder has been added to the tarball, remove the folder from current path.
    shutil.rmtree(masterfold)
