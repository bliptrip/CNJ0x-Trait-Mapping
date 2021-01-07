#!/bin/sh
#Because of snafus I experienced when building R-3.4.1, I've setup a script to deal with my issues.
#
#Note: This is to be manually run from within HTCondor cluster node (not submit server!)
R_VERSION=R-3.1.2
tar -xzf ${R_VERSION}.tar.gz
cd ${R_VERSION}
#Need to specify JAVA_HOME or the compilation stage doesn't appear to work on their CENTOS clusters
./configure JAVA_HOME=/usr/lib/jvm/java-openjdk --prefix=$(pwd)
make
make install
cd ..
#Run the following script to install package dependencies
${R_VERSION}/lib64/R/bin/Rscript --vanilla packagedeps.R
#Modify the R_HOME_DIR as defined in the CHTC webpage
cp ${R_VERSION}/lib64/R/bin/R ${R_VERSION}/lib64/R/bin/R.bak 
awk -F= '{if ($0 ~ /^R_HOME_DIR=/) {print $1"=$(pwd)/R"} else { print $0 }}' < ${R_VERSION}/lib64/R/bin/R.bak > ${R_VERSION}/lib64/R/bin/R
rm ${R_VERSION}/lib64/R/bin/R.bak 
#Now repackage the built R-folder for re-use later
tar -czf ${R_VERSION}.tar.gz ${R_VERSION}
