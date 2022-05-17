#!/bin/bash
#set -x
 
installdir=`pwd`/anaconda
 
###################################
# Python 3 users, run this line:
##################################
bash /collab/usr/gapps/python/blueos_3_ppc64le_ib/conda/Anaconda3-2021.05-Linux-ppc64le.sh -b -f -p $installdir
 
source $installdir/bin/activate
 
conda config --add default_channels https://repo.anaconda.com/pkgs/main
conda config --add default_channels https://repo.anaconda.com/pkgs/r
 
##################################
# https://www.ibm.com/support/knowledgecenter/SS5SF7_1.6.0/navigation/pai_install.html
##################################
 
conda config --set ssl_verify /etc/pki/tls/cert.pem
 
#conda config --prepend channels \
#https://public.dhe.ibm.com/ibmdl/export/pub/software/server/ibm-ai/conda/
 
export IBM_POWERAI_LICENSE_ACCEPT=yes
conda config --prepend channels https://public.dhe.ibm.com/ibmdl/export/pub/software/server/ibm-ai/conda-early-access/
conda create -n tfvenv python=3.7
conda activate ~/.conda/envs/tfvenv
conda install -y scipy opencv absl-py cudatoolkit-dev tensorflow-gpu
