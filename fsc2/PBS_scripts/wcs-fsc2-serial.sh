#!/bin/bash
#PBS -q production
#PBS -N wcs-fsc2
#PBS -l select=1:ncpus=1
#PBS -l place=free
#PBS -V

# change to the working directory 
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
WCSWORKDIR=/scratch/isaac.overcast/WCS-fsc2/nohybrids
cd $WCSWORKDIR
export PATH=./:$PATH
echo $PATH
#/scratch/isaac.overcast/WCS-fsc2/withhybrids/run_fsc.sh -p nomigration -t -m > wcs-fsc.out 2>&1
./run_fsc.sh -p nomigration -m > wcs-fsc.out 2>&1 $EXECDIR/wcs.out

echo ">>>> Begin <job_name> Run ..."
