#!/bin/bash
#PBS -q production
#PBS -N wcs-fsc2
#PBS -l select=1:ncpus=2
#PBS -l place=free
#PBS -V

# change to the working directory 
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
WCSWORKDIR=/scratch/isaac.overcast/WCS-fsc2/hickerlab-repository/WCS-fsc2/withhybrids
cd $WCSWORKDIR
export PATH=./:$PATH
echo $PATH

export OMP_NUM_THREADS=2
#/scratch/isaac.overcast/WCS-fsc2/withhybrids/run_fsc.sh -p nomigration -t -m > wcs-fsc.out 2>&1
./run_fsc.sh -p secondary-asym -m -n 5 > wcs-fsc.out 2>&1 $EXECDIR/wcs.out
#./run_fsc.sh -p symmetricalmigration -m -n 50 > wcs-fsc.out 2>&1 $EXECDIR/wcs.out

echo ">>>> Hoggy!..."
