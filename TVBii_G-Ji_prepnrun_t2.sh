#!/bin/bash

#PBS -N TVBii
#PBS -m abe
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -l vmem=30GB


# ==============================================================================
# DATA:
# ==============================================================================

if [ -z "$PBS_ARRAYID" ]
then
	echo 'ERROR: $PBS_ARRAYID is not set, submit job(s) using "qsub -t <array expression>"'
	exit 1
fi

# Directory where inputs are located (required only when copying from non-unix pc to hpc)
#dos2unix $HOMEDIR/input.txt --quiet

# Make input variable with subIDs
export subID=`sed -n "${PBS_ARRAYID}p" $HOMEDIR/input.txt`
echo "subID: $subID"

# Define TVBii folder
subFolder="$HOMEDIR/subjects/${subID}"
cd $subFolder


# First get TR 
# TODO: Make a simple text file for every subject, with 2.1 (only for PAT01T2) 
# or 2.4 (for all other subjects)
tr=$( cat $subFolder/tr.txt ) 


# Generate parameter files (depending on TR)
LANG=eng_US
it=0

if [ $tr = 2.1 ]; then 
	trm=2100
	simlen=420000
elif [ $tr = 2.4 ]; then 
	trm=2400
	simlen=480000
fi

for G in $(seq 0.01 0.015 3.00)
do
	it=$(($it+1))
	echo "68 $G 0.15 1.40 1.0 0.01 $simlen 10000 $trm 999999.99 42" > param_set_$it
done


# Then copy the TVBii executable
cp $HOMEDIR/tvbii_linux $subFolder/tvbii_linux


# Run TVBii for all parameter files, for all thresholds and distance metrics
# --> for post-op analyses, sufficient to run on thrA, distLog

# thrA - distLog
for sim in {1..200}
do
	echo "Processing iteration $sim"
	./tvbii_linux param_set_${sim} "${subID}"_scale68_thrA_distLog
done

mv $subFolder/output $subFolder/output_thrA_distLog
rm ${subFolder}/param_set_*

