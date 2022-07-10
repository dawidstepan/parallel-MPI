#!/bin/bash
# Number of cores
#SBATCH -c 16
# Runtime of this jobs is less then 40 minutes
#            (hh:mm:ss)
#SBATCH --time=00:40:00
# Clear the environment
module purge > /dev/null 2>&1

rm -f log.txt

for s in $(seq 2 1 16);
do
	mpirun -n $s ./parallel 8192 100 0
done
