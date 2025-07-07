#!/bin/sh 
#PBS -l walltime=72:00:00,select=1:ncpus=32:mpiprocs=32:mem=1000gb 

module load tools/prod
module load OpenFOAM/10-foss-2022a-20230119
source $FOAM_BASH

#copy the solver files to the temporary directory 
cp -r $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4 $TMPDIR 

#compile solver over 1 processor 
wmake $TMPDIR/run_stept03_32cores_4/solver/QRcode/
wmake $TMPDIR/run_stept03_32cores_4/solver/SC_Evolver/

cd $PBS_O_WORKDIR

  log_file="loop_times.log"
# Capture the start time in seconds since epoch
  start_time=$(date +%s)
  start_time_human=$(date +'%Y-%m-%d %H:%M:%S') 

for i in {1..4600}
do
  decomposePar -force
  mpirun -np 32 SC_Evolver -parallel 
  reconstructPar            
  QRcode  
done

# Capture the end time in seconds since epoch
  end_time=$(date +%s)
  end_time_human=$(date +'%Y-%m-%d %H:%M:%S')
  # Calculate the duration of the iteration
  duration=$((end_time - start_time))
  # Write the information to the log file
  echo "Iteration $i took $duration seconds" >> $log_file
  echo "---------------------------" >> $log_file
  echo 
