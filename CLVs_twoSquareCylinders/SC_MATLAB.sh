#!/bin/sh 
#PBS -l walltime=72:00:00,select=1:ncpus=8:mpiprocs=8:mem=1000gb 

module load tools/prod
module load MATLAB/2023a_Update_3

#cd $PBS_O_WORKDIR

cp -r $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/U_D_collection.txt $TMPDIR 
cp -r $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/R_full.txt $TMPDIR 
cp -r $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/SC_CLV_orig.m $TMPDIR

cd $TMPDIR

matlab -nosplash -nodisplay -batch SC_CLV_orig -logfile output.log

cp -r $TMPDIR/CLV1x.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV1y.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV2x.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV2y.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV3x.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV3y.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
#cp -r $TMPDIR/CLV4x.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
#cp -r $TMPDIR/CLV4y.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
#cp -r $TMPDIR/CLV5x.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
#cp -r $TMPDIR/CLV5y.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV6x.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV
cp -r $TMPDIR/CLV6y.txt $HOME/PhD/OpenFOAM/SC_HPC/run_stept03_32cores_4/run_plotCLV

