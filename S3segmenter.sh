#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-1:20
#SBATCH -p short
#SBATCH --mem=64000
module load matlab/2018b
matlab -nodesktop -r "O2batchS3segmenterWrapperR('/n/scratch2/cy101/Denis/MedDev','HPC','true','fileNum',$1,'TissueMaskChan',[1 2 38],'logSigma',[3 60],'mask','tissue','CytoMaskChan',[28])"
