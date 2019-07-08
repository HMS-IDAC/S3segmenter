#!/bin/bash
for fileNum in {1..6}
do
sbatch S3segmenter.sh $fileNum
done

