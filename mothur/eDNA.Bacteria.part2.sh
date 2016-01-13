#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=250gb,walltime=24:00:00
#PBS -M lennonj@indiana.edu,mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2016_eDNA
module load gcc/4.9.2
module load mothur/1.36.1
mothur eDNA.Bacteria.part2.Batch
