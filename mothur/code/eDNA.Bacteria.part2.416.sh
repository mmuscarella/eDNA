#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=250gb,walltime=48:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2016_eDNA
module load gcc/4.9.2
module load boost/1.52.0
module load mothur/1.38.1
mothur eDNA.Bacteria.part2.Batch
