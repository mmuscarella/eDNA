#!/bin/bash

wget http://www.mothur.org/w/images/2/27/Silva.nr_v119.tgz
tar -xzvf Silva.nr_v119.tgz

wget http://mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
unzip -j Silva.gold.bacteria.zip

wget http://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip
unzip -j Trainset9_032012.pds.zip

module load gcc/4.9.2
module load mothur/1.36.1
mothur "#pcr.seqs(fasta=silva.nr_v119.align, start=11894, end=25319, keepdots=F, processors=1)"
mv silva.nr_v119.pcr.align silva.v4.fasta
