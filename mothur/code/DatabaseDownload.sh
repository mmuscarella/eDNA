#!/bin/bash

wget https://www.mothur.org/w/images/b/be/Silva.nr_v123.tgz
tar -xzvf Silva.nr_v123.tgz

wget https://www.mothur.org/w/images/f/f1/Silva.gold.bacteria.zip
unzip -j Silva.gold.bacteria.zip

wget https://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
tar -xzvf  Trainset14_032015.pds.tgz

module load gcc/4.9.2
module load boost/1.52.0
module load mothur/1.38.1

mothur "#pcr.seqs(fasta=silva.nr_v123.align, start=11894, end=25319, keepdots=F, processors=1)"
mv silva.nr_v123.pcr.align silva.v4.fasta
