#!/bin/bash

# read simulator
mkdir tools; cd tools
git clone https://github.com/lh3/wgsim.git 2>/dev/null
cd wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
cd ../..

# sample data
mkdir genome
cd genome
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/genome2/CultivarB.fa 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/genome2/CultivarB.fa.fai 2>/dev/null
cd ../
mkdir reads
cd reads
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads2/high_bulked_read1.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads2/high_bulked_read2.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads2/low_bulked_read1.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads2/low_bulked_read2.fastq 2>/dev/null
# cd ../
# mkdir simulation
# cd simulation
# wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/simulation/mutations.fa 2>/dev/null
