#!/bin/bash

# read simulator
mkdir tools; cd tools
git clone https://github.com/lh3/wgsim.git
cd wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
cd ../..

# sample data
mkdir genome
cd genome
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/genome/CultivarA.fa 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/genome/CultivarA.fa.fai 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/genome/CultivarB.fa 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/genome/CultivarB.fa.fai 2>/dev/null
cd ../
mkdir reads
cd reads
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads/Mutated_Cultivar_read1.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads/Mutated_Cultivar_read2.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads/Mutated_Cultivar2_read1.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads/Mutated_Cultivar2_read2.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads/bulked_read1.fastq 2>/dev/null
wget https://raw.githubusercontent.com/slt666666/FAO_lecture/main/FAO_2024/data/reads/bulked_read2.fastq 2>/dev/null