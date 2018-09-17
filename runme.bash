#!/bin/bash

for C in 1 10 20 30 40 50 60 70 80 90 100; do
    mkdir "Prueba con C = $C"
    ./ARF -DATA GRCh38.fa -I -Q -L 1024 -B 20 -C $C -P0 0.3625
    mv "GRCh38.fastq" "Prueba con C = $C" 
    mv "GRCh38.fastqseq" "Prueba con C = $C" 
    mv "GRCh38.meta" "Prueba con C = $C" 
    mv "GRCh38.align" "Prueba con C = $C" 
done