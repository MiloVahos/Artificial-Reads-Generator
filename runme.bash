#!/bin/bash

for C in 10 20 30 40 50 60 70 80 90 100; do
    mkdir "Prueba con C = $C"
    ./ARF -DATA lambda_virus.fa -I -Q -L 1024 -B 20 -C $C -P0 0.02
    mv "lambda_virus.fastq" "Prueba con C = $C" 
    mv "lambda_virus.fastqseq" "Prueba con C = $C" 
    mv "lambda_virus.meta" "Prueba con C = $C" 
    mv "lambda_virus.align" "Prueba con C = $C" 
done