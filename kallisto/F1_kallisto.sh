#!/bin/bash

### Retrieve the raw data from cabernet to whitney
scp lynnly@cabernet.genomecenter.ucdavis.edu:/share/malooflab/Fastqs/KIAT/F1/415F1-young* .

for f1 in *1.fq.gz
do
  base=${f1%%1.fq.gz}
  
  ### Use trimmomatic to do quality control and trimming on reads
  java -jar \
  ~/../ruijuanli/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
  PE ${base}1.fq.gz ${base}2.fq.gz ${base}1.trimmed.fq.gz \
  ${base}1.unpaired.fq.gz ${base}2.trimmed.fq.gz ${base}2.unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  ### Use kallisto to pseudo-align the reads, using 100 bootstraps for future sleuth analysis
  kallisto quant -i \
  ~/../ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5.cds.19.kai \
  -o ${base}.dir --bootstrap-samples=100 --threads=8 ${base}1.trimmed.fq.gz ${base}2.trimmed.fq.gz
        
done