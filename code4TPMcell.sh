#!/bin/bash
## 质量控制
conda activate kneaddata
mkdir -p temp/qc
sed -n '1,40p' lst.txt|cut -f1|./rush -j 2 \
      "fastp -i seq/{1}_1.fastq -I seq/{1}_2.fastq \
        -o temp/qc/{1}_1.fastq  -O temp/qc/{1}_2.fastq"

mkdir -p temp/fastqc
fastqc temp/qc/*.fastq -t 2 -o temp/fastqc/ > temp/temp.log
    multiqc -d temp/fastqc/ -o result/qc/
## 物种识别（默认使用metaphlan）
# # kraken2
# conda activate kraken2
# mkdir -p temp/kraken2
# db=/media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/db
# sed -n '201,217p' litelst.txt|cut -f1|./rush -j 2 \
#      "kraken2 --db ${db}/kraken2/pluspf --paired temp/qc/{1}_?.fastq \
#      --threads 10 --use-names --report-zero-counts \
#      --report temp/kraken2/{1}.report \
#      --output temp/kraken2/{1}.output"
# metaphlan
conda activate humann3
mkdir -p temp/metaphlan/
sed -n '1,40p' lst.txt|cut -f1|./rush -j 6 \
     'metaphlan temp/qc/{1}_1.fastq,temp/qc/{1}_2.fastq --input_type fastq --bowtie2out temp/metaphlan/{1}.bowtie2out.txt -o temp/metaphlan/{1}_tavg_g_profiled_metagenome.txt --nproc 10'
## TPM/cell数据处理（need run.sh）
conda activate qsc
mkdir -p temp/concat
for i in `sed -n '1,40p' lst.txt|cut -f1`;do
     echo ${i};
     # 双端合并 
     cat temp/qc/${i}_?.fastq > temp/concat/${i}.fq
     echo "cat done"
     # TPM/cell（默认采用kraken+metaphlan的物种结果）
     time bash run.sh ${i} globalsoil /media/dell/55b2b60a-b2b1-4c7d-a4d3-05b1959cfe72/data/QSC;
     rm temp/concat/${i}.fq
     gzip temp/qc/${i}*.fastq
done
