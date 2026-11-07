#!/usr/bin/env bash
# Yeast RNA-seq: wild type vs snf2 knock-out (ENA PRJEB5348), 3 vs 3 replicates
set -euo pipefail
mkdir -p ~/rnaseq/{ref,raw,qc,aln,counts} && cd ~/rnaseq

# 1. Run metadata; map each run to its library and lane through the submitted file name
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB5348&result=read_run&fields=run_accession,sample_alias,experiment_alias,library_name,read_count,fastq_bytes&format=tsv" -o runs2.tsv
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB5348&result=read_run&fields=run_accession,run_alias,submitted_ftp,library_name&format=tsv" -o runs4.tsv
awk -F'\t' 'NR>1 {n=split($3,a,"/"); split(a[n],p,"_"); lib=p[5]; sub(/\.sanfastq\.gz/,"",lib); print $1"\t"$4"\t"lib"\t"p[4]}' runs4.tsv > runs_clean.tsv

# 2. Select 6 runs: lane 1 only, >= 1M reads, first 3 libraries per condition
awk -F'\t' 'NR==FNR {rc[$1]=$5; next} $4==1 && rc[$1]>=1000000 {print $1"\t"$2"\t"$3"\t"rc[$1]}' \
  <(tail -n +2 runs2.tsv) runs_clean.tsv | sort -k2,2 -k3,3 > lane1_candidates.tsv
(grep -P "\twild type\t" lane1_candidates.tsv | head -n 3
 grep -P "\tmutant\t" lane1_candidates.tsv | head -n 3) > samples.tsv

# 3. Download (single-end, 51 bp) and QC
(cd raw && for r in $(cut -f1 ../samples.tsv); do prefetch $r && fasterq-dump $r --split-files --threads 2; done)
fastqc raw/*.fastq -o qc -t 2
multiqc qc -o qc

# 4. Reference (yeast R64) and HISAT2 index
B=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64
(cd ref && wget $B/GCF_000146045.2_R64_genomic.fna.gz $B/GCF_000146045.2_R64_genomic.gtf.gz && gunzip *.gz)
hisat2-build -p 2 ref/GCF_000146045.2_R64_genomic.fna ref/sc_r64

# 5. Align and sort
for r in $(cut -f1 samples.tsv); do
  hisat2 -p 2 -x ref/sc_r64 -U raw/$r.fastq 2> aln/$r.hisat2.log | samtools sort -@ 2 -o aln/$r.bam -
  samtools index aln/$r.bam
done

# 6. Count (unstranded: -s 0 chosen after testing -s 0/1/2 on one sample)
featureCounts -s 0 -t exon -g gene_id -a ref/GCF_000146045.2_R64_genomic.gtf \
  -o counts/counts_raw.txt $(cut -f1 samples.tsv | sed 's|^|aln/|; s|$|.bam|')
cut -f1,7- counts/counts_raw.txt | tail -n +2 | sed '1s|aln/||g; 1s|\.bam||g' > counts/counts_matrix.tsv
awk -F'\t' 'BEGIN{OFS="\t"; print "run","condition","library"} {gsub(" ","_",$2); print $1,$2,$3}' samples.tsv > counts/sample_info.tsv
awk -F'\t' '$3=="gene" {id=""; nm=""; if (match($9,/gene_id "[^"]+"/)) id=substr($9,RSTART+9,RLENGTH-10); if (match($9,/; gene "[^"]+"/)) nm=substr($9,RSTART+8,RLENGTH-9); if (nm=="") nm=id; print id"\t"nm}' ref/*.gtf | sort -u > counts/gene_names.tsv
