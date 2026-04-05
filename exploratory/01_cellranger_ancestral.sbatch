#!/bin/bash
#SBATCH --job-name=cellranger_ancestral
#SBATCH --output=cellranger_ancestral.out
#SBATCH --error=cellranger_ancestral.err
#SBATCH --account=pi-araman
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2000

date
export PATH=/project2/araman/annisa/cellranger-7.1.0:$PATH
which cellranger

# YPD-based conditions
cellranger count --id=ypd --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd --sample=DP-AD-1s-anc-ypd --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-37 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-37 --sample=DP-AD-1s-anc-ypd37 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-diam --sample=DP-AD-1s-anc-ypddiam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-h2o2 --sample=DP-AD-1s-anc-ypdh2o2 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-kcl --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-kcl --sample=DP-AD-1s-anc-ypdkcl --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-sorb --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-sorb --sample=DP-AD-1s-anc-ypdsorb --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-kcl-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-kcl-diam --sample=DP-AD-1s-anc-ypdkcldiam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-kcl-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypd-kcl-h2o2 --sample=DP-AD-1s-anc-ypdkclh2o2 --expect-cells=2000 --localcores=8 --localmem=64

# SDC-based conditions
cellranger count --id=sdc --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc --sample=DP-AD-1s-anc-sdc --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-37 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-37 --sample=DP-AD-1s-anc-sdc37 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-diam --sample=DP-AD-1s-anc-sdcdiam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-h2o2 --sample=DP-AD-1s-anc-sdch2o2 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-kcl --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-kcl --sample=DP-AD-1s-anc-sdckcl --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-sorb --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-sorb --sample=DP-AD-1s-anc-sdcsorb --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-kcl-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-kcl-diam --sample=DP-AD-1s-anc-sdckcldiam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-kcl-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sdc-kcl-h2o2 --sample=DP-AD-1s-anc-sdckclh2o2 --expect-cells=2000 --localcores=8 --localmem=64

# Ethanol/glycerol conditions (ancestral-only)
cellranger count --id=ypeg-30 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypeg-30 --sample=DP-AD-1s-anc-ypeg30 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypeg-37 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/ypeg-37 --sample=DP-AD-1s-anc-ypeg37 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sceg-30 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sceg-30 --sample=DP-AD-1s-anc-sceg30 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sceg-37 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/full_scales_novaseq/fastq/sceg-37 --sample=DP-AD-1s-anc-sceg37 --expect-cells=2000 --localcores=8 --localmem=64

date
