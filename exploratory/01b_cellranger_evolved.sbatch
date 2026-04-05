#!/bin/bash
#SBATCH --job-name=cellranger_evolved
#SBATCH --output=cellranger_evolved.out
#SBATCH --error=cellranger_evolved.err
#SBATCH --account=pi-araman
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2000

date
export PATH=/project2/araman/annisa/cellranger-7.1.0:$PATH
which cellranger

# YPD-based conditions
cellranger count --id=ypd --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-37 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-37 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-diam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-h2o02 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-kcl --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-kcl --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-sorb --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-sorb --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-kcl-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-kcl-diam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypd-kcl-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypd-kcl-h2o02 --expect-cells=2000 --localcores=8 --localmem=64

# SDC-based conditions
cellranger count --id=sdc --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-37 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-37 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-diam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-h2o02 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-kcl --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-kcl --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-sorb --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-sorb --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-kcl-diam --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-kcl-diam --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=sdc-kcl-h2o2 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sdc-kcl-h2o02 --expect-cells=2000 --localcores=8 --localmem=64

# Ethanol/glycerol condition
cellranger count --id=sceg-30 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-sceg-30 --expect-cells=2000 --localcores=8 --localmem=64
cellranger count --id=ypeg-30 --transcriptome=/project2/araman/annisa/Saccharomyces.cerevisiae_genome --fastqs=/project2/araman/annisa/evolved_scales_novaseq/fastq --sample=DP-AD-1s-evo-ypeg-30 --expect-cells=2000 --localcores=8 --localmem=64

date
