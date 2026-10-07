# NGS pipelines for Illumina data: full commands

Companion to the paper *What can the university do with its Illumina sequencer: analyses, Galaxy
pipelines and conditions for their use* (Annual of Assen Zlatarov University, Burgas, 2026).
Table 3 of the paper gives each program with its repository and a short call; this page gives the
complete command line for paired-end Illumina data.

In Galaxy every program below is a tool form whose fields match these options; the command lines are
for running the same step on the server's shell or for checking what a Galaxy tool does.

## Conventions

| Name | Meaning |
|---|---|
| `R1.fq.gz`, `R2.fq.gz` | paired-end reads (gzip) |
| `ref.fa` | reference genome (FASTA), indexed where the tool needs it |
| `s.bam`, `s.vcf.gz` | sorted alignments / variants of one sample |
| `-t 8`, `-p 8`, `--cpus 8` | thread count; in Galaxy this is `${GALAXY_SLOTS}` |

**Verification.** Options marked ✓ were checked on the university Galaxy host (ARM Neoverse-N1,
aarch64, Ubuntu 25.10) on 7 October 2026 against the `--help` output of the Bioconda version given;
the script is [`check_tools.sh`](check_tools.sh). Options marked ◦ were checked against the tool's
own README. The calls show the options; they were not run on real data.

**ARM note.** The server is aarch64. Install tools from Bioconda (it has aarch64 builds); most
BioContainers images are x86-only. NVIDIA Parabricks needs an NVIDIA Grace CPU on ARM and stops
with `SIGILL` on Neoverse-N1.

## 1. Trimming and quality control

```bash
# fastp: adapters, quality, poly-G, filtering and an HTML/JSON report in one pass ◦
fastp -i R1.fq.gz -I R2.fq.gz -o t1.fq.gz -O t2.fq.gz \
      --detect_adapter_for_pe -w 8 -j fastp.json -h fastp.html

# Trimmomatic 0.41 ✓
trimmomatic PE -threads 8 R1.fq.gz R2.fq.gz p1.fq.gz u1.fq.gz p2.fq.gz u2.fq.gz \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

# Trim Galore ◦  and cutadapt ◦ (cutadapt also removes 16S primers)
trim_galore --paired --fastqc R1.fq.gz R2.fq.gz
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -m 30 -o t1.fq.gz -p t2.fq.gz R1.fq.gz R2.fq.gz

# FastQC, Falco ◦ (same report, faster) and MultiQC ◦ (one report for all samples)
fastqc -t 4 R1.fq.gz R2.fq.gz
falco -o falco_out R1.fq.gz
multiqc .
```

## 2. Alignment

```bash
# BWA-MEM2 (BWA has the same syntax: bwa index / bwa mem)
bwa-mem2 index ref.fa
bwa-mem2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -@ 8 -o s.bam -

# Bowtie2 ◦
bowtie2-build ref.fa idx
bowtie2 -p 8 -x idx -1 R1.fq.gz -2 R2.fq.gz -S s.sam

# HISAT2 2.2.3 ✓ (RNA-seq; --dta for transcript assembly)
hisat2 -p 8 -x idx -1 R1.fq.gz -2 R2.fq.gz -S s.sam --dta

# STAR (RNA STAR in Galaxy; about 30 GB RAM for a human index)
STAR --runThreadN 16 --genomeDir star_idx --readFilesIn R1.fq.gz R2.fq.gz \
     --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

# Salmon ◦ (transcript quantification without alignment)
salmon quant -i salmon_idx -l A -1 R1.fq.gz -2 R2.fq.gz -p 8 -o salmon_out

# Snippy ◦ (reads to SNPs for haploid genomes)
snippy --cpus 8 --ref ref.gbk --R1 R1.fq.gz --R2 R2.fq.gz --outdir snippy_out
```

## 3. De novo assembly

```bash
# Shovill ◦ (SPAdes inside) and QUAST
shovill --R1 R1.fq.gz --R2 R2.fq.gz --outdir shovill_out --cpus 8
quast.py shovill_out/contigs.fa -r ref.fa -o quast_out
# Flye ◦ is for long reads only (Nanopore, PacBio), not for Illumina data:
# flye --nano-hq reads.fq.gz -o flye_out -t 16
```

## 4. Sorting and duplicates

```bash
samtools sort -@ 8 -o s.bam in.bam
samtools index s.bam
# Picard 3.5.0 ✓ (new syntax; the old I=/O= form still runs but is deprecated)
java -jar picard.jar SortSam -I in.sam -O s.bam -SO coordinate
java -jar picard.jar MarkDuplicates -I s.bam -O d.bam -M dup_metrics.txt
```

## 5. Variant calling

```bash
# GATK 4.6.2.0 ✓
gatk HaplotypeCaller -R ref.fa -I s.bam -O s.vcf.gz            # add -ERC GVCF for cohorts
gatk Mutect2 -R ref.fa -I tumor.bam -I normal.bam -normal NORMAL_SAMPLE -O somatic.vcf.gz
gatk FilterMutectCalls -R ref.fa -V somatic.vcf.gz -O somatic.filtered.vcf.gz

# bcftools 1.24 ✓
bcftools mpileup -f ref.fa s.bam | bcftools call -mv -Oz -o s.vcf.gz
```

## 6. Annotation, filtering and viewing

```bash
# SnpEff 5.4c ✓ (DB = genome version, e.g. from: snpEff databases)
snpEff ann DB s.vcf.gz > s.ann.vcf

# Funcotator (GATK 4.6.2.0) ✓
gatk Funcotator --variant s.vcf.gz --reference ref.fa --ref-version hg38 \
     --data-sources-path funcotator_dataSources -O s.func.vcf --output-file-format VCF

# ANNOVAR (free for academic use after registration at annovar.openbioinformatics.org; not verified)
table_annovar.pl s.vcf humandb/ -buildver hg38 -out ann -protocol refGene -operation g -vcfinput

# VCFtools 0.1.17 ✓ and bcftools filter ✓
vcftools --gzvcf s.vcf.gz --minQ 30 --max-missing 0.9 --recode --out filtered
bcftools view -i 'QUAL>30 && DP>10' s.vcf.gz -Oz -o filtered.vcf.gz
```

JBrowse 2 (GMOD/jbrowse-components): load the reference, the BAM and the VCF as tracks; in Galaxy use
the JBrowse visualisation of a history dataset.

## 7. Expression

```bash
# featureCounts, Subread 2.1.1 ✓ (--countReadPairs counts fragments, not reads)
featureCounts -T 8 -p --countReadPairs -a genes.gtf -o counts.txt s1.bam s2.bam s3.bam
```

```r
# DESeq2 1.50.2 ✓
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# edgeR (Bioconductor)
library(edgeR)
y <- DGEList(counts = counts, group = coldata$condition)
y <- normLibSizes(y)        # calcNormFactors() in edgeR < 4
```

## 8. Pathways and gene ontology

```r
# clusterProfiler 4.18.4 ✓ (YuLab-SMU/clusterProfiler)
library(clusterProfiler); library(org.Hs.eg.db)
ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
ekg <- enrichKEGG(gene = genes, organism = "hsa")
```

EnrichR: paste the gene list at maayanlab.cloud/Enrichr. IPA (QIAGEN) is commercial and needs a licence.

## 9. Association, miRNA and methylation

```bash
# PLINK 2.0 ✓ (needs hundreds of samples)
plink2 --vcf s.vcf.gz --make-bed --out x
plink2 --bfile x --pheno pheno.txt --glm --out assoc

# miRDeep2 2.0.1.3 ✓ (small-RNA libraries; miRBase files from mirbase.org)
mapper.pl reads.fa -c -j -m -p genome_idx -s reads_collapsed.fa -t reads_vs_genome.arf
quantifier.pl -p precursors.fa -m mature.fa -r reads_collapsed.fa -t hsa

# Bismark 3.1.0 ✓ (bisulfite libraries)
bismark --genome genome_folder -1 R1.fq.gz -2 R2.fq.gz
bismark_methylation_extractor --gzip --bedGraph R1_bismark_bt2_pe.bam
```

## Repositories

| Program | Repository |
|---|---|
| fastp | https://github.com/OpenGene/fastp |
| Trimmomatic | https://github.com/usadellab/Trimmomatic |
| Trim Galore | https://github.com/FelixKrueger/TrimGalore |
| cutadapt | https://github.com/marcelm/cutadapt |
| FastQC | https://github.com/s-andrews/FastQC |
| Falco | https://github.com/smithlabcode/falco |
| MultiQC | https://github.com/MultiQC/MultiQC |
| BWA / BWA-MEM2 | https://github.com/lh3/bwa, https://github.com/bwa-mem2/bwa-mem2 |
| Bowtie2 | https://github.com/BenLangmead/bowtie2 |
| HISAT2 | https://github.com/DaehwanKimLab/hisat2 |
| STAR | https://github.com/alexdobin/STAR |
| Salmon | https://github.com/COMBINE-lab/salmon |
| Snippy | https://github.com/tseemann/snippy |
| Shovill / SPAdes | https://github.com/tseemann/shovill, https://github.com/ablab/spades |
| QUAST | https://github.com/ablab/quast |
| Flye | https://github.com/mikolmogorov/Flye |
| minimap2 | https://github.com/lh3/minimap2 |
| samtools / bcftools | https://github.com/samtools/samtools, https://github.com/samtools/bcftools |
| Picard | https://github.com/broadinstitute/picard |
| GATK (HaplotypeCaller, Mutect2, Funcotator) | https://github.com/broadinstitute/gatk |
| SnpEff | https://github.com/pcingola/SnpEff |
| ANNOVAR | https://annovar.openbioinformatics.org (registration) |
| VCFtools | https://github.com/vcftools/vcftools |
| JBrowse 2 | https://github.com/GMOD/jbrowse-components |
| featureCounts (Subread) | https://github.com/ShiLab-Bioinformatics/subread |
| DESeq2 | https://github.com/thelovelab/DESeq2 |
| edgeR | https://bioconductor.org/packages/edgeR |
| clusterProfiler | https://github.com/YuLab-SMU/clusterProfiler |
| EnrichR | https://maayanlab.cloud/Enrichr |
| PLINK 2 | https://github.com/chrchang/plink-ng |
| miRDeep2 / miRBase | https://github.com/rajewsky-lab/mirdeep2, https://mirbase.org |
| Bismark | https://github.com/FelixKrueger/Bismark |
