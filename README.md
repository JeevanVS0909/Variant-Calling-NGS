Variant-Calling-NGS 🧬

Clinical Somatic Whole Exome Sequencing (WES) Pipeline
Automated end-to-end pipeline for tumor-normal paired WES analysis, including QC, alignment, preprocessing, somatic variant calling, CNV detection, and clinical annotation.

Pipeline Overview

This Python-based pipeline is designed to analyze paired tumor-normal WES samples and produce clinically-relevant variant calls.

Key steps include:

Reference Preparation
Checks and indexes the human reference genome (GRCh38)
Creates sequence dictionary and FASTA index for GATK & BWA
Tools: bwa, samtools, gatk
FASTQ Quality Control
Performs adapter trimming, quality filtering, and generates HTML QC report
Tool: fastp
Alignment
Aligns reads to reference genome using BWA-MEM
Sorts, indexes, and validates BAM files
Tools: bwa, samtools
GATK Preprocessing
MarkDuplicates, Base Quality Score Recalibration (BQSR)
Ensures known-sites indexing for dbSNP and Mills indels
Tool: gatk
Somatic Variant Calling
Runs GATK Mutect2 per chromosome
Generates filtered PASS variants
Tools: gatk
Variant Annotation
Annotates somatic variants using VEP with ClinVar, dbSNP, and other clinical annotations
Produces final annotated VCF
Tool: vep
VCF → MAF Conversion
Converts annotated VCF to MAF format for downstream analysis
Tool: vcf2maf
Copy Number Variation (CNV) Analysis
Runs CNVkit on tumor-normal pair
Produces gene-level CNVs with clinical prioritization
Tool: cnvkit
Clinical Annotation (CIViC)
Maps somatic variants to known cancer biomarkers
Provides clinical-grade variant prioritization
Directory Structure
Variant-Calling-NGS/
│
├─ data/                    # Raw FASTQ files (tumor/normal)
├─ resources/               # Reference genome, known sites, CNV targets, VEP cache
├─ results/                 # Output directory for pipeline results
│   ├─ cnvkit/              # CNV analysis outputs
│   ├─ annotated.vcf        # VEP-annotated VCF
│   ├─ somatic.maf          # Final MAF file
│   └─ civic.merged.maf     # CIViC annotated variants
├─ pipeline.py              # Main Python WES pipeline
└─ README.md
Requirements

System Requirements:

Linux / macOS (tested on Ubuntu 22.04)
≥16 CPU threads recommended
≥64GB RAM for Mutect2 + VEP

Software Dependencies:

Python 3.10+
BWA
Samtools
GATK
Fastp
VEP
CNVkit
vcf2maf
R (optional for downstream analysis)
Usage
1. Clone the repository
git clone https://github.com/JeevanVS0909/Variant-Calling-NGS.git
cd Variant-Calling-NGS
2. Prepare resources

Place reference genome, known-sites, VEP cache, and CNV target files in resources/.

3. Add raw data

Place paired tumor-normal FASTQ files in data/ directory.

4. Run the pipeline
python pipeline.py
The pipeline automatically performs all steps from QC → alignment → preprocessing → somatic calling → CNV → annotation.
Logs are stored in results/log.txt.
Outputs
File	Description
results/*.bam	Aligned and preprocessed BAM files
results/somatic_PASS.vcf.gz	High-confidence somatic variants
results/annotated.vcf	VEP annotated variants
results/somatic.maf	MAF format for analysis
results/cnv_gene_level.csv	Gene-level CNVs
results/civic.merged.maf	Clinical-grade CIViC annotation
results/log.txt	Pipeline execution log
Example
python pipeline.py

Output example:

✔ Skipping FASTP for TUMOR
✔ Skipping alignment for TUMOR
✔ MarkDuplicates completed
✔ BaseRecalibrator completed
✔ Mutect2 Somatic Calling Completed Successfully
✔ VEP Annotation Completed
✔ CNV Gene-level Annotation Completed
✔ CIViC Annotation Completed
References
Van der Auwera, G.A., & O’Connor, B.D. (2020). Genomics in the Cloud. O’Reilly.
Cibulskis, K. et al., (2013). Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnol.
McLaren, W. et al., (2016). The Ensembl Variant Effect Predictor. Genome Biol.
Talevich, E. et al., (2016). CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol.
