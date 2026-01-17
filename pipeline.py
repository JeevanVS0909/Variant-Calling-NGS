import os, subprocess, sys, shutil, csv, gzip
import pandas as pd

# === CONFIGURATION ===
FASTQ_FILE_1 = "data/SRR8261574_1.fastq"
FASTQ_FILE_2 = "data/SRR8261574_2.fastq"
REFERENCE_GENOME = "GCF_000001405.40_GRCh38.p14_genomic.fna"
MODEL_TYPE = "WGS"

SAM_FILE = "aligned_reads.sam"
BAM_FILE = "aligned_reads.bam"
SORTED_BAM_FILE = "sorted_reads.bam"
VCF_OUTPUT = "deepvariant_output.vcf.gz"
GVCF_OUTPUT = "deepvariant_output.g.vcf.gz"
DEEPVARIANT_CMD = "/opt/deepvariant/bin/run_deepvariant"

OUTPUT_DIR = "deepvariant_results"
FASTQC_DIR = os.path.join(OUTPUT_DIR, "fastqc")
VCF_DIR = "VCF"
CSV_OUTPUT = os.path.join(VCF_DIR, "Biomarkers.csv")
PARSED_OUTPUT = os.path.join(VCF_DIR, "parsed_ann_output.csv")
FILTERED_OUTPUT = os.path.join(VCF_DIR, "filtered_biomarkers.csv")

SNPEFF_JAR = "snpEff/snpEff.jar"
SNPEFF_DB = "snpEff/data/GRCh38.105"

# === SETUP ===
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(VCF_DIR, exist_ok=True)

# === UTILITY FUNCTIONS ===
def check_tool(tool_name):
    if shutil.which(tool_name) is None:
        print(f"❌ '{tool_name}' not in PATH.")
        sys.exit(1)

def check_required_tools():
    for tool in ["bwa", "samtools", "java"]:
        check_tool(tool)

def run_command(cmd, description):
    print(f"\n🔹 {description}...")
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("✅ Done.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

# === PIPELINE STEPS ===

def run_fastqc():
    tag = os.path.join(OUTPUT_DIR, "fastqc")
    ensure_dir(tag)
    fq1_html = os.path.join(tag, os.path.basename(FASTQ_FILE_1).replace(".fastq", "_fastqc.html"))
    fq2_html = os.path.join(tag, os.path.basename(FASTQ_FILE_2).replace(".fastq", "_fastqc.html"))
    if os.path.exists(fq1_html) and os.path.exists(fq2_html):
        print("🔹 Skipping FastQC (already done)")
        return
    run_command(f"fastqc {FASTQ_FILE_1} {FASTQ_FILE_2} -o {tag}", "Running FastQC")

def index_reference():
    required_exts = ["bwt", "pac", "ann", "amb", "sa"]
    if not all(os.path.exists(f"{REFERENCE_GENOME}.{ext}") for ext in required_exts):
        run_command(f"bwa index {REFERENCE_GENOME}", "Indexing reference genome")

def index_reference_fasta():
    fai_file = REFERENCE_GENOME + ".fai"
    if not os.path.exists(fai_file):
        run_command(f"samtools faidx {REFERENCE_GENOME}", "Indexing reference FASTA (.fai)")

def align_reads():
    if not os.path.exists(SAM_FILE):
        run_command(f"bwa mem {REFERENCE_GENOME} {FASTQ_FILE_1} {FASTQ_FILE_2} > {SAM_FILE}", "Aligning reads")

def convert_sam_to_bam():
    if not os.path.exists(BAM_FILE):
        run_command(f"samtools view -S -b {SAM_FILE} > {BAM_FILE}", "Converting SAM to BAM")

def sort_bam():
    if not os.path.exists(SORTED_BAM_FILE):
        run_command(f"samtools sort {BAM_FILE} -o {SORTED_BAM_FILE}", "Sorting BAM")

def index_bam():
    bai = SORTED_BAM_FILE + ".bai"
    if not os.path.exists(bai):
        run_command(f"samtools index {SORTED_BAM_FILE}", "Indexing BAM")

import os
import subprocess
import sys

# === CONFIGURATION ===
REFERENCE_GENOME = "GCF_000001405.40_GRCh38.p14_genomic.fna"
SORTED_BAM_FILE = "sorted_reads.bam"
VCF_OUTPUT = "output.vcf.gz"
GVCF_OUTPUT = "output.g.vcf.gz"

def run_deepvariant():
    print("\n🔹 Checking DeepVariant output...")

    output_dir = os.path.join(os.getcwd(), "deepvariant_results")
    os.makedirs(output_dir, exist_ok=True)

    vcf_path = os.path.join(output_dir, VCF_OUTPUT)
    gvcf_path = os.path.join(output_dir, GVCF_OUTPUT)

    if os.path.exists(vcf_path) and os.path.exists(gvcf_path):
        print(f"⏭️  DeepVariant skipped (outputs already exist):\n  - {vcf_path}\n  - {gvcf_path}")
        return

    print("🚀 Running DeepVariant (native installation)...")

    deepvariant_cmd = [
        "run_deepvariant",  
        "--model_type=WGS",
        f"--ref={REFERENCE_GENOME}",
        f"--reads={SORTED_BAM_FILE}",
        f"--output_vcf=deepvariant_results/{VCF_OUTPUT}",
        f"--output_gvcf=deepvariant_results/{GVCF_OUTPUT}",
        "--num_shards=8"
    ]

    try:
        subprocess.run(deepvariant_cmd, check=True)
        print(f"✅ DeepVariant completed:\n  - VCF: deepvariant_results/{VCF_OUTPUT}\n  - gVCF: deepvariant_results/{GVCF_OUTPUT}")
    except subprocess.CalledProcessError as e:
        print(f"❌ DeepVariant failed with error:\n{e}")
        sys.exit(1)

def unzip_vcf(gz_path):
    vcf_path = gz_path.replace(".gz", "")
    if not os.path.exists(vcf_path):
        with gzip.open(gz_path, 'rt') as f_in, open(vcf_path, 'w') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return vcf_path

def rename_chromosomes(vcf):
    renamed = vcf.replace(".vcf", "_renamed.vcf")
    if os.path.exists(renamed):
        return renamed
    rename_map = {
        "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3",
        "NC_000004.12": "4", "NC_000005.10": "5", "NC_000006.12": "6",
        "NC_000007.14": "7", "NC_000008.11": "8", "NC_000009.12": "9",
        "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
        "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15",
        "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18",
        "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9": "21",
        "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y",
        "NC_012920.1": "MT"
    }
    with open(vcf, "r") as fin, open(renamed, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            else:
                parts = line.strip().split("\t")
                parts[0] = rename_map.get(parts[0], parts[0])
                fout.write("\t".join(parts) + "\n")
    return renamed

def annotate_vcf_with_snpeff():
    print("\n🔹 Annotating VCF with SnpEff...")

    snpeff_dir = "snpEff"
    snpeff_jar = os.path.join(snpeff_dir, "snpEff.jar")
    genome_version = "GRCh38.105"

    input_vcf = "deepvariant_results/output_renamed.vcf"
    annotated_vcf = "deepvariant_results/output_renamed_snpeff.vcf"

    if not os.path.exists(input_vcf):
        print(f"❌ Input VCF not found: {input_vcf}")
        sys.exit(1)

    if os.path.exists(annotated_vcf):
        print(f"⚠️ Skipping annotation: {annotated_vcf} already exists.")
        return annotated_vcf

    cmd = [
        "java", "-Xmx4g", "-jar", snpeff_jar,
        genome_version,
        input_vcf
    ]

    try:
        print(f"🛠️ Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        with open(annotated_vcf, "w") as out_f:
            out_f.write(result.stdout)
        print(f"✅ Annotation complete. Output saved to: {annotated_vcf}")
    except subprocess.CalledProcessError as e:
        print("❌ SnpEff annotation failed with the following output:\n")
        print(e.output)
        sys.exit(1)

    return annotated_vcf

def vcf_to_csv_large(input_vcf_path, output_csv_path):
    encodings = ['utf-8', 'utf-8-sig', 'utf-16', 'latin-1']
    info_fields, header_fields = set(), []
    for enc in encodings:
        try:
            with open(input_vcf_path, 'r', encoding=enc) as vcf:
                for line in vcf:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#CHROM'):
                        header_fields = line.lstrip('#').strip().split('\t')[:12]
                        break
                for line in vcf:
                    if not line.startswith('#'):
                        info = line.strip().split('\t')[7]
                        for field in info.split(';'):
                            key = field.split('=')[0] if '=' in field else field
                            info_fields.add(key)
            break
        except UnicodeDecodeError:
            continue

    fieldnames = header_fields + sorted(info_fields)
    with open(input_vcf_path, 'r', encoding='utf-8') as vcf, open(output_csv_path, 'w', newline='', encoding='utf-8') as out_csv:
        writer = csv.DictWriter(out_csv, fieldnames=fieldnames)
        writer.writeheader()
        for line in vcf:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            info_dict = {}
            for item in cols[7].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True
            base = dict(zip(header_fields, cols[:8]))
            writer.writerow({**base, **info_dict})
    print(f"✅ CSV written: {output_csv_path}")

def parse_ann_to_biomarkers_streaming(csv_input, parsed_output, filtered_output):
    import csv
    ann_headers = [
        "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type",
        "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length",
        "CDS.pos / CDS.length", "AA.pos / AA.length", "Distance", "Errors_Warnings_Info", "Genotype_Number"
    ]
    biomarker_genes = [
        "PIK3CA","AKT1","PTEN","ESR1","HER2","BRCA1","BRCA2","FGFR1","FGFR2","FGFR3"
    ]
    # Set up fieldnames/order for outgoing CSVs
    out_headers = ann_headers + ['CHROM', 'POS', 'REF', 'ALT']
    with open(csv_input, encoding='utf-8') as infile, \
         open(parsed_output, 'w', newline='', encoding='utf-8') as out_parsed, \
         open(filtered_output, 'w', newline='', encoding='utf-8') as out_biom:
        reader = csv.DictReader(infile)
        parsed_writer = csv.DictWriter(out_parsed, fieldnames=out_headers)
        biom_writer = csv.DictWriter(out_biom, fieldnames=out_headers)
        parsed_writer.writeheader()
        biom_writer.writeheader()
        for row in reader:
            ann_data = row.get('ANN', '')
            if not ann_data:
                continue
            for entry in ann_data.split(','):
                fields = entry.split('|') + [''] * (len(ann_headers) - len(entry.split('|')))
                record = dict(zip(ann_headers, fields))
                record.update({k: row.get(k, '') for k in ['CHROM', 'POS', 'REF', 'ALT']})
                parsed_writer.writerow(record)
                if record["Gene_Name"] in biomarker_genes:
                    biom_writer.writerow(record)
    print(f"Parsed annotations saved: {parsed_output}")
    print(f"Biomarker variants saved: {filtered_output}")

# === MAIN PIPELINE ===
def main():
    run_fastqc()
    check_required_tools()
    index_reference()
    index_reference_fasta()
    align_reads()
    convert_sam_to_bam()
    sort_bam()
    index_bam()
    run_deepvariant()

    final_vcf = unzip_vcf(os.path.join(OUTPUT_DIR, VCF_OUTPUT))
    renamed_vcf = rename_chromosomes(final_vcf)
    annotated_vcf = annotate_vcf_with_snpeff()
    os.makedirs(os.path.dirname(CSV_OUTPUT), exist_ok=True)
    vcf_to_csv_large(annotated_vcf, CSV_OUTPUT)
    parse_ann_to_biomarkers_streaming(CSV_OUTPUT, PARSED_OUTPUT, FILTERED_OUTPUT)

if __name__ == "__main__":
    main()
