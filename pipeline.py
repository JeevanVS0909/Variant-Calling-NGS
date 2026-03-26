############################################################
# CLINICAL SOMATIC WES PIPELINE 
############################################################

import os
import subprocess
import pandas as pd
import time
import pathlib

############################################################
# GLOBAL CONFIGURATION
############################################################

THREADS = 16
OUTPUT = "results"
os.makedirs(OUTPUT, exist_ok=True)

REFERENCE = "resources/Homo_sapiens_assembly38.fasta"
GERMLINE  = "resources/af-only-gnomad.hg38.vcf.gz"
PON       = "resources/1000g_pon.hg38.vcf.gz"
RESOURCES = "/app/resources"

KNOWN_SITES = [
    "resources/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    "resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
]

TUMOR_R1  = "data/tumor_SRR6656083_1.fastq.gz"
TUMOR_R2  = "data/tumor_SRR6656083_2.fastq.gz"
NORMAL_R1 = "data/normal_SRR5801907_1.fastq.gz"
NORMAL_R2 = "data/normal_SRR5801907_2.fastq.gz"

############################################################
# LOGGING SYSTEM (PERSISTENT SINGLE LOG FILE)
############################################################

LOG_FILE = os.path.join(OUTPUT, "log.txt")

def log(message):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as f:
        f.write(f"[{timestamp}] {message}\n")

# Mark start of a new run
log("\n" + "="*70)
log("🚀 NEW PIPELINE RUN STARTED")
log("="*70)


############################################################
# UTILITIES
############################################################

def run(cmd, step):
    log(f"\n--- STARTING: {step} ---")
    print(f"\n STARTING → {step}")

    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    for line in process.stdout:
        print(line, end="")
        log(line.strip())

    process.wait()

    if process.returncode != 0:
        log(f"❌ FAILED: {step}")
        raise RuntimeError(f"FAILED → {step}")

    log(f"✅ COMPLETED: {step}")
    print(f" COMPLETED → {step}")



def exists(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


############################################################
# 1️⃣ REFERENCE PREPARATION
############################################################

def prepare_reference():
    global REFERENCE

    if not os.path.exists(REFERENCE):
        raise RuntimeError("Reference genome not found.")

    if REFERENCE.endswith(".gz"):
        fasta = REFERENCE[:-3]
        if not os.path.exists(fasta):
            run(f"gunzip -k {REFERENCE}", "UNZIP REFERENCE")
        REFERENCE = fasta

    if not os.path.exists(REFERENCE + ".bwt"):
        run(f"bwa index {REFERENCE}", "BWA INDEX")

    if not os.path.exists(REFERENCE + ".fai"):
        run(f"samtools faidx {REFERENCE}", "FASTA INDEX")

    dict_file = str(pathlib.Path(REFERENCE).with_suffix(".dict"))
    if not os.path.exists(dict_file):
        run(
            f"gatk CreateSequenceDictionary -R {REFERENCE} -O {dict_file}",
            "SEQ DICT"
        )


############################################################
# 2️⃣ FASTP QC
############################################################

def fastp(sample, r1, r2):
    out1 = f"{OUTPUT}/{sample}_R1.fastq.gz"

    if exists(out1):
        print(f" Skipping fastp for {sample}")
        return

    run(f"""
    fastp \
    -i {r1} -I {r2} \
    -o {OUTPUT}/{sample}_R1.fastq.gz \
    -O {OUTPUT}/{sample}_R2.fastq.gz \
    --thread {THREADS} \
    --html {OUTPUT}/{sample}_fastp.html
    """, f"FASTP {sample}")


############################################################
# 3️⃣ ALIGNMENT
############################################################

def align(sample):

    bam = f"{OUTPUT}/{sample}.bam"
    bai = bam + ".bai"
    unsorted = f"{OUTPUT}/{sample}_unsorted.bam"

    r1 = f"{OUTPUT}/{sample}_R1.fastq.gz"
    r2 = f"{OUTPUT}/{sample}_R2.fastq.gz"

    if os.path.exists(bam) and os.path.exists(bai):
        print(f" Skipping alignment for {sample}")
        return

    if not os.path.exists(unsorted):
        run(f"""
        bwa mem -t {THREADS} \
        -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' \
        {REFERENCE} {r1} {r2} | \
        samtools view -Sb - > {unsorted}
        """, f"BWA MEM {sample}")

    if not os.path.exists(bam):
        run(f"""
        samtools sort \
        -@ 4 -m 1G \
        -o {bam} {unsorted}
        """, f"SORT {sample}")

    run(f"samtools quickcheck {bam}", "BAM VALIDATION")

    if not os.path.exists(bai):
        run(f"samtools index {bam}", "BAM INDEX")

    if os.path.exists(unsorted):
        os.remove(unsorted)


############################################################
# 4️⃣ GATK PREPROCESSING 
############################################################

def preprocess(sample):

    final = f"{OUTPUT}/{sample}_final.bam"
    final_bai = final + ".bai"
    dedup = f"{OUTPUT}/{sample}_dedup.bam"
    dedup_bai = dedup + ".bai"
    recal_table = f"{OUTPUT}/{sample}_recal.table"

    print(f"\n========== PREPROCESSING {sample} ==========")

    # -----------------------------------------------------
    # Ensure reference indexes exist (CRITICAL)
    # ------------------------------------------------------
    if not exists(REFERENCE + ".fai"):
        run(f"samtools faidx {REFERENCE}", "Index reference FASTA")

    dict_file = REFERENCE.replace(".fasta", ".dict")
    if not exists(dict_file):
        run(f"gatk CreateSequenceDictionary -R {REFERENCE}", "Create reference dict")

    # ------------------------------------------------------
    # Ensure known-sites are indexed
    # ------------------------------------------------------
    for vcf in KNOWN_SITES:
        if not exists(vcf + ".tbi"):
            run(f"tabix -p vcf {vcf}", f"Index {vcf}")

    # ------------------------------------------------------
    # MarkDuplicates
    # ------------------------------------------------------
    if not exists(dedup):
        run(f"""
        gatk --java-options "-Xmx4G" MarkDuplicates \
        -I {OUTPUT}/{sample}.bam \
        -O {dedup} \
        -M {OUTPUT}/{sample}_metrics.txt \
        --CREATE_INDEX true
        """, "MarkDuplicates")
    else:
        print(" Skipping MarkDuplicates (dedup BAM exists)")

    # Ensure dedup index exists
    if not exists(dedup_bai):
        run(f"samtools index {dedup}", "Index dedup BAM")

    # ------------------------------------------------------
    # BaseRecalibrator
    # ------------------------------------------------------
    if not exists(recal_table):
        known = " ".join([f"--known-sites {k}" for k in KNOWN_SITES])

        run(f"""
        gatk --java-options "-Xmx4G" BaseRecalibrator \
        -I {dedup} \
        -R {REFERENCE} \
        {known} \
        -O {recal_table}
        """, "BaseRecalibrator")
    else:
        print(" Skipping BaseRecalibrator (recal table exists)")

    # ------------------------------------------------------
    # ApplyBQSR
    # ------------------------------------------------------
    if not exists(final):
        run(f"""
        gatk --java-options "-Xmx4G" ApplyBQSR \
        -R {REFERENCE} \
        -I {dedup} \
        --bqsr-recal-file {recal_table} \
        -O {final}
        """, "ApplyBQSR")
    else:
        print(" Skipping ApplyBQSR (final BAM exists)")

    # ------------------------------------------------------
    # Index final BAM
    # ------------------------------------------------------
    if not exists(final_bai):
        run(f"samtools index {final}", "Index final BAM")
    else:
        print(" Final BAM index exists")

    print(f" Preprocessing complete for {sample}")

############################################################
# 5 GATK MUTECT2 SOMATIC CALLING (PER-CHROM SAFE VERSION)
############################################################

def mutect2():

    merged_unfiltered = f"{OUTPUT}/merged.unfiltered.vcf.gz"
    merged_stats = f"{merged_unfiltered}.stats"
    final_vcf = f"{OUTPUT}/somatic.filtered.vcf.gz"
    pass_vcf = f"{OUTPUT}/somatic_PASS.vcf.gz"

    if exists(final_vcf):
        print("✔ Skipping Mutect2 (already completed)")
        return

    os.makedirs(f"{OUTPUT}/tmp", exist_ok=True)
    os.makedirs(f"{OUTPUT}/chr_vcfs", exist_ok=True)
    os.makedirs(f"{OUTPUT}/f1r2", exist_ok=True)

    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    pon_arg = f"--panel-of-normals {PON}" if os.path.exists(PON) else ""

    chr_vcf_files = []
    f1r2_files = []

    ############################################################
    # RUN MUTECT2 PER CHROMOSOME
    ############################################################

    for chrom in chromosomes:

        unfiltered_chr = f"{OUTPUT}/chr_vcfs/{chrom}.unfiltered.vcf.gz"
        f1r2_chr = f"{OUTPUT}/f1r2/{chrom}.f1r2.tar.gz"

        chr_vcf_files.append(unfiltered_chr)
        f1r2_files.append(f1r2_chr)

        if exists(unfiltered_chr):
            print(f"✔ Skipping {chrom}")
            continue

        run(f"""
gatk --java-options "-Xmx8G" Mutect2 \
-R {REFERENCE} \
-I {OUTPUT}/TUMOR_final.bam \
-tumor TUMOR \
-I {OUTPUT}/NORMAL_final.bam \
-normal NORMAL \
--germline-resource {GERMLINE} \
{pon_arg} \
-L {chrom} \
--native-pair-hmm-threads 2 \
--f1r2-tar-gz {f1r2_chr} \
--tmp-dir {OUTPUT}/tmp \
-O {unfiltered_chr}
""", f"Mutect2 {chrom}")

    ############################################################
    # MERGE VCFs
    ############################################################

    inputs = " ".join([f"-I {vcf}" for vcf in chr_vcf_files])

    run(f"""
gatk MergeVcfs \
{inputs} \
-O {merged_unfiltered}
""", "Merge VCFs")

    ############################################################
    # MERGE MUTECT STATS
    ############################################################

    stats_inputs = " ".join([f"-stats {vcf}.stats" for vcf in chr_vcf_files])

    run(f"""
gatk MergeMutectStats \
{stats_inputs} \
-O {merged_stats}
""", "Merge Mutect Stats")

    ############################################################
    # LEARN READ ORIENTATION MODEL (NO GATHER STEP)
    ############################################################

    f1r2_input = " ".join([f"-I {f}" for f in f1r2_files])
    read_orientation_model = f"{OUTPUT}/read-orientation-model.tar.gz"

    run(f"""
gatk LearnReadOrientationModel \
{f1r2_input} \
-O {read_orientation_model}
""", "Learn Orientation Model")

    ############################################################
    # FILTER MUTECT CALLS
    ############################################################

    run(f"""
gatk FilterMutectCalls \
-R {REFERENCE} \
-V {merged_unfiltered} \
--stats {merged_stats} \
--ob-priors {read_orientation_model} \
-O {final_vcf}
""", "FilterMutectCalls")

    ############################################################
    # SELECT PASS VARIANTS
    ############################################################

    run(f"""
gatk SelectVariants \
-R {REFERENCE} \
-V {final_vcf} \
--exclude-filtered true \
-O {pass_vcf}
""", "Select PASS Variants")

    print("✔ Mutect2 Somatic Calling Completed Successfully")

############################################################
# 6️. VEP ANNOTATION (CLINICAL-GRADE FINAL)
############################################################

def vep():

    annotated = f"{OUTPUT}/annotated.vcf"

    cache_root = f"{RESOURCES}/vep_cache"
    species_dir = f"{cache_root}/homo_sapiens"
    cache_version = "115_GRCh38"

    tar_file = f"{species_dir}/{cache_version}.tar.gz"
    extracted_dir = f"{species_dir}/{cache_version}"

    input_vcf = f"{OUTPUT}/somatic_PASS.vcf.gz"

    if exists(annotated):
        print("✔ Skipping VEP (already completed)")
        return

    # ------------------------------------------------------
    # EXTRACT CACHE IF NEEDED
    # ------------------------------------------------------
    if not os.path.exists(extracted_dir):
        if os.path.exists(tar_file):
            run(f"tar -xzf {tar_file} -C {species_dir}", "Extract VEP Cache")
        else:
            raise FileNotFoundError("❌ VEP cache missing")

    # ------------------------------------------------------
    # SAFETY: Ensure ClinVar is indexed
    # ------------------------------------------------------
    clinvar = f"{RESOURCES}/clinvar.vcf.gz"
    if os.path.exists(clinvar) and not os.path.exists(clinvar + ".tbi"):
        run(f"tabix -p vcf {clinvar}", "Index ClinVar")

    # ------------------------------------------------------
    # RUN VEP
    # ------------------------------------------------------
    run(f"""
vep \
-i {input_vcf} \
-o {annotated} \
--cache \
--dir_cache {cache_root} \
--species homo_sapiens \
--assembly GRCh38 \
--offline \
--fork {THREADS} \
--buffer_size 10000 \
--vcf \
--everything \
--hgvs \
--symbol \
--numbers \
--canonical \
--protein \
--biotype \
--uniprot \
--af \
--max_af \
--pubmed \
--variant_class \
--check_existing \
--custom {RESOURCES}/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT \
--custom {RESOURCES}/Homo_sapiens_assembly38.dbsnp138.vcf.gz,dbSNP,vcf,exact,0,ID
""", "VEP")

    print("✔ VEP Annotation Completed")

############################################################
# 7️. VCF → MAF
############################################################

def vcf_to_maf():

    maf = f"{OUTPUT}/somatic.maf"
    input_vcf = f"{OUTPUT}/annotated.vcf"

    if exists(maf):
        print("✔ Skipping vcf2maf")
        return

    run(f"""
vcf2maf.pl \
--input-vcf {input_vcf} \
--output-maf {maf} \
--tumor-id TUMOR \
--normal-id NORMAL \
--ref-fasta {REFERENCE} \
--inhibit-vep
""", "VCF2MAF")
    
############################################################
# 8. CNV CALLING (CNVkit - WES TUMOR/NORMAL)
############################################################

def cnv():

    print("\n========== CNV ANALYSIS (CNVkit - CLINICAL) ==========")

    TARGET_BED = "resources/Twist_Exome_Target_hg38.bed"
    ANNOTATION = "resources/refFlat.txt"

    # 🔥 ADD THIS (CRITICAL)
    ANTITARGET_BED = f"{OUTPUT}/cnvkit/antitarget.bed"

    tumor_bam = f"{OUTPUT}/TUMOR_final.bam"
    normal_bam = f"{OUTPUT}/NORMAL_final.bam"

    cnvkit_dir = f"{OUTPUT}/cnvkit"
    os.makedirs(cnvkit_dir, exist_ok=True)

    reference = f"{cnvkit_dir}/reference.cnn"
    tumor_cns = f"{cnvkit_dir}/TUMOR_final.cns"

    # ------------------------------------------------------
    # 0️⃣ CREATE ANTITARGET (WES FIX)
    # ------------------------------------------------------
    if not exists(ANTITARGET_BED):
        run(f"""
        cnvkit.py antitarget {TARGET_BED} \
        --output {ANTITARGET_BED}
        """, "CNVkit Antitarget Creation")

    # ------------------------------------------------------
    # 1️⃣ CREATE REFERENCE
    # ------------------------------------------------------
    if not exists(reference):
        run(f"""
        cnvkit.py batch \
        {normal_bam} \
        --normal \
        --targets {TARGET_BED} \
        --antitargets {ANTITARGET_BED} \
        --fasta {REFERENCE} \
        --annotate {ANNOTATION} \
        --output-reference {reference} \
        --output-dir {cnvkit_dir}
        """, "CNVkit Reference")

    else:
        print("✔ Skipping CNVkit reference")

    # ------------------------------------------------------
    # 2️⃣ RUN TUMOR CNV
    # ------------------------------------------------------
    if not exists(tumor_cns):
        run(f"""
        cnvkit.py batch \
        {tumor_bam} \
        --reference {reference} \
        --output-dir {cnvkit_dir}
        """, "CNVkit Tumor")

    else:
        print("✔ Skipping CNVkit tumor")

    print("✔ CNVkit Analysis Completed")

############################################################
# 8B. CNV → GENE LEVEL (CLINICAL + PRIORITIZED VERSION)
############################################################

def cnv_to_genes():

    import pandas as pd
    import numpy as np
    import os

    cns_file = f"{OUTPUT}/cnvkit/TUMOR_final.cns"
    gene_file = f"{OUTPUT}/cnvkit/TUMOR_gene_metrics.txt"
    final_output = f"{OUTPUT}/cnv_gene_level.csv"

    # ------------------------------------------------------
    # SKIP IF DONE
    # ------------------------------------------------------
    if exists(final_output) and os.path.getsize(final_output) > 0:
        print("✔ Skipping CNV gene-level annotation (already exists)")
        return

    if not exists(cns_file):
        print("⚠ CNVkit segments file not found")
        return

    # ------------------------------------------------------
    # 1️⃣ RUN GENEMETRICS
    # ------------------------------------------------------
    if not exists(gene_file):
        print(" Running CNVkit genemetrics...")
        run(f"""
        cnvkit.py genemetrics {cns_file} \
        --min-probes 3 \
        --output {gene_file}
        """, "CNVkit Gene Metrics")

    # ------------------------------------------------------
    # 2️⃣ LOAD DATA
    # ------------------------------------------------------
    try:
        df = pd.read_csv(gene_file, sep="\t")
    except Exception:
        print("⚠ Failed to read gene metrics")
        pd.DataFrame().to_csv(final_output, index=False)
        return

    # ------------------------------------------------------
    # SAFETY CHECK
    # ------------------------------------------------------
    for col in ["gene", "log2"]:
        if col not in df.columns:
            raise RuntimeError(f" Missing column: {col}")

    # ------------------------------------------------------
    # 3️⃣ CLEAN DATA
    # ------------------------------------------------------
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["log2"])
    df = df[df["log2"].between(-5, 5)]

    if "probes" in df.columns:
        df = df[df["probes"] >= 3]

    if "depth" in df.columns:
        df = df[df["depth"] >= 20]

    # ------------------------------------------------------
    # 🔥 REMOVE NOISY GENES (KEY FIX)
    # ------------------------------------------------------
    noise_patterns = [
        "RGPD","POTEB","USP","ANKRD","OR","KRTAP",
        "ZNF","RPL","RPS","LOC","LINC","MIR"
    ]
    df = df[~df["gene"].str.upper().str.contains("|".join(noise_patterns), na=False)]

    # ------------------------------------------------------
    # 4️⃣ CNV CLASSIFICATION
    # ------------------------------------------------------
    df["CNV_TYPE"] = df["log2"].apply(
        lambda x: "AMPLIFICATION" if x >= 1.0 else (
            "DELETION" if x <= -1.0 else "NEUTRAL"
        )
    )

    df["COPY_NUMBER"] = (2 * (2 ** df["log2"])).round(2)

    # ------------------------------------------------------
    # 5️⃣ REMOVE NEUTRAL
    # ------------------------------------------------------
    df = df[df["CNV_TYPE"] != "NEUTRAL"]

    # ------------------------------------------------------
    # 🔥 BIOMARKER PRIORITY (NO FILTERING)
    # ------------------------------------------------------
    important_genes = [
        "EGFR","KRAS","BRAF","ALK","ROS1","MET","RET","ERBB2",
        "TP53","RB1","PTEN","CDKN2A",
        "MYC","CCND1","CDK4","CDK6","MDM2",
        "PIK3CA","FGFR1","FGFR2","FGFR3"
    ]

    df["IMPORTANT_GENE"] = df["gene"].str.upper().isin(important_genes)

    # ------------------------------------------------------
    # 6️⃣ PRIORITY SCORING
    # ------------------------------------------------------
    df["ABS_LOG2"] = df["log2"].abs()

    df["PRIORITY_SCORE"] = (
        df["ABS_LOG2"] * 2 +
        df["IMPORTANT_GENE"].astype(int) * 5
    )

    # ------------------------------------------------------
    # 7️⃣ SORT BY CLINICAL RELEVANCE
    # ------------------------------------------------------
    df = df.sort_values(by="PRIORITY_SCORE", ascending=False)

    # ------------------------------------------------------
    # 8️⃣ TOP RESULTS LIMIT (OPTIONAL BUT USEFUL)
    # ------------------------------------------------------
    df = df.head(50)

    # ------------------------------------------------------
    # 9️⃣ HANDLE EMPTY
    # ------------------------------------------------------
    if df.empty:
        print("⚠ No CNVs detected after filtering")

        pd.DataFrame(columns=[
            "gene","chromosome","start","end",
            "log2","CNV_TYPE","COPY_NUMBER","IMPORTANT_GENE"
        ]).to_csv(final_output, index=False)

        return

    # ------------------------------------------------------
    # 🔟 PRINT IMPORTANT CNVs
    # ------------------------------------------------------
    important_df = df[df["IMPORTANT_GENE"]]

    print("\n⭐ IMPORTANT CNVs:")
    if important_df.empty:
        print("⚠ No known driver CNVs found")
    else:
        print(important_df[["gene","CNV_TYPE","COPY_NUMBER"]].head(10))

    # ------------------------------------------------------
    # FINAL SAVE
    # ------------------------------------------------------
    df.drop(columns=["ABS_LOG2","PRIORITY_SCORE"], inplace=True)

    df.to_csv(final_output, index=False)

    print(f"✔ CNV Gene-level events: {len(df)}")
    print(f"⭐ Important genes found: {df['IMPORTANT_GENE'].sum()}")
    print("✔ CNV Gene-level Annotation Completed")

############################################################
# 9. CIViC (FINAL CLINICAL-GRADE VERSION)
############################################################

def civic():
    import pandas as pd
    import requests
    import os
    import time
    import re
    from tqdm import tqdm
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    # -------------------- FILE PATHS ------------------------
    input_maf = f"{OUTPUT}/somatic.maf"
    output_file = f"{OUTPUT}/civic.merged.maf"
    civic_tsv = r"D:\Variant calling\Lung Cancer\WES\results\civic_variants.tsv"

    # -------------------- SAFETY ----------------------------
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        print("✔ Skipping CIViC (already completed)")
        return

    if not os.path.exists(input_maf):
        print("⚠ MAF file not found")
        return

    print("🌐 Running CIViC Annotation (FINAL CLINICAL-GRADE)")

    maf = pd.read_csv(input_maf, sep="\t", comment="#", low_memory=False)

    # -------------------- DETECT GENE ------------------------
    gene_col = next((c for c in ["Hugo_Symbol","Gene","SYMBOL"] if c in maf.columns), None)
    if gene_col is None:
        print("❌ Gene column missing")
        return

    maf["GENE"] = maf[gene_col].astype(str).str.upper()

    # -------------------- FILTER HIGH-QUALITY VARIANTS -------
    if "Variant_Classification" in maf.columns:
        maf = maf[
            maf["Variant_Classification"].isin([
                "Missense_Mutation","Nonsense_Mutation",
                "Frame_Shift_Del","Frame_Shift_Ins",
                "Splice_Site"
            ])
        ]

    if "FILTER" in maf.columns:
        maf = maf[maf["FILTER"] == "PASS"]

    # -------------------- AMINO ACID MAP ---------------------
    aa_map = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
        "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
        "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
        "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
    }

    def convert_aa(p):
        p = str(p).upper()

        # Convert 3-letter → 1-letter
        for k, v in aa_map.items():
            p = p.replace(k, v)

        p = p.replace("P.", "")
        p = re.sub(r'[^A-Z0-9]', '', p)

        if len(p) >= 3:
            return "P." + p
        return ""

    # -------------------- EXTRACT PROTEIN --------------------
    def extract_protein(row):
        for col in ["HGVSp_Short","HGVSp","Protein_Change"]:
            if col in row and pd.notna(row[col]):
                return convert_aa(row[col])
        return ""

    maf["PROTEIN"] = maf.apply(extract_protein, axis=1)

    print("\n🔍 Sample variants:")
    print(maf[["GENE","PROTEIN"]].head(10))

    # -------------------- DRIVER GENE FILTER -----------------
    important_genes = [
        "EGFR","KRAS","BRAF","TP53","ALK","ROS1",
        "MET","ERBB2","PIK3CA","RET","NTRK1","NTRK2","NTRK3"
    ]

    maf_filtered = maf[maf["GENE"].isin(important_genes)].copy()

    if len(maf_filtered) == 0:
        print("⚠ No driver genes → skipping CIViC")
        pd.DataFrame(columns=["GENE","PROTEIN","VARIANT","SOURCE"]).to_csv(output_file, sep="\t", index=False)
        return

    print(f"✔ Driver gene variants: {len(maf_filtered)}")

    # -------------------- LOAD CIViC TSV ---------------------
    if os.path.exists(civic_tsv):
        civic_db = pd.read_csv(civic_tsv, sep="\t", low_memory=False)
        civic_db["gene"] = civic_db["gene"].str.upper()

        civic_db["ALL"] = (
            civic_db["variant"].fillna("") + "|" +
            civic_db["variant_aliases"].fillna("") + "|" +
            civic_db["hgvs_descriptions"].fillna("")
        ).str.upper()

        civic_db.set_index("gene", inplace=True)
    else:
        print("⚠ CIViC TSV missing → API mode")
        civic_db = None

    # -------------------- API SETUP --------------------------
    session = requests.Session()
    retry = Retry(total=3, backoff_factor=1,
                  status_forcelist=[429,500,502,503,504])
    session.mount('https://', HTTPAdapter(max_retries=retry))
    gene_cache = {}

    results = []

    # -------------------- MAIN LOOP --------------------------
    for _, row in tqdm(maf_filtered.iterrows(), total=len(maf_filtered), desc="CIViC"):

        gene = row["GENE"]
        protein = row["PROTEIN"]

        if not protein:
            continue

        matched = False

        # -------- OFFLINE MATCH --------
        if civic_db is not None and gene in civic_db.index:

            matches = civic_db.loc[[gene]] if isinstance(civic_db.loc[gene], pd.DataFrame) else pd.DataFrame([civic_db.loc[gene]])

            m = matches[matches["ALL"].str.contains(protein, na=False)]

            if not m.empty:
                m = m.iloc[0]
                results.append({
                    "GENE": gene,
                    "PROTEIN": protein,
                    "VARIANT": m["variant"],
                    "SOURCE": "CIViC_OFFLINE"
                })
                matched = True

        # -------- API MATCH --------
        if not matched:
            try:
                if gene not in gene_cache:
                    url = f"https://civicdb.org/api/variants?entrez_name={gene}"
                    r = session.get(url, timeout=10)
                    gene_cache[gene] = r.json().get("records", []) if r.status_code == 200 else []
                    time.sleep(0.1)

                for var in gene_cache[gene]:
                    name = str(var.get("name","")).upper()
                    clean_name = re.sub(r'[^A-Z0-9]', '', name)

                    if protein.replace("P.","") in clean_name:
                        results.append({
                            "GENE": gene,
                            "PROTEIN": protein,
                            "VARIANT": name,
                            "SOURCE": "CIViC_API"
                        })
                        matched = True
                        break

            except Exception:
                continue

    # -------------------- EMPTY ------------------------------
    if len(results) == 0:
        print("⚠ No CIViC matches found")
        pd.DataFrame(columns=["GENE","PROTEIN","VARIANT","SOURCE"]).to_csv(output_file, sep="\t", index=False)
        return

    df = pd.DataFrame(results).drop_duplicates()

    # -------------------- SAVE -------------------------------
    df.to_csv(output_file, sep="\t", index=False)

    print(f"✔ CIViC matches: {len(df)}")
    print("✔ CIViC Annotation Completed (FINAL CLINICAL-GRADE)")
    
############################################################
# 10. AMP / ASCO / CAP Clinical Classification
############################################################

def amp_asco_classification():

    input_file = f"{OUTPUT}/civic.maf"
    output_file = f"{OUTPUT}/clinical_scored_variants.csv"

    # 🔥 FIX: check empty file
    if not exists(input_file) or os.path.getsize(input_file) == 0:
        print("⚠ CIViC file empty → skipping AMP classification")
        return

    try:
        df = pd.read_csv(input_file, sep="\t")
    except Exception:
        print("⚠ CIViC file unreadable → skipping AMP classification")
        return

    if df.empty:
        print("⚠ No data for AMP classification")
        return

    # -------------------------------
    # Tier Assignment
    # -------------------------------
    def assign_tier(row):

        level = str(row.get("EVIDENCE_LEVEL", "")).upper()
        significance = str(row.get("CLINICAL_SIGNIFICANCE", "")).lower()

        if level == "A":
            return "Tier I"

        if level == "B":
            if "sensitivity" in significance or "responsive" in significance:
                return "Tier I"
            return "Tier II"

        if level == "C":
            return "Tier II"

        if level == "D":
            return "Tier III"

        if level == "E":
            return "Tier IV"

        return "Tier III"

    df["AMP_TIER"] = df.apply(assign_tier, axis=1)

    # -------------------------------
    # Clinical Interpretation Label
    # -------------------------------
    def interpretation(tier):

        return {
            "Tier I": "Strong clinical significance",
            "Tier II": "Potential clinical significance",
            "Tier III": "Variant of unknown significance",
            "Tier IV": "Likely benign"
        }.get(tier, "Unknown")

    df["CLINICAL_INTERPRETATION"] = df["AMP_TIER"].apply(interpretation)

    df.to_csv(output_file, index=False)

    print("✔ AMP/ASCO Classification Completed")

############################################################
# 11. Actionable Drugs
############################################################

def actionable_drugs():

    input_maf = f"{OUTPUT}/civic.maf"

    strict_output = f"{OUTPUT}/actionable_drugs_strict.csv"
    extended_output = f"{OUTPUT}/actionable_drugs_extended.csv"
    resistance_output = f"{OUTPUT}/resistance_variants.csv"

    if not exists(input_maf):
        print("⚠ CIViC MAF not found")
        return

    df = pd.read_csv(input_maf, sep="\t")

    if df.empty:
        print("⚠ No actionable variants found")
        return

    df["EVIDENCE_DIRECTION"] = df.apply(
        lambda x: "Resistant" if "resistance" in str(x["CLINICAL_SIGNIFICANCE"]).lower()
        else "Sensitive",
        axis=1
    )

    strict = df[df["ONCOKB_LEVEL"].isin(["LEVEL_1", "LEVEL_2", "LEVEL_3A"])]
    extended = df[df["ONCOKB_LEVEL"].isin(["LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B"])]
    resistance = df[df["EVIDENCE_DIRECTION"] == "Resistant"]

    clinical_cols = [
        "GENE",
        "PROTEIN_CHANGE",
        "VARIANT_CLASSIFICATION",
        "TUMOR_TYPE",
        "ONCOKB_LEVEL",
        "DRUG",
        "EVIDENCE_DIRECTION"
    ]

    strict[clinical_cols].drop_duplicates().to_csv(strict_output, index=False)
    extended[clinical_cols].drop_duplicates().to_csv(extended_output, index=False)
    resistance[clinical_cols].drop_duplicates().to_csv(resistance_output, index=False)

    print(f"✔ Strict actionable variants: {len(strict)}")
    print(f"✔ Extended actionable variants: {len(extended)}")
    print(f"✔ Resistance variants: {len(resistance)}")

############################################################
# 12. TUMOR MUTATIONAL BURDEN (TMB)
############################################################

def tmb():

    vcf_file = f"{OUTPUT}/somatic.filtered.vcf.gz"
    capture_size_mb = 30  # <-- change if your WES kit differs

    if not exists(vcf_file):
        print("⚠ Filtered VCF not found")
        return

    import gzip

    count = 0
    with gzip.open(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[6] == "PASS":
                count += 1

    tmb_value = count / capture_size_mb

    with open(f"{OUTPUT}/TMB.txt", "w") as out:
        out.write(f"{tmb_value:.2f} mutations/Mb")

    print(f"✔ TMB: {tmb_value:.2f} mutations/Mb")

############################################################
# 13. Final HTML report
############################################################

def report():

    report_file = f"{OUTPUT}/clinical_report.html"

    # --------------------------------------------------
    # LOAD DATA
    # --------------------------------------------------

    tmb_value = open(f"{OUTPUT}/TMB.txt").read() if exists(f"{OUTPUT}/TMB.txt") else "Not Available"

    # SNV
    snv_file = f"{OUTPUT}/actionable_drugs_strict.csv"
    if exists(snv_file):
        snv_df = pd.read_csv(snv_file)
        snv_table = snv_df.to_html(index=False, classes="table table-striped")
        snv_count = len(snv_df)
    else:
        snv_table = "<p>No clinically actionable SNVs identified.</p>"
        snv_count = 0

    # CNV
    cnv_file = f"{OUTPUT}/cnv_gene_level.csv"
    if exists(cnv_file):
        cnv_df = pd.read_csv(cnv_file)
        amp = cnv_df[cnv_df["CNV_TYPE"] == "AMPLIFICATION"]
        dele = cnv_df[cnv_df["CNV_TYPE"] == "DELETION"]

        amp_table = amp.head(20).to_html(index=False, classes="table")
        del_table = dele.head(20).to_html(index=False, classes="table")

        amp_count = len(amp)
        del_count = len(dele)
    else:
        amp_table = "<p>No amplifications detected</p>"
        del_table = "<p>No deletions detected</p>"
        amp_count = del_count = 0

    # Resistance
    resistance_file = f"{OUTPUT}/resistance_variants.csv"
    if exists(resistance_file):
        res_df = pd.read_csv(resistance_file)
        res_table = res_df.to_html(index=False, classes="table")
    else:
        res_table = "<p>No resistance-associated variants detected</p>"

    # AMP Classification
    amp_file = f"{OUTPUT}/clinical_scored_variants.csv"
    if exists(amp_file):
        amp_df = pd.read_csv(amp_file)
        tier_counts = amp_df["AMP_TIER"].value_counts().to_dict()
    else:
        tier_counts = {}

    # CancerVar
    cancervar_file = f"{OUTPUT}/cancervar_results.txt"
    cancervar_section = "<p>CancerVar results available as supplementary output.</p>"
    if exists(cancervar_file):
        cancervar_section = "<p>Variants interpreted using CancerVar evidence-based scoring.</p>"

    # --------------------------------------------------
    # TMB INTERPRETATION
    # --------------------------------------------------
    try:
        tmb_numeric = float(tmb_value.split()[0])
        if tmb_numeric >= 10:
            tmb_label = "High (Potential immunotherapy benefit)"
        elif tmb_numeric >= 5:
            tmb_label = "Intermediate"
        else:
            tmb_label = "Low"
    except:
        tmb_label = "Not Interpretable"

    # --------------------------------------------------
    # HTML REPORT
    # --------------------------------------------------

    html = f"""
    <html>
    <head>
        <title>Clinical Oncology Report</title>
        <style>
            body {{ font-family: Arial; margin: 40px; background: #f4f6f7; }}
            h1 {{ color: #1f618d; }}
            h2 {{ border-bottom: 2px solid #ccc; padding-bottom: 5px; }}
            .section {{ margin-bottom: 30px; background: white; padding: 20px; border-radius: 10px; box-shadow: 0 0 8px #ccc; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ccc; padding: 8px; }}
            th {{ background-color: #2e86c1; color: white; }}
            .highlight {{ color: #c0392b; font-weight: bold; }}
        </style>
    </head>

    <body>

        <h1>🧬 Clinical Genomic Oncology Report</h1>

        <!-- SUMMARY -->
        <div class="section">
            <h2>📊 Clinical Summary</h2>
            <p><strong>Tumor Mutational Burden:</strong> {tmb_value} 
            (<span class="highlight">{tmb_label}</span>)</p>
            <p><strong>Actionable SNVs:</strong> {snv_count}</p>
            <p><strong>CNV Amplifications:</strong> {amp_count}</p>
            <p><strong>CNV Deletions:</strong> {del_count}</p>
        </div>

        <!-- SNV -->
        <div class="section">
            <h2>🧪 Clinically Actionable Variants (SNV/Indel)</h2>
            {snv_table}
        </div>

        <!-- CNV -->
        <div class="section">
            <h2>📈 Copy Number Alterations</h2>

            <h3>Amplifications</h3>
            {amp_table}

            <h3>Deletions</h3>
            {del_table}
        </div>

        <!-- RESISTANCE -->
        <div class="section">
            <h2>⚠️ Drug Resistance Markers</h2>
            {res_table}
        </div>

        <!-- AMP -->
        <div class="section">
            <h2>🏥 Clinical Evidence Classification (AMP/ASCO/CAP)</h2>
            <ul>
                <li>Tier I (Strong Clinical Significance): {tier_counts.get("Tier I", 0)}</li>
                <li>Tier II (Potential Clinical Significance): {tier_counts.get("Tier II", 0)}</li>
                <li>Tier III (Unknown Significance): {tier_counts.get("Tier III", 0)}</li>
                <li>Tier IV (Benign/Likely Benign): {tier_counts.get("Tier IV", 0)}</li>
            </ul>
        </div>

        <!-- CancerVar -->
        <div class="section">
            <h2>🧠 Variant Interpretation</h2>
            {cancervar_section}
        </div>

        <!-- METHODS -->
        <div class="section">
            <h2>⚙️ Bioinformatics Pipeline</h2>
            <ul>
                <li>Alignment: BWA-MEM</li>
                <li>Variant Calling: GATK Mutect2</li>
                <li>Variant Annotation: VEP (Ensembl)</li>
                <li>Clinical Annotation: CIViC + CancerVar</li>
                <li>Copy Number Analysis: CNVkit</li>
                <li>Variant Filtering: PASS somatic mutations</li>
            </ul>
        </div>

        <!-- INTERPRETATION -->
        <div class="section">
            <h2>🧾 Clinical Interpretation</h2>
            <p>
            This analysis identifies somatic alterations with potential clinical relevance. 
            Variants classified as Tier I/II may guide targeted therapy selection, while 
            resistance variants indicate possible therapy failure mechanisms.
            </p>
        </div>

        <!-- DISCLAIMER -->
        <div class="section">
            <h2>⚠️ Disclaimer</h2>
            <p>
            This report is intended for research and clinical decision support only. 
            Final clinical decisions should be made by a qualified oncologist 
            in conjunction with patient-specific data.
            </p>
        </div>

    </body>
    </html>
    """

    with open(report_file, "w") as f:
        f.write(html)

    print("✔ Clinical-Grade Oncology Report Generated")

############################################################
# MAIN EXECUTION
############################################################

def main():

    start = time.time()

    # -----------------------------
    # PREPROCESSING
    # -----------------------------
    prepare_reference()

    fastp("TUMOR", TUMOR_R1, TUMOR_R2)
    fastp("NORMAL", NORMAL_R1, NORMAL_R2)

    align("TUMOR")
    align("NORMAL")

    preprocess("TUMOR")
    preprocess("NORMAL")

    # -----------------------------
    # VARIANT CALLING
    # -----------------------------
    mutect2()

    # -----------------------------
    # ANNOTATION
    # -----------------------------
    vep()
    vcf_to_maf()

    # -----------------------------
    # CNV
    # -----------------------------
    cnv()
    cnv_to_genes()

    #  CIViC
    civic()
    
    # AMP / ASCO / CAP Clinical Classification
    amp_asco_classification()

    # -----------------------------
    # CLINICAL INTERPRETATION
    # -----------------------------
    actionable_drugs()

    # -----------------------------
    # DOWNSTREAM ANALYSIS
    # -----------------------------
    tmb()
    report()

    # -----------------------------
    # LOGGING
    # -----------------------------
    total_hours = (time.time() - start) / 3600

    log(f"\n PIPELINE COMPLETE")
    log(f"Total runtime: {total_hours:.2f} hours")

    print("\n PIPELINE COMPLETE")
    print("⏱ Hours:", total_hours)


if __name__ == "__main__":
    main()
