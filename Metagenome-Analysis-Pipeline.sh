#!/bin/bash
#SBATCH --job-name=metagenome_pipeline
#SBATCH --cpus-per-task=40
#SBATCH --mem=500G
#SBATCH --time=72:00:00

# --- Path and Environment Configuration ---
DATA_DIR="/path/to/raw"
OUT_DIR="/path/to/output"
DB_PATH="/path/to/DB"
THREADS=40
MEMORY=500

# Create output directory structure
mkdir -p $OUT_DIR/{1-clean,2-assembly,3-gene_catalog,4-tax_kegg,5-ARGs_OAP_reads,6-BIN}

#------------ [ Step 1: Quality Control (Fastp) ] ------------
# Filtering: N > 10%, low quality bases > 50%
while IFS=$'\t' read -r SAMPLE R1 R2; do
    fastp -i "$R1" -I "$R2" \
          -o "$OUT_DIR/1-clean/${SAMPLE}_R1.fq.gz" -O "$OUT_DIR/1-clean/${SAMPLE}_R2.fq.gz" \
          --n_base_limit 10 --unqualified_percent_limit 50 --qualified_quality_phred 10 \
          -w $THREADS
    echo -e "$SAMPLE\t$OUT_DIR/1-clean/${SAMPLE}_R1.fq.gz\t$OUT_DIR/1-clean/${SAMPLE}_R2.fq.gz" >> "$OUT_DIR/1-clean/clean.path"
done < sample_list.txt

#------------ [ Step 2: Read-based Analysis (ARGs-OAP 3.2) ] ------------
# Stage_one: Extract ARG-like reads and estimate Cell Number based on essential single-copy genes
mkdir -p "$OUT_DIR/5-ARGs_OAP_reads/stage_one"
args_oap stage_one -i "$OUT_DIR/1-clean" -o "$OUT_DIR/5-ARGs_OAP_reads/stage_one" -f fastq -t $THREADS

# Stage_two: Annotation for SARG, VFDB, and MobileOG databases
# Thresholds: Identity >=80%, Query Coverage >=75%, E-value <=1e-7
for DB in "SARG" "VFDB" "MobileOG"; do
    diamond blastx -d "$DB_PATH/${DB}.dmnd" \
                   -q "$OUT_DIR/5-ARGs_OAP_reads/stage_one/extracted.fa" \
                   --evalue 1e-7 --id 80 --query-cover 75 \
                   -f 6 qseqid sseqid pident length evalue bitscore \
                   -o "$OUT_DIR/5-ARGs_OAP_reads/${DB}_reads_hits.txt" -p $THREADS
done

#------------ [ Step 3: Assembly & Gene Catalog ] ------------
# Contig assembly using MEGAHIT
while IFS=$'\t' read -r SAMPLE R1 R2; do
    megahit -1 "$R1" -2 "$R2" -o "$OUT_DIR/2-assembly/$SAMPLE" -t $THREADS --min-contig-len 500
    cat "$OUT_DIR/2-assembly/$SAMPLE"/*.contigs.fa >> "$OUT_DIR/2-assembly/all_combined.fa"
done < "$OUT_DIR/1-clean/clean.path"

# ORF prediction (Prodigal) and Gene De-replication (CD-HIT)
# Thresholds: 95% Identity, 90% Coverage
prodigal -i "$OUT_DIR/2-assembly/all_combined.fa" -p meta \
         -a "$OUT_DIR/3-gene_catalog/all_pro.faa" \
         -d "$OUT_DIR/3-gene_catalog/all_nuc.fna"

cd-hit-est -i "$OUT_DIR/3-gene_catalog/all_nuc.fna" \
           -o "$OUT_DIR/3-gene_catalog/gene_catalog.fna" \
           -c 0.95 -aS 0.9 -M 0 -T $THREADS

#------------ [ Step 4: Taxonomy & KEGG (Contig-based) ] ------------
# 1. Taxonomic identification via NCBI NR and MEGAN LCA
diamond blastp -d "$DB_PATH/nr.dmnd" -q "$OUT_DIR/3-gene_catalog/all_pro.faa" \
               --evalue 1e-5 --max-target-seqs 50 -f 100 \
               -o "$OUT_DIR/4-tax_kegg/nr.daa" -p $THREADS

# 2. Functional annotation via KEGG
# Thresholds: E-value <= 1e-5, HSP score > 60 bits
diamond blastp -d "$DB_PATH/kegg.dmnd" -q "$OUT_DIR/3-gene_catalog/all_pro.faa" \
               --evalue 1e-5 --min-score 60 -f 6 \
               -o "$OUT_DIR/4-tax_kegg/kegg.txt" -p $THREADS

# 3. TPM Quantification using Bowtie2 and SAMtools
bowtie2-build "$OUT_DIR/3-gene_catalog/gene_catalog.fna" "$OUT_DIR/3-gene_catalog/bt2_idx"
while IFS=$'\t' read -r SAMPLE R1 R2; do
    bowtie2 -x "$OUT_DIR/3-gene_catalog/bt2_idx" -1 "$R1" -2 "$R2" -p $THREADS | \
    samtools sort -@ $THREADS -o "$OUT_DIR/3-gene_catalog/${SAMPLE}.bam"
    samtools idxstats "$OUT_DIR/3-gene_catalog/${SAMPLE}.bam" > "$OUT_DIR/3-gene_catalog/${SAMPLE}.counts"
done < "$OUT_DIR/1-clean/clean.path"

#------------ [ Step 5: Genome Binning (MetaWRAP) ] ------------
# Standard MetaWRAP pipeline: Binning -> Refinement -> GTDB-Tk
metawrap binning -a "$OUT_DIR/2-assembly/all_combined.fa" -o "$OUT_DIR/6-BIN/bins" \
                 -t $THREADS -m $MEMORY --metabat2 --maxbin2 --concoct "$OUT_DIR/1-clean"/*.fq.gz

# Refinement: Completeness > 50%, Contamination < 10%
metawrap bin_refinement -o "$OUT_DIR/6-BIN/refined" -t $THREADS -c 50 -x 10 \
                        -A "$OUT_DIR/6-BIN/bins/metabat2_bins" \
                        -B "$OUT_DIR/6-BIN/bins/maxbin2_bins" \
                        -C "$OUT_DIR/6-BIN/bins/concoct_bins"

# Taxonomic classification of MAGs using GTDB-Tk
gtdbtk classify_wf --genome_dir "$OUT_DIR/6-BIN/refined/metawrap_50_10_bins" \
                   --out_dir "$OUT_DIR/6-BIN/gtdbtk" --cpus $THREADS -x fna

# ARG annotation within MAGs (SARG database)
# Mode: Sensitive, Thresholds: ID >= 80%, Coverage >= 70%, E-value <= 1e-7
for mag in "$OUT_DIR/6-BIN/refined/metawrap_50_10_bins"/*.fna; do
    bn=$(basename "$mag" .fna)
    prodigal -i "$mag" -a "$OUT_DIR/6-BIN/refined/${bn}.faa" -p meta
    diamond blastp --sensitive -d "$DB_PATH/SARG.dmnd" -q "$OUT_DIR/6-BIN/refined/${bn}.faa" \
                   --evalue 1e-7 --id 80 --query-cover 70 \
                   -o "$OUT_DIR/6-BIN/refined/${bn}_sarg.txt" -p $THREADS
done

#------------ [ Step 6: Final Data Summarization ] ------------
echo "Summarizing results into a Master Matrix..."
python3 summarize_metagenome.py