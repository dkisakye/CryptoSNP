# Define all samples (replace with your actual sample names)
SAMPLES = ['Strain001', 'Strain002', 'Strain003']

# Rule: All final output files
rule all:
    input:
        expand("{OUT_DIR}/fastqc_reports/{sample}_fastqc.html", sample=SAMPLES),
        expand("{BAM_DIR}/{sample}.mapq20.RGs.bam.bai", sample=SAMPLES),
        expand("{OUT_DIR_VCF}/Fv_{contig}_Raw_Variants.vcf", contig=CONTIGS)

# Define variables
READS_DIR = "/path/to/raw/reads"
OUT_DIR = "/path/to/output"
BWA_THREADS = 8
BAM_DIR = "/path/to/output/bam"
UNMAPPED = "/path/to/output/unmapped_reads"
STATS_DIR = "/path/to/output/bam_stats"
GENOME_IDX = "/path/to/genome/index.mmi"
REF = "/path/to/reference.fasta"
OUT_DIR_VCF = "/path/to/output/VCF"
CTG_LIST = "/path/to/chromosome/list.txt"
BAM_FOF = "/path/to/bam/file-of-files.txt"
ADAPTERS = "/path/to/adapters.fa"

N_PAR = 4  # Adjust the number of parallel jobs

# Rule: Run FastQC
rule fastqc:
    input:
        $(find "${READS_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*.fastq.gz')
    output:
        directory("{OUT_DIR}/fastqc_reports")
    shell:
        """
        mkdir -p {OUT_DIR}/fastqc_reports
        fastqc -f {input} -o {output} -t {threads} {input}
        """

# Rule: Run MultiQC
rule multiqc:
    input:
        fastqc_reports=expand("{OUT_DIR}/fastqc_reports/{sample}_fastqc.{ext}",
                            sample=SAMPLES,
                            ext=["zip", "html"])
    output:
        multiqc_report="{OUT_DIR}/multiqc_report.html"
    shell:
        """
        multiqc --outdir {OUT_DIR} --export --flat {OUT_DIR}/fastqc_reports
        """

# Rule: Trim reads
rule trim_reads:
    input:
        fwd_read=lambda wildcards: get_fwd_read(wildcards.sample, READS_DIR),
        rev_read=lambda wildcards: get_rev_read(wildcards.sample, READS_DIR)
    output:
        fwd_trimmed="{OUT_DIR}/{sample}_1P.fastq.gz",
        fwd_untrimmed="{OUT_DIR}/{sample}_1U.fastq.gz",
        rev_trimmed="{OUT_DIR}/{sample}_2P.fastq.gz",
        rev_untrimmed="{OUT_DIR}/{sample}_2U.fastq.gz",
        summary="{LOG_DIR}/{sample}_summary.txt"
    params:
        adapters=ADAPTERS,
        trim_threads=TRIM_THREADS
    shell:
        """
        java -jar {TRIMMOMATIC}/trimmomatic-0.39.jar PE \
            -threads {params.trim_threads} \
            -summary {output.summary} \
            {input.fwd_read} \
            {input.rev_read} \
            {output.fwd_trimmed} \
            {output.fwd_untrimmed} \
            {output.rev_trimmed} \
            {output.rev_untrimmed} \
            ILLUMINACLIP:{params.adapters}:4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
        """

# Helper functions to get forward and reverse read paths
def get_fwd_read(wildcards, READS_DIR):
    fwd_read = os.path.join(READS_DIR, f"{wildcards.sample}_R1.fastq.gz")
    return fwd_read

def get_rev_read(wildcards, READS_DIR):
    rev_read = os.path.join(READS_DIR, f"{wildcards.sample}_R2.fastq.gz")
    return rev_read

# Rule: Map reads using bwa mem
rule map_reads:
    input:
        fwd_read="{READS_DIR}/{sample}_1P.fastq.gz",
        rev_read="{READS_DIR}/{sample}_2P.fastq.gz"
    output:
        bam="{BAM_DIR}/{sample}.bam"
    params:
        genome_idx=GENOME_IDX,
        bwa_threads=BWA_THREADS
    shell:
        """
        bwa mem \
            -t {params.bwa_threads} \
            -k 15 \
            -r 1.0 \
            -M \
            {params.genome_idx} \
            {input.fwd_read} {input.rev_read} \
            | samtools view -hbu \
            | samtools sort -o {output.bam} -
        """

# Rule: Remove unmapped reads
rule remove_unmapped:
    input:
        bam="{BAM_DIR}/{sample}.bam"
    output:
        unmapped_fwd="{UNMAPPED}/{sample}_Unmapped_R1.fastq.gz",
        unmapped_rev="{UNMAPPED}/{sample}_Unmapped_R2.fastq.gz",
        unmapped_single="{UNMAPPED}/{sample}_Unmapped_U.fastq.gz"
    params:
        bwa_threads=BWA_THREADS
    shell:
        """
        samtools view -f 4 -hbu {input.bam} \
            | samtools sort -@ {params.bwa_threads} -n -o /{UNMAPPED}/{wildcards.sample}.unmapped.qsort.bam -
        samtools fastq \
            -1 {output.unmapped_fwd} \
            -2 {output.unmapped_rev} \
            -s {output.unmapped_single} \
            -0 /dev/null \
            -@ {params.bwa_threads} \
            /{UNMAPPED}/{wildcards.sample}.unmapped.qsort.bam
        """

# Rule: Filter by mapping quality
rule filter_mapq:
    input:
        bam="{BAM_DIR}/{sample}.bam"
    output:
        mapq_filtered="{BAM_MAPQ}/{sample}.mapq20.bam"
    shell:
        """
        samtools view \
            -hb -q 20 \
            {input.bam} \
            > {output.mapq_filtered}
        """

# Rule: Index BAM files
rule index_bam:
    input:
        bam="{BAM_DIR}/{sample}.bam"
    output:
        "{input.bam}.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Rule: Add read groups to BAM files
rule add_read_groups:
    input:
        bam="{BAM_DIR}/{sample}.mapq20.bam"
    output:
        bam="{BAM_DIR}/{sample}.mapq20.RGs.bam",
        bai="{BAM_DIR}/{sample}.mapq20.RGs.bam.bai"
    shell:
        """
        java -jar picard.jar AddOrReplaceReadGroups \
        I={input.bam} \
        O={output.bam} \
        RGID={sample} \
        RGLB={sample} \
        RGPL=ILLUMINA \
        RGPU={sample} \
        RGSM={sample} \
        CREATE_INDEX=true
        """

# Rule: Mark duplicates
rule mark_duplicates:
    input:
        bam="{BAM_DIR}/{sample}.mapq20.bam"
    output:
        marked_bam="{BAM_DIR}/{sample}.mapq20.MarkDups.bam",
        metrics="{STATS_DIR}/{sample}.MarkDups.metrics.txt"
    
    shell:
        """
        picard MarkDuplicates \
            I={input.bam} \
            O={output.marked_bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false
        """
# Rule:Generate MLST consensus sequences
# Define MLST loci
CAP59 = "CP003820.1:1879601-1880160"
GPD1 = "CP003826.1:523911-524454"
IGS1 = "CP003821.1:280331-281053"
LAC1 = "CP003827.1:1022213-1022682"
PLB1 = "CP003831.1:260896-261428"
SOD1 = "CP003824.1:1416715-1417250"
URA5 = "CP003827.1:315161-315797"

rule generate_consensus:
    input:
        bam="{BAM_DIR}/{sample}.mapq20.RGs.MarkDups.bam"
    output:
        consensus_fasta="MLST/{sample}_consensus.fasta"
    params:
        cap59=CAP59,
        gpd1=GPD1,
        igs1=IGS1,
        lac1=LAC1,
        plb1=PLB1,
        sod1=SOD1,
        ura5=URA5

    shell:
        """
        # Get the consensus sequences for each locus
        cap59=$(samtools consensus -r {params.cap59} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')
        gpd1=$(samtools consensus -r {params.gpd1} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')
        igs1=$(samtools consensus -r {params.igs1} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')
        lac1=$(samtools consensus -r {params.lac1} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')
        plb1=$(samtools consensus -r {params.plb1} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')
        sod1=$(samtools consensus -r {params.sod1} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')
        ura5=$(samtools consensus -r {params.ura5} -m "simple" -c 0.5 -H 0.5 -a {input.bam} | tail -n +2 | tr -d '\n')

        # Print the consensus sequence to the output FASTA file
        echo ">{wildcards.sample}" > {output.consensus_fasta}
        echo "${cap59}${gpd1}${igs1}${lac1}${plb1}${sod1}${ura5}" >> {output.consensus_fasta}
        """


# Rule: Variant calling with Freebayes
rule variant_calling:
    input:
        bam=expand("{BAM_DIR}/{sample}.mapq20.RGs.bam", sample=SAMPLES)
    output:
        vcf="{OUT_DIR_VCF}/Fv_{contig}_Raw_Variants.vcf"
    params:
        bam_fof=BAM_FOF,
        ref=REF,
        out_dir=OUT_DIR_VCF,
        ctg_list=CTG_LIST
    shell:
        """
        module load parallel/20210822
        module load samtools/1.14
        module load java/openjdk-17.0.2

        mkdir -p {params.out_dir}

        find {BAM_DIR} -mindepth 1 -maxdepth 1 -type f -name '*.bam' | \
        grep -vE {params.exclude} | sort -V > {params.bam_fof}

        for contig in $(cat {params.ctg_list})
        do
            freebayes \
            -f {params.ref} \
            -L {params.bam_fof} \
            -r $contig \
            -v {output.vcf} \
            --strict-vcf \
            --use-best-n-alleles 4 \
            -E -1 \
            -q 20 \
            -R 40 \
            --min-coverage 50 \
            --limit-coverage 1000 \
            --no-population-priors
        done
        """
