import os
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/config.sh; set -eo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "%s/config.yaml" % SNAKEMAKE_DIR

REFERENCE = config["reference"]

BARCODE_FASTA = config["barcode_fasta"]
BARCODE_PREFIX = config["barcode_prefix"]

MRSFAST_BINARY = config["mrsfast_path"]
MRSFAST_OPTS = config["mrsfast_opts"]
MRSFAST_INDEX = config[REFERENCE]["mrsfast_index"]
CONTIGS = config[REFERENCE]["contigs"]
SUNK_MASK = config[REFERENCE]["sunk_mask"]

MANIFEST_FILE = config["manifest"]

MANIFEST = pd.read_table(MANIFEST_FILE)
MANIFEST.index = MANIFEST.sample_name

if not os.path.exists("log"):
    os.makedirs("log")

def get_well_split_fq_from_sample(wildcards):
    well = MANIFEST.loc[wildcards.sample, "well"]
    return ["mapping/{well}/{well}/fastq_split/{well}.{num}_part0.fastq.gz".format(well=well, num=num) for num in [1, 2]]

localrules: all

rule all:
    input: "clone_locations.bed"

rule get_pileup_locations:
    input: expand("sunk_pileup/{sample}.sorted.bam_sunk.bw", sample=MANIFEST.sample_name)
    output: "clone_locations.bed"
    params: sge_opts = "-l mfree=1G -l h_rt=1:00:00"
    benchmark: "benchmarks/clone_loc/{sample}.txt"
    run:
        for file in input:
            fn = os.path.basename(file)
            shell("""bigWigToBedGraph {file} /dev/stdout \
            | awk 'OFS="\t" {{ print $1,$2,$3,".",$4 }}' \
            | bedtools merge -i stdin -d 100000 -c 5 -o sum | sort -k 4,4rn \
            | head -n 1 | awk 'OFS="\t" {{ print $1,$2,$3,$4,"{fn}" }}' >> {output}""")

rule make_bw_pileup:
    input: "bam/{sample}.sorted.bam"
    output: "sunk_pileup/{sample}.sorted.bam_sunk.bw"
    params: sge_opts="-N plp_{sample} -l h_rt=0:20:00 -l mfree=2G"
    benchmark: "benchmarks/pileup/{sample}.txt"
    shell:
        "python /net/eichler/vol2/local/inhousebin/sunks/pileups/sam_to_bw_pileup.py "
        "--inSam {SNAKEMAKE_DIR}/{input} --contigs {CONTIGS} "
        "--outdir {SNAKEMAKE_DIR}/sunk_pileup "
        "--sunk_mask {SUNK_MASK} --track_url test"

rule make_bam:
    input: "mapping/{sample}/{sample}/mrfast_out/{sample}.sam.gz"
    output: "bam/{sample}.sorted.bam", "bam/{sample}.sorted.bam.bai"
    params: sge_opts="-N bam_{sample} -l mfree=1G -l h_rt=0:30:00",
            output_prefix="bam/{sample}.sorted"
    benchmark: "benchmarks/make_bam/{sample}.txt"
    shell:
        "zcat {input} | samtools view -b -t {CONTIGS} - -S | samtools sort - -T $TMPDIR/{wildcards.sample} -o {output[0]}; "
        "samtools index {output[0]}"

rule map:
    input: get_well_split_fq_from_sample
    output: "mapping/{sample}/{sample}/mrfast_out/{sample}.sam.gz"
    params: sge_opts="-N map_{sample} -l mfree=4G -l h_rt=1:0:0:0",
            output_prefix="mapping/{sample}/{sample}/mrfast_out/{sample}"
    benchmark: "benchmarks/map/{sample}.txt"
    shell:
        "zcat {input} | {MRSFAST_BINARY} --search {MRSFAST_INDEX} {MRSFAST_OPTS} --seq /dev/stdin -o {params.output_prefix} --disable-nohit"

rule split_fastq:
    input: "merged_megapool_lanes_by_barcode/{well}/{well}.{num}.fastq.gz"
    output: "mapping/{well}/{well}/fastq_split/{well}.{num}_part0.fastq.gz"
    params: sge_opts="-N split_{well} -q eichler-short.q -l h_rt=6:00:00 -pe orte 5-10 -l disk_free=10G", 
            input_dir="%s/split_barcodes/{well}/{well}/fastq" % SNAKEMAKE_DIR,
            output_dir="%s/mapping/{well}/{well}/fastq_split" % SNAKEMAKE_DIR
    benchmark: "benchmarks/split_fastq/{sample}.txt"
    run: 
        fn = os.path.basename(input[0])
        if input[0].endswith(".bam"):
            infile = "$TMPDIR/%s" % (fn.replace(".bam", ".fastq"))
            shell("""bamToFastq -i {input} -fq /dev/stdout -fq2 /dev/stdout > {infile}; bgzip {infile}""")
            infile += ".gz"
        else:
            infile = "$TMPDIR/%s" % fn
            shell("rsync --bwlimit=50000 {input} {infile}")
        of = os.path.basename(output[0])
        shell("""mpirun -x PATH -x LD_LIBRARY_PATH --prefix $MPIBASE -mca plm ^rshd -mca btl ^openib 
                 /net/eichler/vol4/home/a5ko/bin/readSplit -s 36 -k 36 -n 1000000 -i {infile} -o $TMPDIR; 
                 rsync --bwlimit=50000 $TMPDIR/{of} {output}""")

rule demultiplex_fastq:
    input: ["fastq/%s%d.fq.gz" % (BARCODE_PREFIX, num) for num in [1,2,3]]
    output: "split_barcodes/{well}/{well}/fastq/{well}.1.fastq.gz", "split_barcodes/{well}/{well}/fastq/{well}.2.fastq.gz"
    params: sge_opts=""
    shell:
        "/net/eichler/vol7/home/psudmant/EEE_Lab/projects/fastq_utils/code/parse_barcodes_new/parse_barcodes.py "
        "--barcodes_fa {BARCODE_FASTA} --input_fastq {input[0]} --input_fastq_pe2 {input[2]} "
        "--barcode_fastq {input[1]} --output_dir split_barcodes/"
