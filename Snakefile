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
SAMTOOLS_CONTIGS = config[REFERENCE]["samtools_contigs"]

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
    input: expand("bam/{sample}.sorted.{ext}", sample=MANIFEST.sample_name, ext=["bam", "bam.bai"])

rule make_bam:
    input: "mapping/{sample}/{sample}/mrfast_out/{sample}.sam.gz"
    output: "bam/{sample}.sorted.bam", "bam/{sample}.sorted.bam.bai"
    params: sge_opts="-l mfree=1G -l h_rt=0:30:00", output_prefix="bam/{sample}.sorted"
    shell:
        "zcat {input} | samtools view -b -t {SAMTOOLS_CONTIGS} - -S | samtools sort - {params.output_prefix}; "
        "samtools index {output[0]}"

rule map:
    input: get_well_split_fq_from_sample
    output: "mapping/{sample}/{sample}/mrfast_out/{sample}.sam.gz"
    params: sge_opts="-N map_{sample} -l mfree=4G -l h_rt=1:0:0:0", output_prefix="mapping/{sample}/{sample}/mrfast_out/{sample}"
    shell:
        "zcat {input} | {MRSFAST_BINARY} --search {MRSFAST_INDEX} {MRSFAST_OPTS} --seq /dev/stdin -o {output} --disable-nohit"

rule split_fastq:
    output: ["mapping/{well}/{well}/fastq_split/{well}.%d_part0.fastq.gz" % num for num in [1, 2]] 
    params: sge_opts="-N split_{well} -q eichler-short.q -l h_rt=6:00:00 -pe orte 2", 
            input_dir="%s/split_barcodes/{well}/{well}/fastq" % SNAKEMAKE_DIR,
            output_dir="%s/mapping/{well}/{well}/fastq_split" % SNAKEMAKE_DIR
    shell:
        "~jlhudd/pipelines/read_depth/scripts/do_splitup_fastqs.sh -i {params.input_dir} -o {params.output_dir}"

rule demultiplex_fastq:
    input: ["fastq/%s%d.fq.gz" % (BARCODE_PREFIX, num) for num in [1,2,3]]
    output: "split_barcodes/{well}/{well}/fastq/{well}.1.fastq.gz", "split_barcodes/{well}/{well}/fastq/{well}.2.fastq.gz"
    params:
    shell:
        "/net/eichler/vol7/home/psudmant/EEE_Lab/projects/fastq_utils/code/parse_barcodes_new/parse_barcodes.py "
        "--barcodes_fa {BARCODE_FASTA} --input_fastq {input[0]} --input_fastq_pe2 {input[2]} "
        "--barcode_fastq {input[1]} --output_dir split_barcodes/"
