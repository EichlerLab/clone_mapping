import os
import pandas as pd
import tempfile

TMPDIR = tempfile.gettempdir()

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/config.sh; set -eo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "%s/config.yaml" % SNAKEMAKE_DIR

REFERENCE = config["reference"]

BARCODE_FASTA = config["barcode_fasta"]
BARCODE_PREFIX = config["barcode_prefix"]

TRACK_OUTPUT_DIR = config["track_output_dir"]
TRACK_URL = config["track_url"]

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

def get_sample_and_num(wildcards):
    reads = MANIFEST.loc[wildcards.sample, "reads"].split(",")
    return ["read_counts/raw/{sample}.{num}.txt".format(sample=wildcards.sample, num=num) for num in range(1, len(reads)+1)]

def get_well_split_fq_from_sample(wildcards):
    reads = MANIFEST.loc[wildcards.sample, "reads"].split(",")
    return ["mapping/{sample}/{sample}/fastq_split/{sample}.{num}_part0.fastq.gz".format(sample=wildcards.sample, num=num) for num in range(1, len(reads)+1)]

localrules: all, make_tracks, make_bams, clean

rule all:
    input: "clone_locations.bed", "clone_mapping_tracklist.txt"

rule clean:
    shell:
        "rm bam/* mapping/*/*/mrsfast_out/* mapping/*/*/fastq_split/*"

rule get_mapping_stats:
    input: "clone_locations.bed", "bam/bamlist.txt", "read_counts/combined.tab"
    output: "clone_locations.annotated.tab"
    params: sge_opts = "-l mfree=4G -l h_rt=01:00:00", sunks=config[REFERENCE]["sunk_bed"], cores=config[REFERENCE]["cores"]
    shell:
        "python get_clone_mapping_stats.py {input[0]} {params.sunks} {output} --cores {params.cores} --bamlist {input[1]} --read_counts {input[2]}"

rule collect_pileup_locations:
    input: expand("clone_locations/{sample}.bed", sample=MANIFEST.sample_name)
    output: "clone_locations.bed"
    params: sge_opts = "-l mfree=1G -l h_rt=00:10:00"
    benchmark: "benchmarks/clone_loc.txt"
    shell: "sort -k 5,5V {input} > {output}"

rule get_pileup_locations:
    input: "sunk_pileup/{sample}.sorted.bam_sunk.bw"
    output: "clone_locations/{sample}.bed"
    params: sge_opts = "-l mfree=1G -l h_rt=00:10:00"
    benchmark: "benchmarks/clone_locations/{sample}.txt"
    shell: """bigWigToBedGraph {input} /dev/stdout \
            | awk 'OFS="\\t" {{ print $1,$2,$3,".",$4 }}' \
            | bedtools merge -i stdin -d 100000 -c 5 -o sum | sort -k 4,4rn \
            | awk 'OFS="\\t" {{ if (NR == 1) {{ print $1,$2,$3,$4,"{wildcards.sample}" }} }}' > {output}"""

rule make_tracks:
    input: expand("sunk_pileup/{sample}.sorted.bam_sunk.bw.trackdef", sample=MANIFEST.sample_name)
    output: "clone_mapping_tracklist.txt"
    shell:
        "sort -V {input} > {output}; "
        "mkdir -p {TRACK_OUTPUT_DIR}; "
        "chmod 755 {TRACK_OUTPUT_DIR}; "
        "rsync -arv --bwlimit=70000 sunk_pileup/*.bw {TRACK_OUTPUT_DIR}; "
        "rsync {output[0]} {TRACK_OUTPUT_DIR}; "
        "chmod 644 {TRACK_OUTPUT_DIR}/*"

rule make_bw_pileup:
    input: "bam/{sample}.sorted.bam"
    output: "sunk_pileup/{sample}.sorted.bam_sunk.bw", "sunk_pileup/{sample}.sorted.bam_sunk.bw.trackdef"
    params: sge_opts="-N plp_{sample} -l h_rt=0:20:00 -l mfree=2G"
    benchmark: "benchmarks/pileup/{sample}.txt"
    run:
        if "EICHLERLAB_PASSWORD" not in os.environ:
            print("Error: Shell variable 'EICHLERLAB_PASSWORD' not set. Exitting...")
            sys.exit(1)
        PASSWORD = os.environ["EICHLERLAB_PASSWORD"]
        shell("python /net/eichler/vol2/local/inhousebin/sunks/pileups/sam_to_bw_pileup.py "
        "--inSam {SNAKEMAKE_DIR}/{input} --contigs {CONTIGS} "
        "--outdir {SNAKEMAKE_DIR}/sunk_pileup "
        "--sunk_mask {SUNK_MASK} --track_url https://eee:{PASSWORD}@{TRACK_URL}")

rule make_bams:
    input: expand("bam/{sample}.sorted.{ext}", sample=MANIFEST.sample_name, ext=["bam", "bam.bai"])
    output: "bam/bamlist.txt"
    params: sge_opts="-l mfree=1G -l h_rt=1:0:0"
    run:
        with open(output[0], "w") as outfile:
            for clone in MANIFEST.sample_name:
                cwd = os.getcwd()
                fn = "%s/bam/%s.sorted.bam" % (cwd, clone)
                print(clone, fn, sep="\t", file=outfile)
            
rule make_bam:
    input: "mapping/{sample}/{sample}/mrsfast_out/{sample}.sam.gz"
    output: "bam/{sample}.sorted.bam", "bam/{sample}.sorted.bam.bai"
    params: sge_opts="-N bam_{sample} -l mfree=4G -l h_rt=02:00:00",
            output_prefix="bam/{sample}.sorted"
    benchmark: "benchmarks/make_bam/{sample}.txt"
    shell:
        "zcat {input} | samtools view -b -t {CONTIGS} - -S | samtools sort - -T {TMPDIR}/{wildcards.sample} -o {output[0]} -m 2G; "
        "samtools index {output[0]}"

rule map:
    input: get_well_split_fq_from_sample
    output: "mapping/{sample}/{sample}/mrsfast_out/{sample}.sam.gz"
    params: sge_opts="-N map_{sample} -l mfree=6G -l h_rt=04:00:00",
            output_prefix="mapping/{sample}/{sample}/mrsfast_out/{sample}"
    benchmark: "benchmarks/map/{sample}.txt"
    shell:
        "zcat {input} | {MRSFAST_BINARY} --search {MRSFAST_INDEX} {MRSFAST_OPTS} --seq /dev/stdin -o {params.output_prefix} --disable-nohit"

rule merge_read_counts:
    input: expand("read_counts/merged/{sample}.txt", sample=MANIFEST.sample_name)
    output: "read_counts/combined.tab"
    params: sge_opts="-l mfree=1G -l h_rt=1:0:0"
    shell:
        """echo -e "clone\\treads" > {output}
           cat {input} >> {output}"""

rule merge_read_counts_by_sample:
    input: get_sample_and_num
    output: "read_counts/merged/{sample}.txt"
    params: sge_opts="-l mfree=1G -l h_rt=1:0:0"
    run:
        count = 0
        for fn in input:
            with open(fn, "r") as infile:
                sn_num, reads = infile.readlines()[0].rstrip().split()
                count += int(reads)
        with open(output[0], "w") as outfile:
            print(wildcards.sample, count, sep="\t", file=outfile)

rule split_fastq:
    input: lambda wildcards: MANIFEST.loc[wildcards.sample, "reads"].split(",")[int(wildcards.num)-1]
    output: "mapping/{sample}/{sample}/fastq_split/{sample}.{num}_part0.fastq.gz", "read_counts/raw/{sample}.{num}.txt"
    params: sge_opts="-N split_{sample} -l h_rt=6:00:00 -l disk_free=10G -l mfree=1G",
            split_read_length=36
    benchmark: "benchmarks/split_fastq/{sample}.{num}.txt"
    shell:
        "python split_reads.py --full_length_only {input} --read_counts {output[1]} --clone_name {wildcards.sample}.{wildcards.num} {params.split_read_length} | "
        "bgzip -c | pv -L 50000K > {output[0]}"

#rule demultiplex_fastq:
#    input: ["fastq/%s%d.fq.gz" % (BARCODE_PREFIX, num) for num in [1,2,3]]
#    output: "split_barcodes/{sample}/{sample}/fastq/{sample}.1.fastq.gz", "split_barcodes/{sample}/{sample}/fastq/{sample}.2.fastq.gz"
#    params: sge_opts=""
#    shell:
#        "/net/eichler/vol7/home/psudmant/EEE_Lab/projects/fastq_utils/code/parse_barcodes_new/parse_barcodes.py "
#        "--barcodes_fa {BARCODE_FASTA} --input_fastq {input[0]} --input_fastq_pe2 {input[2]} "
#        "--barcode_fastq {input[1]} --output_dir split_barcodes/"
