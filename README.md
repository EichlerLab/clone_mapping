# Clone mapping

Pipeline for generating SUNK pileups from Nextera-based clone sequence data.

## Quick start

### Manifest generation
There is a helper script is included in this repository to help with proper formation of the manifest.
Usage information can be found with 

```bash
./sunk_manifest.sh -h.  
```

Add your samples to the manifest (```manifest.tab```) and check that the settings in ```config.yaml``` make sense for your analysis, then run:
```bash
export DRMAA_LIBRARY_PATH=/opt/uge/lib/lx-amd64/libdrmaa.so.1.0
snakesub -j 20 -w 60 -kT
```
where snakesub is aliased to:
```bash
snakemake --drmaa " -V -cwd -w n -e ./log -o ./log {params.sge_opts} -S /bin/bash"
```

Additionally, if you want to make these changes to your environment, you can edit your ~/.bash_profile to include the following lines
```bash
alias snakesub='snakemake --drmaa " -V -cwd -w n -e ./log -o ./log {params.sge_opts} -S /bin/bash"'
export DRMAA_LIBRARY_PATH=/opt/uge/lib/lx-amd64/libdrmaa.so.1.0
```

To implement changes, type
```bash
source ~/.bash_profile
```

After that, you can just run the following command each time. 
```bash
snakesub -j 20 -w 60 -kT
```

## Full stats

To run the pipeline with additional statistics (SUNK coverage and core hits), add ```get_mapping_stats``` to the snakemake command: 
```bash
snakesub -j 20 -w 60 -kT get_mapping_stats
```
Make sure the reference and reference files settings in ```config.yaml``` make sense for your analysis.

## Troubleshooting

Most jobs are killed by our cluster management either for using too much memory
or running too long. All of the rules have log files in the log/ directory with
names like ```[rule]_[sample].e[jobid]```. To figure out why a job died, first find
the log for the job by the rule and sample name that failed. For example, if a
"map" rule failed for the sample named "BC1_H2", the error log will be in the
file named "log/map_81_BC1_H2.e19174994" where the trailing number is the job id
on the cluster.

View the log file with less or cat. If you see the word "Killed" in the log, you
can check the job's cluster usage statistics with qacct.

```bash
qacct -j 19174994
```

Check the values for the maximum memory used (named "ru_maxrss" and measured in
kilobytes) and job run time (named "ru_wallclock" and measured in seconds). If,
for instance, the max memory used for this job was "4155396" (~4 GB) and the map
job requests exactly 4 GB of memory, we can infer that the job was killed for
using too much memory.
