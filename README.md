# Clone mapping

Pipeline for generating SUNK pileups from Nextera-based clone sequence data.

## Troubleshooting

Most jobs are killed by our cluster management either for using too much memory
or running too long. All of the rules have log files in the log/ directory with
names like: "<rule>_<sample>.e<jobid>". To figure out why a job died, first find
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
