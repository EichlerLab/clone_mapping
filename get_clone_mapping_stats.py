from __future__ import print_function
from __future__ import division

import pandas as pd
import pybedtools
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="bed file with total sunk hits per clone")
    parser.add_argument("sunks", help="bed file with genome-wide sunk positions")
    parser.add_argument("--cores", default=None, help="Bed file with pseudocore locations as in Dennis et al. 2017")
    parser.add_argument("--bamlist", default=None, help="Tab-delimited file with clone name and and bam columns. Required for --cores")
    parser.add_argument("output")

    args = parser.parse_args()

    if args.cores is not None and args.bamlist is None:
        sys.exit("Bamlist must be specified if --cores is specified.")

    dat = pd.read_table(args.input, header=None, names=["chr", "start", "end", "sunk_hits", "name"])
    dat.index = dat.name

    sunks = pybedtools.BedTool(args.sunks)

    clones = pybedtools.BedTool(dat[["chr", "start", "end", "name"]].to_string(header=None, index=False), from_string=True)

    sunks_per_clone = {name: 0 for name in dat.name}
    for entry in clones.intersect(sunks, wao=True):
        name = entry[3]
        bases = int(entry[7])
        sunks_per_clone[name] += bases

    sunk_dat = pd.DataFrame.from_dict(sunks_per_clone, orient="index")
    sunk_dat.columns = ["sunk_bases"]
    merged = dat.merge(sunk_dat, how="left", left_index=True, right_index=True)
    merged["sunk_depth"] = merged.sunk_hits / merged.sunk_bases

    if args.cores is not None:
        cores = pybedtools.BedTool(args.cores)
        bamfiles = pd.read_table(args.bamlist, header=None, names=["clone", "bamfile"])
        bamfiles.index = bamfiles.clone
        total_cov = {name: 0 for name in merged.name}
        for clone in merged.name:
            bamfile = bamfiles.loc[clone, "bamfile"]
            bambed = pybedtools.BedTool(bamfile)
            cov = bambed.coverage(b=cores)
            total_cov[clone] = sum([int(entry[4]) for entry in cov])
        cov_dat = pd.DataFrame.from_dict(total_cov, orient="index")
        cov_dat.columns = ["core_hits"]
        final = merged.merge(cov_dat, how="left", left_index=True, right_index=True)
    else:
        final = merged
    final.to_csv(args.output, sep="\t", index=False)


