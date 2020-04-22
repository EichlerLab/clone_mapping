from __future__ import print_function
from __future__ import division

import pandas as pd
import pybedtools
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="bed file with total sunk hits per clone")
    parser.add_argument("sunks", help="bed file with genome-wide sunk positions")
    parser.add_argument("--cores", default=None, help="Tab-delimited table with clone names and core hits")
    parser.add_argument("--read_counts", default=None, help="Tab-delimited file with clone and reads columns")
    parser.add_argument("output")

    args = parser.parse_args()

    dat = pd.read_csv(args.input, header=None, names=["chr", "start", "end", "sunk_hits", "name"], sep='\t')
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

    if args.read_counts is not None:
        rc = pd.read_csv(args.read_counts, sep='\t')
        rc.index = rc.clone
        merged = merged.merge(rc, left_index=True, right_index=True)
        merged = merged[[col for col in merged.columns if col not in ["clone"]]]

    if args.cores is not None:
        cores = pd.read_csv(args.cores, header=None, names=["name", "core_hits"], sep='\t', index_col=0)
        merged = merged.merge(cores, how="left", left_index=True, right_index=True)
    
    merged['length'] = merged['end'] - merged['start']
    merged[['chr','start','end','length','reads','sunk_bases','sunk_depth','sunk_hits','core_hits']].to_csv(args.output, sep="\t", index=True, index_label='clone')


