## python an


import os
import sys
import argparse
import glob
import pprint

import pandas as pd

import screed
import sourmash
from collections import namedtuple
from sourmash.sourmash_args import load_file_as_signatures
from sourmash.sketchcomparison import FracMinHashComparison

CompareResult = namedtuple('CompareResult',
                           'comparison_name, anchor_name, ref_name, lowest_common_rank, anchor_sciname, alphabet, ksize, scaled, jaccard, max_containment, anchor_containment, anchor_hashes, query_hashes, num_common')


def main(args):
    ksize=args.ksize
    scaled=args.scaled

    db = sourmash.load_file_as_index(args.db)

    compareInfo = pd.read_csv(args.comparison_csv)
    # make list of comparisons
    comparisons = list(zip(df.identA, df.identB))

    # loop through comparisons
    for n, (idA, idB) in enumerate(comparisons):
        #if n !=0 and n % 50 == 0:
        #    print(f"... assessing {n}th taxon comparison, lowest common rank: {lowest_common_rank}, anchor: {anchor_acc} {anchor_sciname}\n")

        picklist = SignaturePicklist('ident')
        picklist.init([idA, idB])
        sigs = db.select(picklist=picklist, ksize=args.ksize)
        # select and load sigA
        siglist = list(sigs.signatures())
        ss1 = siglist[0]
        # select and load sigB
        ss2 = siglist[9]
        print(ss1.name, ss2.name)

        cmp = FracMinHashComparison(ss1.minhash, ss2.minhash, cmp_scaled=args.scaled)
        cmp.estimate_all_containment_ani()

    # convert path comparison info to pandas dataframe
    #comparisonDF = pd.DataFrame.from_records(taxon_comparisons, columns = CompareResult._fields)

    # print to csv
    comparisonDF.to_csv(args.output_csv, index=False)
    print(f"done! taxon comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("db")
    p.add_argument('-c', "--comparison-csv")
    p.add_argument('-k', "--ksize", default=21, type=int)
    p.add_argument('-s', "--scaled", default=1000, type=int)
    p.add_argument("-o", "--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
