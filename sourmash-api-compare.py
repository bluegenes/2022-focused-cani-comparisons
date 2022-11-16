## python api comparisons

import os
import sys
import argparse
import glob
import pprint

import pandas as pd

import sourmash
from collections import namedtuple
from sourmash.sourmash_args import load_file_as_signatures
from sourmash.sketchcomparison import FracMinHashComparison
from sourmash.picklist import SignaturePicklist
from sourmash.logging import notify

CompareResult = namedtuple('CompareResult', ['comparison_name', 'identA', 'identB',
                                             'ksize', 'scaled', 'avg_cANI', 'jaccard',
                                             'idA_in_idB', 'idA_in_idB_cANI', 'idB_in_idA',
                                             'idB_in_idA_cANI','idA_hashes', 'idB_hashes',
                                             'num_common'])


def main(args):
    ksize=args.ksize
    scaled=args.scaled

    db = sourmash.load_file_as_index(args.db)

    compareInfo = pd.read_csv(args.comparison_csv)
    # make list of comparisons
    comparisons = list(zip(compareInfo.identA, compareInfo.identB))

    # loop through comparisons
    results = []
    for n, (idA, idB) in enumerate(comparisons):
        if n !=0 and n % 50 == 0:
            notify(f"... assessing {n}th comparison\n")

        comparison_name = f"{idA}_x_{idB}"
        picklist = SignaturePicklist('ident')
        picklist.init([idA, idB])
        sigs = db.select(picklist=picklist, ksize=args.ksize)
        # select and load sigA
        siglist = list(sigs.signatures())
        ss1 = siglist[0]
        # select and load sigB
        ss2 = siglist[1]
        print(ss1.name, ss2.name)

        cmp = FracMinHashComparison(ss1.minhash, ss2.minhash, cmp_scaled=args.scaled)
        cmp.estimate_all_containment_ani()
        idA_sc_hashes = len(cmp.mh1_cmp)
        idB_sc_hashes = len(cmp.mh2_cmp)
        contain1 = cmp.mh1_containment_in_mh2
        cANI_1 = cmp.ani_from_mh2_containment_in_mh1
        contain2 = cmp.mh2_containment_in_mh1
        cANI_2 = cmp.ani_from_mh1_containment_in_mh2
        res = CompareResult(comparison_name, idA, idB, ksize, scaled, \
                             cmp.avg_containment_ani, cmp.jaccard, contain1, \
                             cANI_1, contain2, cANI_2, idA_sc_hashes, idB_sc_hashes, \
                             cmp.total_unique_intersect_hashes)
        results.append(res)

    # convert comparison info to pandas dataframe
    comparisonDF = pd.DataFrame.from_records(results, columns = CompareResult._fields)

    # print to csv
    comparisonDF.to_csv(args.output_csv, index=False)
    print(f"done! comparison info written to {args.output_csv}")

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
