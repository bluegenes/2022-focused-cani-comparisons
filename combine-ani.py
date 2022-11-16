import os
import sys
import argparse
import pandas as pd
from functools import partial

from sourmash.tax import tax_utils
from sourmash.lca import lca_utils
from sourmash.logging import notify


def find_lca(row):
    idA = row['identA']
    linA = tax_assign[idA]
    idB = row['identB']
    linB = tax_assign[idB]
    lintree = build_tree(linA, linB)
    lca, node_len = lca_utils.find_lca(lintree)
    row['lca_lineage'] = lca_utils.display_lineage(lca)
    row['lca_rank'] = lca_utils.display_lineage(lca[-1].rank)
    return row

def main(args):

    pyani_res = pd.read_csv(args.pyani_csv)
    fastani_res = pd.read_csv(args.fastani_csv)
    sourmash_res = pd.read_csv(args.sourmash)

    ani_dfs = [pyani_res, fastani_res, sourmash_res]

    #id_cols = ["comparison_name", "identA", "identB"]
    outer_merge = partial(pd.merge, how='outer')
    aniDF = reduce(outer_merge, ani_dfs)

    if args.taxonomy:
        # could add LCA info here if we want it.
        tax_assign = tax_utils.MultiLineageDB.load([args.taxonomy],
                                                   keep_full_identifiers=False,
                                                   keep_identifier_versions=False)
        aniDF = aniDF.apply(find_lca, axis=1)

    # print to csv
    aniDF.to_csv(args.output_csv, index=False)
    print(f"done! Combined comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("-p", "--pyani-csv")
    p.add_argument("-f", "--fastani-csv")
    p.add_argument("-s", "--sourmash-csv")
    p.add_argument("-o", "--output-csv", required=True)
    p.add_argument("-t", "--taxonomy", default="gtdb-rs207.taxonomy.csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

