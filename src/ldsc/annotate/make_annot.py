#!/usr/bin/env python

import pandas as pd
import numpy as np
from pybedtools import BedTool
import logging

from ldsc.logger import LDSCLogger

logger: logging.Logger = LDSCLogger.get_logger(__name__)


def gene_set_to_bed(args) -> BedTool:
    print("making gene set bed file")
    GeneSet = pd.read_csv(args.gene_set_file, header=None, names=["GENE"])
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace=True)
    df = pd.merge(GeneSet, all_genes, on="GENE", how="inner")
    df["START"] = np.maximum(1, df["START"] - args.windowsize)
    df["END"] = df["END"] + args.windowsize
    iter_df = [
        ["chr" + (str(x1).lstrip("chr")), x2 - 1, x3]
        for (x1, x2, x3) in np.array(df[["CHR", "START", "END"]])
    ]
    return BedTool(iter_df).sort().merge()


def preprocess_args(args) -> BedTool:
    """check some of the arguments to appropriately generate
    a Bedtool object
    """
    if args.gene_set_file is not None:
        bed_for_annot = gene_set_to_bed(args)
    else:
        bed_for_annot = BedTool(args.bed_file).sort()
        if not args.nomerge:
            bed_for_annot = bed_for_annot.merge()

    return bed_for_annot


def make_annot_files(args) -> None:

    bed_for_annot = preprocess_args(args)
    logger.info("making annot file")
    df_bim = pd.read_csv(
        args.bimfile,
        delim_whitespace=True,
        usecols=[0, 1, 2, 3],
        names=["CHR", "SNP", "CM", "BP"],
    )
    iter_bim = [
        ["chr" + str(x1), x2 - 1, x2] for (x1, x2) in np.array(df_bim[["CHR", "BP"]])
    ]
    bimbed = BedTool(iter_bim)
    annotbed = bimbed.intersect(bed_for_annot)
    bp = [x.start + 1 for x in annotbed]
    df_int = pd.DataFrame({"BP": bp, "ANNOT": 1})
    df_annot = pd.merge(df_bim, df_int, how="left", on="BP")
    df_annot.fillna(0, inplace=True)
    df_annot = df_annot[["ANNOT"]].astype(int)

    df_annot.to_csv(args.annot_file, sep="\t", index=False, compression="infer")
