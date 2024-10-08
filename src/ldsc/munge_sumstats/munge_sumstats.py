#!/usr/bin/env python

from typing import Any
import pandas as pd
import numpy as np
from .check_munge_args import check_munge_args

import sys

from scipy.stats import chi2
from xopen import xopen
from ldsc.ldscore import sumstats

# from ldsc import MASTHEAD, Logger, sec_to_str
import time
from ldsc.logger import LDSCLogger
import logging


logger: logging.Logger = LDSCLogger.get_logger(__name__)


null_values = {"LOG_ODDS": 0, "BETA": 0, "OR": 1, "Z": 0}

default_cnames = {
    # RS NUMBER
    "SNP": "SNP",
    "MARKERNAME": "SNP",
    "SNPID": "SNP",
    "RS": "SNP",
    "RSID": "SNP",
    "RS_NUMBER": "SNP",
    "RS_NUMBERS": "SNP",
    # NUMBER OF STUDIES
    "NSTUDY": "NSTUDY",
    "N_STUDY": "NSTUDY",
    "NSTUDIES": "NSTUDY",
    "N_STUDIES": "NSTUDY",
    # P-VALUE
    "P": "P",
    "PVALUE": "P",
    "P_VALUE": "P",
    "PVAL": "P",
    "P_VAL": "P",
    "GC_PVALUE": "P",
    # ALLELE 1
    "A1": "A1",
    "ALLELE1": "A1",
    "ALLELE_1": "A1",
    "EFFECT_ALLELE": "A1",
    "REFERENCE_ALLELE": "A1",
    "INC_ALLELE": "A1",
    "EA": "A1",
    # ALLELE 2
    "A2": "A2",
    "ALLELE2": "A2",
    "ALLELE_2": "A2",
    "OTHER_ALLELE": "A2",
    "NON_EFFECT_ALLELE": "A2",
    "DEC_ALLELE": "A2",
    "NEA": "A2",
    # N
    "N": "N",
    "NCASE": "N_CAS",
    "CASES_N": "N_CAS",
    "N_CASE": "N_CAS",
    "N_CASES": "N_CAS",
    "N_CONTROLS": "N_CON",
    "N_CAS": "N_CAS",
    "N_CON": "N_CON",
    "N_CASE": "N_CAS",
    "NCONTROL": "N_CON",
    "CONTROLS_N": "N_CON",
    "N_CONTROL": "N_CON",
    "WEIGHT": "N",  # metal does this. possibly risky.
    # SIGNED STATISTICS
    "ZSCORE": "Z",
    "Z-SCORE": "Z",
    "GC_ZSCORE": "Z",
    "Z": "Z",
    "OR": "OR",
    "B": "BETA",
    "BETA": "BETA",
    "LOG_ODDS": "LOG_ODDS",
    "EFFECTS": "BETA",
    "EFFECT": "BETA",
    "SIGNED_SUMSTAT": "SIGNED_SUMSTAT",
    # INFO
    "INFO": "INFO",
    # MAF
    "EAF": "FRQ",
    "FRQ": "FRQ",
    "MAF": "FRQ",
    "FRQ_U": "FRQ",
    "F_U": "FRQ",
}

describe_cname = {
    "SNP": "Variant ID (e.g., rs number)",
    "P": "p-Value",
    "A1": "Allele 1, interpreted as ref allele for signed sumstat.",
    "A2": "Allele 2, interpreted as non-ref allele for signed sumstat.",
    "N": "Sample size",
    "N_CAS": "Number of cases",
    "N_CON": "Number of controls",
    "Z": "Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)",
    "OR": "Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)",
    "BETA": "[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)",
    "LOG_ODDS": "Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)",
    "INFO": "INFO score (imputation quality; higher --> better imputation)",
    "FRQ": "Allele frequency",
    "SIGNED_SUMSTAT": "Directional summary statistic as specified by --signed-sumstats.",
    "NSTUDY": "Number of studies in which the SNP was genotyped.",
}

numeric_cols = [
    "P",
    "N",
    "N_CAS",
    "N_CON",
    "Z",
    "OR",
    "BETA",
    "LOG_ODDS",
    "INFO",
    "FRQ",
    "SIGNED_SUMSTAT",
    "NSTUDY",
]


def read_header(fh):
    """Read the first line of a file and returns a list with the column names."""

    with xopen(fh, "r") as opened_file:

        return [x.rstrip("\n") for x in opened_file.readline().split()]


def get_cname_map(flag, default, ignore):
    """
    Figure out which column names to use.

    Priority is
    (1) ignore everything in ignore
    (2) use everything in flags that is not in ignore
    (3) use everything in default that is not in ignore or in flags

    The keys of flag are cleaned. The entries of ignore are not cleaned. The keys of defualt
    are cleaned. But all equality is modulo clean_header().

    """
    clean_ignore = [clean_header(x) for x in ignore]
    cname_map = {x: flag[x] for x in flag if x not in clean_ignore}
    cname_map.update(
        {x: default[x] for x in default if x not in clean_ignore + list(flag.keys())}
    )
    return cname_map


def clean_header(header):
    """
    For cleaning file headers.
    - convert to uppercase
    - replace dashes '-' with underscores '_'
    - replace dots '.' (as in R) with underscores '_'
    - remove newlines ('\n')
    """
    return header.upper().replace("-", "_").replace(".", "_").replace("\n", "")


def filter_pvals(P):
    """Remove out-of-bounds P-values"""
    ii = (P > 0) & (P <= 1)
    bad_p = (~ii).sum()
    if bad_p > 0:
        logger.warning(
            f"WARNING: {bad_p} SNPs had P outside of (0,1]. The P column may be mislabeled."
        )

    return ii


def filter_info(info, args):
    """Remove INFO < args.info_min (default 0.9) and complain about out-of-bounds INFO."""
    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= args.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = ((info > 2.0) & info.notnull()).any(axis=1) | (
            (info < 0) & info.notnull()
        ).any(axis=1)
        ii = info.sum(axis=1) >= args.info_min * (len(info.columns))
    else:
        raise ValueError("Expected pd.DataFrame or pd.Series.")

    bad_info = jj.sum()
    if bad_info > 0:

        logger.info(
            f"WARNING: {bad_info} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled."
        )

    return ii


def filter_frq(frq, args):
    """
    Filter on MAF. Remove MAF < args.maf_min and out-of-bounds MAF.
    """
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    if bad_frq > 0:
        logger.warning(
            f"WARNING: {bad_frq} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled."
        )

    frq = np.minimum(frq, 1 - frq)
    ii = frq > args.maf_min
    return ii & ~jj


def filter_alleles(a):
    """Remove alleles that do not describe strand-unambiguous SNPs"""
    return a.isin(sumstats.VALID_SNPS)


def parse_dat(dat_gen, convert_colname, merge_alleles, args):
    """Parse and filter a sumstats file chunk-wise"""
    tot_snps = 0
    dat_list = []

    logger.info(
        f"Reading sumstats from {args.sumstats} into memory {args.chunksize} SNPs at a time."
    )
    drops = {"NA": 0, "P": 0, "INFO": 0, "FRQ": 0, "A": 0, "SNP": 0, "MERGE": 0}
    for chunk_number, dat in enumerate(dat_gen, start=1):
        logger.info(f"reading in chunk #{chunk_number}")
        tot_snps += len(dat)
        old = len(dat)
        dat = dat.dropna(
            axis=0, how="any", subset=[x for x in dat.columns if x != "INFO"]
        ).reset_index(drop=True)
        drops["NA"] += old - len(dat)
        dat.columns = [convert_colname[x] for x in dat.columns]

        wrong_types = [
            c
            for c in dat.columns
            if c in numeric_cols and not np.issubdtype(dat[c].dtype, np.number)
        ]
        if len(wrong_types) > 0:
            raise ValueError(f"Columns {wrong_types} are expected to be numeric")

        ii = np.array([True for i in range(len(dat))])
        if args.merge_alleles:
            old = ii.sum()
            ii = dat.SNP.isin(merge_alleles.SNP)
            drops["MERGE"] += old - ii.sum()
            if ii.sum() == 0:
                continue

            dat = dat[ii].reset_index(drop=True)
            ii = np.array([True for i in range(len(dat))])

        if "INFO" in dat.columns:
            old = ii.sum()
            ii &= filter_info(dat["INFO"], args)
            new = ii.sum()
            drops["INFO"] += old - new
            old = new

        if "FRQ" in dat.columns:
            old = ii.sum()
            ii &= filter_frq(dat["FRQ"], args)
            new = ii.sum()
            drops["FRQ"] += old - new
            old = new

        old = ii.sum()
        if args.keep_maf:
            dat.drop([x for x in ["INFO"] if x in dat.columns], inplace=True, axis=1)
        else:
            dat.drop(
                [x for x in ["INFO", "FRQ"] if x in dat.columns], inplace=True, axis=1
            )
        ii &= filter_pvals(dat.P)
        new = ii.sum()
        drops["P"] += old - new
        old = new
        if not args.no_alleles:
            dat.A1 = dat.A1.str.upper()
            dat.A2 = dat.A2.str.upper()
            ii &= filter_alleles(dat.A1 + dat.A2)
            new = ii.sum()
            drops["A"] += old - new
            old = new

        if ii.sum() == 0:
            continue

        dat_list.append(dat[ii].reset_index(drop=True))

    # logger.info(" done\n")
    dat = pd.concat(dat_list, axis=0).reset_index(drop=True)
    logger.info(f"Read {tot_snps} SNPs from --sumstats file.")
    if args.merge_alleles:
        logger.info(f"Removed {drops['MERGE']} SNPs not in --merge-alleles.")

    logger.info(f"Removed {drops['NA']} SNPs with missing values.")
    logger.info(f"Removed {drops['INFO']} SNPs with INFO <= {args.info_min}.")
    logger.info(f"Removed {drops['FRQ']} SNPs with MAF <= {args.maf_min}.")
    logger.info(f"Removed {drops['P']} SNPs with out-of-bounds p-values.")
    logger.info(
        f"Removed {drops['A']} variants that were not SNPs or were strand-ambiguous."
    )
    logger.info(f"{len(dat)} SNPs remain.")
    return dat


def process_n(dat, args):
    """Determine sample size from --N* flags or N* columns. Filter out low N SNPs.s"""
    if all(i in dat.columns for i in ["N_CAS", "N_CON"]):
        N = dat.N_CAS + dat.N_CON
        P = dat.N_CAS / N
        dat["N"] = N * P / P[N == N.max()].mean()
        dat.drop(["N_CAS", "N_CON"], inplace=True, axis=1)
        # NB no filtering on N done here -- that is done in the next code block

    if "N" in dat.columns:
        n_min = args.n_min if args.n_min else dat.N.quantile(0.9) / 1.5
        old = len(dat)
        dat = dat[dat.N >= n_min].reset_index(drop=True)
        new = len(dat)
        logger.info(f"Removed {old-new} SNPs with N < {n_min} ({new} SNPs remain).")

    elif "NSTUDY" in dat.columns and "N" not in dat.columns:
        nstudy_min = args.nstudy_min if args.nstudy_min else dat.NSTUDY.max()
        old = len(dat)
        dat = (
            dat[dat.NSTUDY >= nstudy_min]
            .drop(["NSTUDY"], axis=1)
            .reset_index(drop=True)
        )
        new = len(dat)
        logger.info(
            f"Removed {old - new} SNPs with NSTUDY < {nstudy_min} ({new} SNPs remain)."
        )

    if "N" not in dat.columns:
        if args.N:
            dat["N"] = args.N
            logger.info(f"Using N = {args.N}")
        elif args.N_cas and args.N_con:
            dat["N"] = args.N_cas + args.N_con
            if args.daner is None:
                logger.info(f"Using N_cas = {args.N_cas}; N_con = {args.N_con}")
        else:
            raise ValueError(
                "Cannot determine N. This message indicates a bug.\n"
                "N should have been checked earlier in the program."
            )

    return dat


def p_to_z(P, N):
    """Convert P-value and N to standardized beta."""
    return np.sqrt(chi2.isf(P, 1))


def check_median(x, expected_median, tolerance, name):
    """Check that median(x) is within tolerance of expected_median."""
    m = np.median(x)
    if np.abs(m - expected_median) > tolerance:
        msg = "WARNING: median value of {F} is {V} (should be close to {M}). This column may be mislabeled."
        raise ValueError(msg.format(F=name, M=expected_median, V=round(m, 2)))
    else:
        msg = "Median value of {F} was {C}, which seems sensible.".format(C=m, F=name)

    return msg


def parse_flag_cnames(args):
    """
    Parse flags that specify how to interpret nonstandard column names.

    flag_cnames is a dict that maps (cleaned) arguments to internal column names
    """
    cname_options = [
        [args.nstudy, "NSTUDY", "--nstudy"],
        [args.snp, "SNP", "--snp"],
        [args.N_col, "N", "--N"],
        [args.N_cas_col, "N_CAS", "--N-cas-col"],
        [args.N_con_col, "N_CON", "--N-con-col"],
        [args.a1, "A1", "--a1"],
        [args.a2, "A2", "--a2"],
        [args.p, "P", "--P"],
        [args.frq, "FRQ", "--nstudy"],
        [args.info, "INFO", "--info"],
    ]
    flag_cnames = {clean_header(x[0]): x[1] for x in cname_options if x[0] is not None}
    if args.info_list:
        try:
            flag_cnames.update(
                {clean_header(x): "INFO" for x in args.info_list.split(",")}
            )
        except ValueError:
            logger.critical(
                "The argument to --info-list should be a comma-separated list of column names."
            )
            raise

    null_value = None
    if args.signed_sumstats:
        try:
            cname, null_value = args.signed_sumstats.split(",")
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = "SIGNED_SUMSTAT"
        except ValueError:
            logger.critical(
                "The argument to --signed-sumstats should be column header comma number."
            )
            raise ValueError()

    return [flag_cnames, null_value]


def allele_merge(dat, alleles):
    """
    WARNING: dat now contains a bunch of NA's~
    Note: dat now has the same SNPs in the same order as --merge alleles.
    """
    dat = pd.merge(alleles, dat, how="left", on="SNP", sort=False).reset_index(
        drop=True
    )
    ii = dat.A1.notnull()
    a1234 = dat.A1[ii] + dat.A2[ii] + dat.MA[ii]
    match = a1234.apply(lambda y: y in sumstats.MATCH_ALLELES)
    jj = pd.Series(np.zeros(len(dat), dtype=bool))
    jj.loc[ii] = match
    old = ii.sum()
    n_mismatch = (~match).sum()
    if n_mismatch < old:
        logger.info(
            "Removed {M} SNPs whose alleles did not match --merge-alleles ({N} SNPs remain).".format(
                M=n_mismatch, N=old - n_mismatch
            )
        )
    else:
        logger.critical("All SNPs have alleles that do not match --merge-alleles.")
        raise ValueError("All SNPs have alleles that do not match --merge-alleles.")

    dat.loc[~jj.astype("bool"), [i for i in dat.columns if i != "SNP"]] = float("nan")
    dat.drop(["MA"], axis=1, inplace=True)
    return dat


# set p = False for testing in order to prevent printing
def munge_sumstats(args):

    check_munge_args(args)
    file_cnames = read_header(args.sumstats)  # note keys not cleaned
    flag_cnames, signed_sumstat_null = parse_flag_cnames(args)
    if args.ignore:
        ignore_cnames = [clean_header(x) for x in args.ignore.split(",")]
    else:
        ignore_cnames = []

    # remove LOG_ODDS, BETA, Z, OR from the default list
    if args.signed_sumstats is not None or args.a1_inc:
        mod_default_cnames = {
            x: default_cnames[x]
            for x in default_cnames
            if default_cnames[x] not in null_values
        }
    else:
        mod_default_cnames = default_cnames

    cname_map = get_cname_map(flag_cnames, mod_default_cnames, ignore_cnames)
    if args.daner:
        frq_u = [x for x in file_cnames if x.startswith("FRQ_U_")][0]
        frq_a = [x for x in file_cnames if x.startswith("FRQ_A_")][0]
        N_cas = float(frq_a[6:])
        N_con = float(frq_u[6:])
        logger.info(
            "Inferred that N_cas = {N1}, N_con = {N2} from the FRQ_[A/U] columns.".format(
                N1=N_cas, N2=N_con
            )
        )
        args.N_cas = N_cas
        args.N_con = N_con
        # drop any N, N_cas, N_con or FRQ columns
        for c in ["N", "N_CAS", "N_CON", "FRQ"]:
            for d in [x for x in cname_map if cname_map[x] == "c"]:
                del cname_map[d]

        cname_map[frq_u] = "FRQ"

    if args.daner_n:
        frq_u = [x for x in file_cnames if x.startswith("FRQ_U_")][0]
        cname_map[frq_u] = "FRQ"
        try:
            dan_cas = clean_header(file_cnames[file_cnames.index("Nca")])
        except ValueError:
            raise ValueError("Could not find Nca column expected for daner-n format")

        try:
            dan_con = clean_header(file_cnames[file_cnames.index("Nco")])
        except ValueError:
            raise ValueError("Could not find Nco column expected for daner-n format")

        cname_map[dan_cas] = "N_CAS"
        cname_map[dan_con] = "N_CON"

    cname_translation = {
        x: cname_map[clean_header(x)]
        for x in file_cnames
        if clean_header(x) in cname_map
    }  # note keys not cleaned
    cname_description = {
        x: describe_cname[cname_translation[x]] for x in cname_translation
    }
    if args.signed_sumstats is None and not args.a1_inc:
        sign_cnames = [
            x for x in cname_translation if cname_translation[x] in null_values
        ]
        if len(sign_cnames) > 1:
            raise ValueError(
                "Too many signed sumstat columns. Specify which to ignore with the --ignore flag."
            )
        if len(sign_cnames) == 0:
            raise ValueError("Could not find a signed summary statistic column.")

        sign_cname = sign_cnames[0]
        signed_sumstat_null = null_values[cname_translation[sign_cname]]
        cname_translation[sign_cname] = "SIGNED_SUMSTAT"
    else:
        sign_cname = "SIGNED_SUMSTATS"

    # check that we have all the columns we need
    if not args.a1_inc:
        req_cols = ["SNP", "P", "SIGNED_SUMSTAT"]
    else:
        req_cols = ["SNP", "P"]

    for c in req_cols:
        if c not in list(cname_translation.values()):
            raise ValueError("Could not find {C} column.".format(C=c))

    # check aren't any duplicated column names in mapping
    for field in cname_translation:
        numk = file_cnames.count(field)
        if numk > 1:
            raise ValueError(
                "Found {num} columns named {C}".format(C=field, num=str(numk))
            )

    # check multiple different column names don't map to same data field
    for head in list(cname_translation.values()):
        numc = list(cname_translation.values()).count(head)
        if numc > 1:
            raise ValueError(
                "Found {num} different {C} columns".format(C=head, num=str(numc))
            )

    if (
        (not args.N)
        and (not (args.N_cas and args.N_con))
        and ("N" not in list(cname_translation.values()))
        and (any(x not in list(cname_translation.values()) for x in ["N_CAS", "N_CON"]))
    ):
        raise ValueError("Could not determine N.")
    if (
        "N" in list(cname_translation.values())
        or all(x in list(cname_translation.values()) for x in ["N_CAS", "N_CON"])
    ) and "NSTUDY" in list(cname_translation.values()):
        nstudy = [x for x in cname_translation if cname_translation[x] == "NSTUDY"]
        for x in nstudy:
            del cname_translation[x]
    if not args.no_alleles and not all(
        x in list(cname_translation.values()) for x in ["A1", "A2"]
    ):
        raise ValueError("Could not find A1/A2 columns.")

    logger.info("Interpreting column names as follows:")
    logger.info(
        "\n".join([x + ":\t" + cname_description[x] for x in cname_description]) + "\n"
    )

    if args.merge_alleles:
        logger.info(
            "Reading list of SNPs for allele merge from {F}".format(
                F=args.merge_alleles
            )
        )

        merge_alleles = pd.read_csv(
            args.merge_alleles,
            header=0,
            sep="\s+",
            na_values=".",
        )
        if any(x not in merge_alleles.columns for x in ["SNP", "A1", "A2"]):
            raise ValueError("--merge-alleles must have columns SNP, A1, A2.")

        logger.info("Read {N} SNPs for allele merge.".format(N=len(merge_alleles)))
        merge_alleles["MA"] = (merge_alleles.A1 + merge_alleles.A2).apply(
            lambda y: y.upper()
        )
        merge_alleles.drop(
            [x for x in merge_alleles.columns if x not in ["SNP", "MA"]],
            axis=1,
            inplace=True,
        )
    else:
        merge_alleles = None

    # figure out which columns are going to involve sign information, so we can ensure
    # they're read as floats
    signed_sumstat_cols = [
        k for k, v in list(cname_translation.items()) if v == "SIGNED_SUMSTAT"
    ]
    dat_gen = pd.read_csv(
        args.sumstats,
        sep="\s+",
        header=0,
        usecols=list(cname_translation.keys()),
        na_values=[".", "NA"],
        iterator=True,
        chunksize=args.chunksize,
        dtype={c: np.float64 for c in signed_sumstat_cols},
    )

    dat = parse_dat(dat_gen, cname_translation, merge_alleles, args)
    if len(dat) == 0:
        raise ValueError("After applying filters, no SNPs remain.")

    old = len(dat)
    dat = dat.drop_duplicates(subset="SNP").reset_index(drop=True)
    new = len(dat)
    logger.info(
        "Removed {M} SNPs with duplicated rs numbers ({N} SNPs remain).".format(
            M=old - new, N=new
        )
    )
    # filtering on N cannot be done chunkwise
    dat = process_n(dat, args)
    dat.P = p_to_z(dat.P, dat.N)
    dat.rename(columns={"P": "Z"}, inplace=True)
    if not args.a1_inc:
        logger.info(
            check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, 0.1, sign_cname)
        )
        dat.Z *= (-1) ** (dat.SIGNED_SUMSTAT < signed_sumstat_null)
        dat.drop("SIGNED_SUMSTAT", inplace=True, axis=1)
    # do this last so we don't have to worry about NA values in the rest of
    # the program
    if args.merge_alleles:
        dat = allele_merge(dat, merge_alleles)

    out_fname = args.out.parent / f"{args.out.name}.sumstats.gz"

    print_colnames = [c for c in dat.columns if c in ["SNP", "N", "Z", "A1", "A2"]]
    if args.keep_maf and "FRQ" in dat.columns:
        print_colnames.append("FRQ")

    logger.info(
        f"Writing summary statistics for {len(dat)} SNPs ({dat.N.notnull().sum()} with nonmissing beta) to {out_fname}."
    )

    dat.to_csv(
        out_fname,
        sep="\t",
        index=False,
        columns=print_colnames,
        float_format="%.3f",
        compression="gzip",
    )

    logger.info("\nMetadata:")
    CHISQ = dat.Z**2
    mean_chisq = CHISQ.mean()
    logger.info("Mean chi^2 = " + str(round(mean_chisq, 3)))
    if mean_chisq < 1.02:
        logger.info("WARNING: mean chi^2 may be too small.")

    logger.info("Lambda GC = " + str(round(CHISQ.median() / 0.4549, 3)))
    logger.info("Max chi^2 = " + str(round(CHISQ.max(), 3)))
    logger.info(
        "{N} Genome-wide significant SNPs (some may have been removed by filtering).".format(
            N=(CHISQ > 29).sum()
        )
    )

    logger.info("\nConversion finished at {T}".format(T=time.ctime()))
    return dat

    # except Exception:
    #     log.log("\nERROR converting summary statistics:\n")
    #     log.log(traceback.format_exc())
    #     raise
    # finally:
    #     log.log("\nConversion finished at {T}".format(T=time.ctime()))
