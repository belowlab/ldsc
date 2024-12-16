"""
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

This module deals with getting all the data needed for LD Score regression from files
into memory and checking that the input makes sense. There is no math here. LD Score
regression is implemented in the regressions module.
"""

import logging
from pathlib import Path
from typing import Union
import numpy as np
import pandas as pd
from scipy import stats
import itertools as it
from . import parse as ps
from . import regressions as reg
import sys
import traceback
import copy
import functools
from .ldsc_check_args import check_args
from ldsc.logger import LDSCLogger

logger: logging.Logger = LDSCLogger.get_logger(__name__)


_N_CHR = 22
# complementary bases
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
# bases
BASES = list(COMPLEMENT.keys())
# true iff strand ambiguous
STRAND_AMBIGUOUS = {
    "".join(x): x[0] == COMPLEMENT[x[1]]
    for x in it.product(BASES, BASES)
    if x[0] != x[1]
}
# SNPS we want to keep (pairs of alleles)
VALID_SNPS = {
    x
    for x in ["".join(y) for y in it.product(BASES, BASES)]
    if x[0] != x[1] and not STRAND_AMBIGUOUS[x]
}
# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {
    x
    for x in ["".join(y) for y in it.product(VALID_SNPS, VALID_SNPS)]
    # strand and ref match
    if ((x[0] == x[2]) and (x[1] == x[3])) or
    # ref match, strand flip
    ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
    # ref flip, strand match
    ((x[0] == x[3]) and (x[1] == x[2]))
    or ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
}  # strand and ref flip
# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {
    "".join(x): ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
    # strand flip
    ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
    for x in MATCH_ALLELES
}


def _splitp(fstr: Union[str, Path]):
    if isinstance(fstr, str):
        flist = fstr.split(",")
        flist = [Path(x) for x in flist]
    else:
        flist = [fstr]

    return flist


def _select_and_log(x, ii):
    """Fiter down to rows that are True in ii. Log # of SNPs removed."""
    new_len = ii.sum()
    msg = f"{new_len if new_len != 0 else 0} SNPs with valid alleles."

    if new_len == 0:
        raise ValueError(msg)
    else:
        x = x[ii]

        logger.info(msg)
    return x


def smart_merge(x, y):
    """Check if SNP columns are equal. If so, save time by using concat instead of merge."""
    if len(x) == len(y) and (x.index == y.index).all() and (x.SNP == y.SNP).all():
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True).drop("SNP", axis=1)
        out = pd.concat([x, y], axis=1)
    else:
        out = pd.merge(x, y, how="inner", on="SNP")
    return out


def _read_ref_ld(args):
    """Read reference LD Scores."""
    ref_ld = _read_chr_split_files(
        args.ref_ld_chr,
        args.ref_ld,
        "reference panel LD Score",
        ps.ldscore_fromlist,
    )
    logger.info(f"Read reference panel LD Scores for {len(ref_ld)} SNPs.")
    return ref_ld


def _read_annot(
    args,
):
    """Read annot matrix."""
    try:
        if args.ref_ld is not None:
            overlap_matrix, M_tot = _read_chr_split_files(
                args.ref_ld_chr,
                args.ref_ld,
                "annot matrix",
                ps.annot,
                frqfile=args.frqfile,
            )
        elif args.ref_ld_chr is not None:
            overlap_matrix, M_tot = _read_chr_split_files(
                args.ref_ld_chr,
                args.ref_ld,
                "annot matrix",
                ps.annot,
                frqfile=args.frqfile_chr,
            )
    except Exception:
        logger.critical("Error parsing .annot file.")
        raise

    return overlap_matrix, M_tot


def _read_M(args, n_annot):
    """Read M (--M, --M-file, etc)."""
    if args.M:
        try:
            M_annot = [float(x) for x in args.M.split(",")]
        except ValueError as e:
            raise ValueError("Could not cast --M to float: " + str(e.args))
    else:
        if args.ref_ld:
            M_annot = ps.M_fromlist(_splitp(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld_chr), _N_CHR, common=(not args.not_M_5_50)
            )

    try:
        M_annot = np.array(M_annot).reshape((1, n_annot))
    except ValueError as e:
        raise ValueError(
            "# terms in --M must match # of LD Scores in --ref-ld.\n" + str(e.args)
        )

    return M_annot


def _read_w_ld(args):
    """Read regression SNP LD."""
    if (args.w_ld and "," in args.w_ld.name) or (
        args.w_ld_chr and "," in args.w_ld_chr.name
    ):
        raise ValueError("--w-ld must point to a single fileset (no commas allowed).")

    w_ld = _read_chr_split_files(
        args.w_ld_chr, args.w_ld, "regression weight LD Score", ps.ldscore_fromlist
    )
    if len(w_ld.columns) != 2:
        raise ValueError("--w-ld may only have one LD Score column.")
    w_ld.columns = ["SNP", "LD_weights"]  # prevent colname conflicts w/ ref ld
    logger.info(f"Read regression weight LD Scores for {len(w_ld)} SNPs.")
    return w_ld


def _read_chr_split_files(chr_arg, not_chr_arg, noun, parsefunc, **kwargs):
    """Read files split across 22 chromosomes (annot, ref_ld, w_ld)."""

    try:
        if not_chr_arg:
            logger.info(f"Reading {noun} from {not_chr_arg} ... ({parsefunc.__name__})")
            out = parsefunc(_splitp(not_chr_arg), **kwargs)
        elif chr_arg:

            f = ps.sub_chr(chr_arg, "[1-22]")
            logger.info(f"Reading {noun} from {f} ... ({parsefunc.__name__})")

            out = parsefunc(_splitp(chr_arg), _N_CHR, **kwargs)

    except ValueError as e:
        logger.critical(f"Error parsing {noun}.")
        raise e

    return out


def _read_sumstats(args: list[any], fh: Path, alleles=False, dropna=False):
    """Parse summary statistics."""
    logger.info(f"Reading summary statistics from {fh} ...")
    sumstats = ps.sumstats(fh, alleles=alleles, dropna=dropna)
    logger.info(f"Read summary statistics for {len(sumstats)} SNPs.")
    m = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset="SNP")
    if m > len(sumstats):
        logger.info(f"Dropped {m - len(sumstats)} SNPs with duplicated rs numbers.")

    return sumstats


def _check_ld_condnum(args, ref_ld):
    """Check condition number of LD Score matrix."""
    if len(ref_ld.shape) >= 2:
        cond_num = int(np.linalg.cond(ref_ld))
        if cond_num > 100000:
            if args.invert_anyway:
                logger.warning(
                    f"WARNING: LD Score matrix condition number is {cond_num}. Inverting anyway because the --invert-anyway flag is set."
                )
            else:
                logger.warning(
                    f"WARNING: LD Score matrix condition number is {cond_num}. Remove collinear LD Scores. "
                )
                raise ValueError(
                    f"WARNING: LD Score matrix condition number is {cond_num}. Remove collinear LD Scores. "
                )


def _check_variance(M_annot, ref_ld):
    """Remove zero-variance LD Scores."""
    ii = ref_ld.iloc[:, 1:].var() == 0  # NB there is a SNP column here
    if ii.all():
        logger.critical("All LD Scores have zero variance.")
        raise ValueError("All LD Scores have zero variance.")
    else:
        logger.info("Removing partitioned LD Scores with zero variance.")
        ii_snp = np.array([True] + list(~ii))
        ii_m = np.array(~ii)
        ref_ld = ref_ld.iloc[:, ii_snp]
        M_annot = M_annot[:, ii_m]

    return M_annot, ref_ld, ii


def _warn_length(sumstats):
    if len(sumstats) < 200000:
        logger.warning(
            "WARNING: number of SNPs less than 200k; this is almost always bad."
        )


def _print_cov(ldscore_reg, ofh):
    """Prints covariance matrix of slopes."""
    logger.info("Printing covariance matrix of the estimates to {ofh}.")
    np.savetxt(ofh, ldscore_reg.coef_cov)


def _print_delete_values(ldscore_reg, ofh):
    """Prints block jackknife delete-k values"""
    logger.info("Printing block jackknife delete values to {ofh}.")
    np.savetxt(ofh, ldscore_reg.tot_delete_values)


def _print_part_delete_values(ldscore_reg, ofh):
    """Prints partitioned block jackknife delete-k values"""
    logger.info("Printing partitioned block jackknife delete values to {ofh}.")
    np.savetxt(ofh, ldscore_reg.part_delete_values)


def _merge_and_log(ld, sumstats, noun):
    """Wrap smart merge with log messages about # of SNPs."""
    sumstats = smart_merge(ld, sumstats)
    msg = f"After merging with {noun}, {len(sumstats)} SNPs remain."
    if len(sumstats) == 0:
        logger.critical(msg)
        raise ValueError(msg.format(N=len(sumstats), F=noun))
    else:
        logger.info(msg)

    return sumstats


def _read_ld_sumstats(args, fh: Path, alleles=False, dropna=True):
    sumstats = _read_sumstats(args, fh, alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args)
    n_annot = len(ref_ld.columns) - 1
    M_annot = _read_M(args, n_annot)
    M_annot, ref_ld, novar_cols = _check_variance(M_annot, ref_ld)
    w_ld = _read_w_ld(args)
    sumstats = _merge_and_log(ref_ld, sumstats, "reference panel LD")
    sumstats = _merge_and_log(sumstats, w_ld, "regression SNP LD")
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.columns[1 : len(ref_ld.columns)]
    return M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols


def cell_type_specific(args):
    """Cell type specific analysis"""
    check_args(args)
    args = copy.deepcopy(args)
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1

    M_annot_all_regr, w_ld_cname, ref_ld_cnames_all_regr, sumstats, novar_cols = (
        _read_ld_sumstats(args, args.h2_cts)
    )
    M_tot = np.sum(M_annot_all_regr)
    _check_ld_condnum(args, ref_ld_cnames_all_regr)
    _warn_length(sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    if args.chisq_max is None:
        chisq_max = max(0.001 * sumstats.N.max(), 80)
    else:
        chisq_max = args.chisq_max

    ii = np.ravel(sumstats.Z**2 < chisq_max)
    sumstats = sumstats.iloc[ii, :]
    logger.info(
        f"Removed {n_snp - np.sum(ii)} SNPs with chi^2 > {chisq_max} ({np.sum(ii)} SNPs remain)"
    )
    n_snp = np.sum(ii)  # lambdas are late-binding, so this works
    ref_ld_all_regr = np.array(sumstats[ref_ld_cnames_all_regr]).reshape(
        (len(sumstats), -1)
    )
    chisq = np.array(sumstats.Z**2)
    keep_snps = sumstats[["SNP"]]

    results_columns = [
        "Name",
        "Coefficient",
        "Coefficient_std_error",
        "Coefficient_P_value",
    ]
    results_data = []
    for name, ct_ld_chr in [x.split() for x in open(args.ref_ld_chr_cts).readlines()]:
        ref_ld_cts_allsnps = _read_chr_split_files(
            ct_ld_chr, None, "cts reference panel LD Score", ps.ldscore_fromlist
        )
        logger.info("Performing regression.")
        ref_ld_cts = np.array(
            pd.merge(keep_snps, ref_ld_cts_allsnps, on="SNP", how="left").iloc[:, 1:]
        )
        if np.any(np.isnan(ref_ld_cts)):
            raise ValueError(
                "Missing some LD scores from cts files. Are you sure all SNPs in ref-ld-chr are also in ref-ld-chr-cts"
            )

        ref_ld = np.hstack([ref_ld_cts, ref_ld_all_regr])
        M_cts = ps.M_fromlist(_splitp(ct_ld_chr), _N_CHR, common=(not args.not_M_5_50))
        M_annot = np.hstack([M_cts, M_annot_all_regr])
        hsqhat = reg.Hsq(
            reshape_array(chisq, n_snp),
            ref_ld,
            reshape_array(sumstats[w_ld_cname], n_snp),
            reshape_array(sumstats.N, n_snp),
            M_annot,
            n_blocks=n_blocks,
            intercept=args.intercept_h2,
            twostep=None,
            old_weights=True,
        )
        coef, coef_se = hsqhat.coef[0], hsqhat.coef_se[0]
        results_data.append((name, coef, coef_se, stats.norm.sf(coef / coef_se)))
        if args.print_all_cts:
            for i in range(1, len(ct_ld_chr.split(","))):
                coef, coef_se = hsqhat.coef[i], hsqhat.coef_se[i]
                results_data.append(
                    (name + "_" + str(i), coef, coef_se, stats.norm.sf(coef / coef_se))
                )

    df_results = pd.DataFrame(data=results_data, columns=results_columns)
    df_results.sort_values(by="Coefficient_P_value", inplace=True)
    df_results.to_csv(str(args.out) + ".cell_type_results.txt", sep="\t", index=False)
    logger.info(f"Results printed to {args.out}.cell_type_results.tx")


def reshape_array(matrix, n_snps: int) -> np.array:
    """method to reshape the provided matrix"""
    return np.array(matrix).reshape((n_snps, 1))


def estimate_h2(args) -> reg.Hsq:
    """Estimate h2 and partitioned h2."""

    check_args(args)
    args = copy.deepcopy(args)
    # if args.samp_prev is not None and args.pop_prev is not None:
    #     args.samp_prev, args.pop_prev = list(
    #         map(float, [args.samp_prev, args.pop_prev])
    #     )
    # if args.intercept_h2 is not None:
    #     args.intercept_h2 = float(args.intercept_h2)
    # if args.no_intercept:
    #     args.intercept_h2 = 1
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, args.h2
    )
    ref_ld = np.array(sumstats[ref_ld_cnames])
    _check_ld_condnum(args, ref_ld_cnames)
    _warn_length(sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    chisq_max = args.chisq_max
    old_weights = False
    if n_annot == 1:
        if args.two_step is None and args.intercept_h2 is None:
            args.two_step = 30
    else:
        old_weights = True
        if args.chisq_max is None:
            chisq_max = max(0.001 * sumstats.N.max(), 80)

    chisq = reshape_array(sumstats.Z**2, n_snp)
    if chisq_max is not None:
        ii = np.ravel(chisq < chisq_max)
        sumstats = sumstats.iloc[ii, :]
        logger.info(
            f"Removed {n_snp - np.sum(ii)} SNPs with chi^2 > {chisq_max} ({np.sum(ii)} SNPs remain)"
        )
        n_snp = np.sum(ii)  # lambdas are late-binding, so this works
        ref_ld = np.array(sumstats[ref_ld_cnames])
        chisq = chisq[ii].reshape((n_snp, 1))

    if args.two_step is not None:
        logger.info(f"Using two-step estimator with cutoff at {args.two_step}.")

    hsqhat = reg.Hsq(
        chisq,
        ref_ld,
        reshape_array(sumstats[w_ld_cname], n_snp),
        reshape_array(sumstats.N, n_snp),
        M_annot,
        n_blocks=n_blocks,
        intercept=args.intercept_h2,
        twostep=args.two_step,
        old_weights=old_weights,
    )

    if args.print_cov:
        _print_cov(hsqhat, ps.sub_chr(args.out, ".cov"))
    if args.print_delete_vals:
        _print_delete_values(hsqhat, ps.sub_chr(args.out, ".delete"))
        _print_part_delete_values(hsqhat, ps.sub_chr(args.out, ".part_delete"))

    hsqhat.summary(
        ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, overlap=args.overlap_annot
    )

    # write output to a file
    hsqhat.write_h2_results(args.out)

    if args.overlap_annot:
        overlap_matrix, M_tot = _read_annot(args)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = hsqhat._overlap_output(
            ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients
        )
        df_results.to_csv(str(args.out) + ".results", sep="\t", index=False)
        logger.info(f"Results printed to {args.out}.results")

    return hsqhat


def estimate_rg(args):
    """Estimate rg between trait 1 and a list of other traits."""
    check_args(args)

    args = copy.deepcopy(args)

    n_pheno = len(args.rg)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev = list(
        map(
            functools.partial(_split_or_none, phenotype_count=n_pheno),
            (args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev),
        )
    )
    list(
        map(
            functools.partial(_check_arg_len, phenotype_count=n_pheno),
            (
                (args.intercept_h2, "--intercept-h2"),
                (args.intercept_gencov, "--intercept-gencov"),
                (args.samp_prev, "--samp-prev"),
                (args.pop_prev, "--pop-prev"),
            ),
        )
    )
    if args.no_intercept:
        args.intercept_h2 = [1 for _ in range(n_pheno)]
        args.intercept_gencov = [0 for _ in range(n_pheno)]
    p1 = args.rg[0]
    out_prefix = Path(args.out)
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(
        args, p1, alleles=True, dropna=True
    )
    RG = []
    n_annot = M_annot.shape[1]
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    if args.two_step is not None:
        logger.info(f"Using two-step estimator with cutoff at {args.two_step}.")

    for i, p2 in enumerate(args.rg[1:n_pheno]):
        logger.info(f"Computing rg for phenotype {i + 2}/{len(args.rg)}")
        try:
            loop = _read_other_sumstats(args, p2, sumstats, ref_ld_cnames)
            rghat = _rg(loop, args, M_annot, ref_ld_cnames, w_ld_cname, i)
            RG.append(rghat)
            _print_gencor(args, rghat, ref_ld_cnames, i, args.rg, i == 0)
            if args.print_cov:
                _print_rg_cov(rghat, f"{out_prefix}.cov")
            if args.print_delete_vals:
                _print_rg_delete_values(rghat, f"{out_prefix}.delete_vals")

        except Exception:  # keep going if phenotype 50/100 causes an error

            logger.critical(
                f"ERROR computing rg for phenotype {i + 2}/{len(args.rg)}, from file {args.rg[i + 1]}."
            )
            ex_type, ex, tb = sys.exc_info()
            logger.critical(traceback.format_exc(ex) + "\n")
            if len(RG) <= i:  # if exception raised before appending to RG
                RG.append(None)

    logger.info(
        "\nSummary of Genetic Correlation Results\n"
        + _get_rg_table(args.rg, RG, out_prefix, args)
    )
    return RG


def _read_other_sumstats(args, p2, sumstats, ref_ld_cnames):
    loop = _read_sumstats(args, p2, alleles=True, dropna=False)
    loop = _merge_sumstats_sumstats(args, sumstats, loop)
    loop = loop.dropna(how="any")
    alleles = loop.A1 + loop.A2 + loop.A1x + loop.A2x
    if not args.no_check_alleles:
        loop = _select_and_log(loop, _filter_alleles(alleles))
        loop["Z2"] = _align_alleles(loop.Z2, alleles)

    loop = loop.drop(["A1", "A1x", "A2", "A2x"], axis=1)
    _check_ld_condnum(args, loop[ref_ld_cnames])
    _warn_length(loop)
    return loop


def _get_rg_table(rg_paths: list[Path], RG, output_prefix: Path, args):
    """Print a table of genetic correlations."""
    output_path = output_prefix.parent / f"{output_prefix.name}.rg_results"

    t = lambda attr: lambda obj: getattr(obj, attr, "NA")
    x = pd.DataFrame()
    x["p1"] = [rg_paths[0] for i in range(1, len(rg_paths))]
    x["p2"] = rg_paths[1 : len(rg_paths)]
    x["rg"] = list(map(t("rg_ratio"), RG))
    x["se"] = list(map(t("rg_se"), RG))
    x["z"] = list(map(t("z"), RG))
    x["p"] = list(map(t("p"), RG))
    if (
        args.samp_prev is not None
        and args.pop_prev is not None
        and all((i is not None for i in args.samp_prev))
        and all((it is not None for it in args.pop_prev))
    ):

        c = list(
            map(
                lambda x, y: reg.h2_obs_to_liab(1, x, y),
                args.samp_prev[1:],
                args.pop_prev[1:],
            )
        )
        x["h2_liab"] = list(
            map(lambda x, y: x * y, c, list(map(t("tot"), list(map(t("hsq2"), RG)))))
        )
        x["h2_liab_se"] = list(
            map(lambda x, y: x * y, c, list(map(t("tot_se"), list(map(t("hsq2"), RG)))))
        )
    else:
        x["h2_obs"] = list(map(t("tot"), list(map(t("hsq2"), RG))))
        x["h2_obs_se"] = list(map(t("tot_se"), list(map(t("hsq2"), RG))))

    x["h2_int"] = list(map(t("intercept"), list(map(t("hsq2"), RG))))
    x["h2_int_se"] = list(map(t("intercept_se"), list(map(t("hsq2"), RG))))
    x["gcov_int"] = list(map(t("intercept"), list(map(t("gencov"), RG))))
    x["gcov_int_se"] = list(map(t("intercept_se"), list(map(t("gencov"), RG))))
    x.to_csv(output_path, sep="\t", index=None)
    return x.to_string(header=True, index=False) + "\n"


def _print_gencor(args, rghat, ref_ld_cnames, i, rg_paths, print_hsq1):
    l = lambda x: x + "".join(["-" for i in range(len(x.replace("\n", "")))])
    P = [args.samp_prev[0], args.samp_prev[i + 1]]
    K = [args.pop_prev[0], args.pop_prev[i + 1]]
    if args.samp_prev is None and args.pop_prev is None:
        args.samp_prev = [None, None]
        args.pop_prev = [None, None]
    if print_hsq1:
        logger.info(l("\nHeritability of phenotype 1\n"))
        logger.info(rghat.hsq1.summary(ref_ld_colnames=ref_ld_cnames, P=P[0], K=K[0]))

    # write the heritability of phenotype to a file
    logger.info(
        l("\nHeritability of phenotype {I}/{N}\n".format(I=i + 2, N=len(rg_paths)))
    )
    rghat.hsq2.summary(ref_ld_colnames=ref_ld_cnames, P=P[1], K=K[1])
    # write the genetic correlation to a file
    logger.info(l("\nGenetic Covariance\n"))
    rghat.gencov.summary(ref_ld_colnames=ref_ld_cnames, P=P, K=K)
    logger.info(l("\nGenetic Correlation\n"))
    rghat.summary()


def _merge_sumstats_sumstats(args, sumstats1, sumstats2):
    """Merge two sets of summary statistics."""
    sumstats1.rename(columns={"N": "N1", "Z": "Z1"}, inplace=True)
    sumstats2.rename(
        columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2"}, inplace=True
    )
    x = _merge_and_log(sumstats1, sumstats2, "summary statistics")
    return x


def _filter_alleles(alleles):
    """Remove bad variants (mismatched alleles, non-SNPs, strand ambiguous)."""
    ii = alleles.apply(lambda y: y in MATCH_ALLELES)
    return ii


def _align_alleles(z, alleles):
    """Align Z1 and Z2 to same choice of ref allele (allowing for strand flip)."""
    try:
        z *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
    except KeyError as e:
        msg = "Incompatible alleles in .sumstats files: %s. " % e.args
        msg += "Did you forget to use --merge-alleles with munge_sumstats.py?"
        raise KeyError(msg)
    return z


def _rg(sumstats, args, M_annot, ref_ld_cnames, w_ld_cname, i):
    """Run the regressions."""
    n_snp = len(sumstats)

    if args.chisq_max is not None:
        ii = sumstats.Z1**2 * sumstats.Z2**2 < args.chisq_max**2
        n_snp = np.sum(ii)  # lambdas are late binding, so this works
        sumstats = sumstats[ii]
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats[ref_ld_cnames].to_numpy()
    intercepts = [
        args.intercept_h2[0],
        args.intercept_h2[i + 1],
        args.intercept_gencov[i + 1],
    ]
    rghat = reg.RG(
        reshape_array(sumstats.Z1, n_snp),
        reshape_array(sumstats.Z2, n_snp),
        ref_ld,
        reshape_array(sumstats[w_ld_cname], n_snp),
        reshape_array(sumstats.N1, n_snp),
        reshape_array(sumstats.N2, n_snp),
        M_annot,
        intercept_hsq1=intercepts[0],
        intercept_hsq2=intercepts[1],
        intercept_gencov=intercepts[2],
        n_blocks=n_blocks,
        twostep=args.two_step,
    )

    return rghat


def _print_rg_delete_values(rg, fh):
    """Print block jackknife delete values."""
    _print_delete_values(rg.hsq1, fh + ".hsq1.delete")
    _print_delete_values(rg.hsq2, fh + ".hsq2.delete")
    _print_delete_values(rg.gencov, fh + ".gencov.delete")


def _print_rg_cov(rghat, fh):
    """Print covariance matrix of estimates."""
    _print_cov(rghat.hsq1, fh + ".hsq1.cov")
    _print_cov(rghat.hsq2, fh + ".hsq2.cov")
    _print_cov(rghat.gencov, fh + ".gencov.cov")


def _split_or_none(x, phenotype_count: int):
    if x is not None:
        y = list(map(float, x.replace("N", "-").split(",")))
    else:
        y = [None for _ in range(phenotype_count)]
    return y


def _check_arg_len(x, phenotype_count: int):
    x, m = x
    if len(x) != phenotype_count:
        raise ValueError(
            "{M} must have the same number of arguments as --rg/--h2.".format(M=m)
        )
