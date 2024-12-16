from pathlib import Path
from itertools import product
import logging
import numpy as np
import bitarray as ba
import pandas as pd
import ldsc.ldscore.parse as ps
import ldsc.ldscore.sumstats as sumstats
import ldsc.ldscore.regressions as reg
from ldsc.logger import LDSCLogger
from .ldsc_check_args import check_args


logger: logging.Logger = LDSCLogger.get_logger(__name__)


def getBlockLefts(coords, max_dist):
    """
    Converts coordinates + max block length to the a list of coordinates of the leftmost
    SNPs to be included in blocks.

    Parameters
    ----------
    coords : array
        Array of coordinates. Must be sorted.
    max_dist : float
        Maximum distance between SNPs included in the same window.

    Returns
    -------
    block_left : 1D np.ndarray with same length as block_left
        block_left[j] :=  min{k | dist(j, k) < max_dist}.

    """
    M = len(coords)
    j = 0
    block_left = np.zeros(M)
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1

        block_left[i] = j

    return block_left


def block_left_to_right(block_left):
    """
    Converts block lefts to block rights.

    Parameters
    ----------
    block_left : array
        Array of block lefts.

    Returns
    -------
    block_right : 1D np.ndarray with same length as block_left
        block_right[j] := max {k | block_left[k] <= j}

    """
    M = len(block_left)
    j = 0
    block_right = np.zeros(M)
    for i in range(M):
        while j < M and block_left[j] <= i:
            j += 1

        block_right[i] = j

    return block_right


class __GenotypeArrayInMemory__(object):
    """
    Parent class for various classes containing interfaces for files with genotype
    matrices, e.g., plink .bed files, etc
    """

    def __init__(
        self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None
    ):
        self.m = len(snp_list.IDList)
        self.n = n
        self.keep_snps = keep_snps
        self.keep_indivs = keep_indivs
        self.df = np.array(snp_list.df[["CHR", "SNP", "BP", "CM"]])
        self.colnames = ["CHR", "SNP", "BP", "CM"]
        self.mafMin = mafMin if mafMin is not None else 0
        self._currentSNP = 0
        (self.nru, self.geno) = self.__read__(fname, self.m, n)
        # filter individuals
        if keep_indivs is not None:
            keep_indivs = np.array(keep_indivs, dtype="int")
            if np.any(keep_indivs > self.n):
                raise ValueError("keep_indivs indices out of bounds")

            (self.geno, self.m, self.n) = self.__filter_indivs__(
                self.geno, keep_indivs, self.m, self.n
            )

            if self.n > 0:
                print("After filtering, {n} individuals remain".format(n=self.n))
            else:
                raise ValueError("After filtering, no individuals remain")

        # filter SNPs
        if keep_snps is not None:
            keep_snps = np.array(keep_snps, dtype="int")
            if np.any(keep_snps > self.m):  # if keep_snps is None, this returns False
                raise ValueError("keep_snps indices out of bounds")

        (self.geno, self.m, self.n, self.kept_snps, self.freq) = (
            self.__filter_snps_maf__(self.geno, self.m, self.n, self.mafMin, keep_snps)
        )

        if self.m > 0:
            print("After filtering, {m} SNPs remain".format(m=self.m))
        else:
            raise ValueError("After filtering, no SNPs remain")

        self.df = self.df[self.kept_snps, :]
        self.maf = np.minimum(self.freq, np.ones(self.m) - self.freq)
        self.sqrtpq = np.sqrt(self.freq * (np.ones(self.m) - self.freq))
        self.df = np.c_[self.df, self.maf]
        self.colnames.append("MAF")

    def __read__(self, fname, m, n):
        raise NotImplementedError

    def __filter_indivs__(geno, keep_indivs, m, n):
        raise NotImplementedError

    def __filter_maf_(geno, m, n, maf):
        raise NotImplementedError

    def ldScoreVarBlocks(self, block_left, c, annot=None):
        """Computes an unbiased estimate of L2(j) for j=1,..,M."""
        # func = lambda x: self.__l2_unbiased__(x)
        snp_getter = self.nextSNPs
        return self.__corSumVarBlocks__(
            block_left, c, self.__l2_unbiased__, snp_getter, annot
        )

    def ldScoreBlockJackknife(self, block_left, c, annot=None, jN=10):
        # func = lambda x: np.square(x)
        snp_getter = self.nextSNPs
        return self.__corSumBlockJackknife__(
            block_left, c, np.square, snp_getter, annot, jN
        )

    def __l2_unbiased__(self, x):
        denom = self.n - 2 if self.n > 2 else self.n  # allow n<2 for testing purposes
        sq = np.square(x)
        return sq - (1 - sq) / denom

    # general methods for calculating sums of Pearson correlation coefficients
    def __corSumVarBlocks__(self, block_left, c, func, snp_getter, annot=None):
        """
        Parameters
        ----------
        block_left : np.ndarray with shape (M, )
            block_left[i] = index of leftmost SNP included in LD Score of SNP i.
            if c > 1, then only entries that are multiples of c are examined, and it is
            assumed that block_left[a*c+i] = block_left[a*c], except at
            the beginning of the chromosome where the 0th SNP is included in the window.

        c : int
            Chunk size.
        func : function
            Function to be applied to the genotype correlation matrix. Before dotting with
            annot. Examples: for biased L2, np.square. For biased L4,
            lambda x: np.square(np.square(x)). For L1, lambda x: x.
        snp_getter : function(int)
            The method to be used to get the next SNPs (normalized genotypes? Normalized
            genotypes with the minor allele as reference allele? etc)
        annot: numpy array with shape (m,n_a)
            SNP annotations.

        Returns
        -------
        cor_sum : np.ndarray with shape (M, num_annots)
            Estimates.

        """
        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c) * c
        if annot is None:
            annot = np.ones((m, 1))
        else:
            annot_m = annot.shape[0]
            if annot_m != self.m:
                raise ValueError("Incorrect number of SNPs in annot")

        n_a = annot.shape[1]  # number of annotations
        cor_sum = np.zeros((m, n_a))
        # b = index of first SNP for which SNP 0 is not included in LD Score
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b / c) * c)  # round up to a multiple of c
        if b > m:
            c = 1
            b = m
        l_A = 0  # l_A := index of leftmost SNP in matrix A
        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        # chunk inside of block
        for l_B in range(0, b, c):  # l_B := index of leftmost SNP in matrix B
            B = A[:, l_B : l_B + c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A : l_A + b, :] += np.dot(rfuncAB, annot[l_B : l_B + c, :])
        # chunk to right of block
        b0 = b
        md = int(c * np.floor(m / c))
        end = md + 1 if md != m else md
        for l_B in range(b0, end, c):
            # check if the annot matrix is all zeros for this block + chunk
            # this happens w/ sparse categories (i.e., pathways)
            # update the block
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                # block_size can't increase more than c
                # block_size can't be less than c unless it is zero
                # both of these things make sense
                A = np.hstack((A[:, old_b - b + c : old_b], B))
                l_A += old_b - b + c
            elif l_B == b0 and b > 0:
                A = A[:, b0 - b : b0]
                l_A = b0 - b
            elif b == 0:  # no SNPs to left in window, e.g., after a sequence gap
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))

            B = snp_getter(c)
            p1 = np.all(annot[l_A : l_A + b, :] == 0)
            p2 = np.all(annot[l_B : l_B + c, :] == 0)
            if p1 and p2:
                continue

            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A : l_A + b, :] += np.dot(rfuncAB, annot[l_B : l_B + c, :])
            cor_sum[l_B : l_B + c, :] += np.dot(annot[l_A : l_A + b, :].T, rfuncAB).T
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = func(rfuncBB)
            cor_sum[l_B : l_B + c, :] += np.dot(rfuncBB, annot[l_B : l_B + c, :])

        return cor_sum


class PlinkBEDFile(__GenotypeArrayInMemory__):
    """
    Interface for Plink .bed format
    """

    def __init__(
        self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None
    ):
        self._bedcode = {
            2: ba.bitarray("11"),
            9: ba.bitarray("10"),
            1: ba.bitarray("01"),
            0: ba.bitarray("00"),
        }

        __GenotypeArrayInMemory__.__init__(
            self,
            fname,
            n,
            snp_list,
            keep_snps=keep_snps,
            keep_indivs=keep_indivs,
            mafMin=mafMin,
        )

    def __read__(self, fname, m, n):
        fname = Path(fname)
        if fname.suffix != ".bed":
            raise ValueError(".bed filename must end in .bed")

        with open(fname, "rb") as fh:
            magicNumber = ba.bitarray(endian="little")
            magicNumber.fromfile(fh, 2)
            bedMode = ba.bitarray(endian="little")
            bedMode.fromfile(fh, 1)
            e = (4 - n % 4) if n % 4 != 0 else 0
            nru = n + e
            self.nru = nru
            # check magic number
            if magicNumber != ba.bitarray("0011011011011000"):
                raise IOError("Magic number from Plink .bed file not recognized")

            if bedMode != ba.bitarray("10000000"):
                raise IOError("Plink .bed file must be in default SNP-major mode")

            # check file length
            self.geno = ba.bitarray(endian="little")
            self.geno.fromfile(fh)
            self.__test_length__(self.geno, self.m, self.nru)
            return (self.nru, self.geno)

    def __test_length__(self, geno, m, nru):
        exp_len = 2 * m * nru
        real_len = len(geno)
        if real_len != exp_len:
            s = "Plink .bed file has {n1} bits, expected {n2}"
            raise IOError(s.format(n1=real_len, n2=exp_len))

    def __filter_indivs__(self, geno, keep_indivs, m, n):
        n_new = len(keep_indivs)
        e = (4 - n_new % 4) if n_new % 4 != 0 else 0
        nru_new = n_new + e
        nru = self.nru
        z = ba.bitarray(m * 2 * nru_new, endian="little")
        z.setall(0)
        for e, i in enumerate(keep_indivs):
            z[2 * e :: 2 * nru_new] = geno[2 * i :: 2 * nru]
            z[2 * e + 1 :: 2 * nru_new] = geno[2 * i + 1 :: 2 * nru]

        self.nru = nru_new
        return (z, m, n_new)

    def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
        """
        Credit to Chris Chang and the Plink2 developers for this algorithm
        Modified from plink_filter.c
        https://github.com/chrchang/plink-ng/blob/master/plink_filter.c

        Genotypes are read forwards (since we are cheating and using endian="little")

        A := (genotype) & 1010...
        B := (genotype) & 0101...
        C := (A >> 1) & B

        Then

        a := A.count() = missing ct + hom major ct
        b := B.count() = het ct + hom major ct
        c := C.count() = hom major ct

        Which implies that

        missing ct = a - c
        # of indivs with nonmissing genotype = n - a + c
        major allele ct = b + c
        major allele frequency = (b+c)/(2*(n-a+c))
        het ct + missing ct = a + b - 2*c

        Why does bitarray not have >> ????

        """
        nru = self.nru
        m_poly = 0
        y = ba.bitarray()
        if keep_snps is None:
            keep_snps = range(m)
        kept_snps = []
        freq = []
        for e, j in enumerate(keep_snps):
            z = geno[2 * nru * j : 2 * nru * (j + 1)]
            A = z[0::2]
            a = A.count()
            B = z[1::2]
            b = B.count()
            c = (A & B).count()
            major_ct = b + c  # number of copies of the major allele
            n_nomiss = n - a + c  # number of individuals with nonmissing genotypes
            f = major_ct / (2 * n_nomiss) if n_nomiss > 0 else 0
            het_miss_ct = (
                a + b - 2 * c
            )  # remove SNPs that are only either het or missing
            if np.minimum(f, 1 - f) > mafMin and het_miss_ct < n:
                freq.append(f)
                y += z
                m_poly += 1
                kept_snps.append(j)

        return (y, m_poly, n, kept_snps, freq)

    def nextSNPs(self, b, minorRef=None):
        """
        Unpacks the binary array of genotypes and returns an n x b matrix of floats of
        normalized genotypes for the next b SNPs, where n := number of samples.

        Parameters
        ----------
        b : int
            Number of SNPs to return.
        minorRef: bool, default None
            Should we flip reference alleles so that the minor allele is the reference?
            (This is useful for computing l1 w.r.t. minor allele).

        Returns
        -------
        X : np.array with dtype float64 with shape (n, b), where n := number of samples
            Matrix of genotypes normalized to mean zero and variance one. If minorRef is
            not None, then the minor allele will be the positive allele (i.e., two copies
            of the minor allele --> a positive number).

        """

        try:
            b = int(b)
            if b <= 0:
                raise ValueError("b must be > 0")
        except TypeError:
            raise TypeError("b must be an integer")

        if self._currentSNP + b > self.m:
            s = "{b} SNPs requested, {k} SNPs remain"
            raise ValueError(s.format(b=b, k=(self.m - self._currentSNP)))

        c = self._currentSNP
        n = self.n
        nru = self.nru
        slice = self.geno[2 * c * nru : 2 * (c + b) * nru]
        X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T
        X = X[0:n, :]
        Y = np.zeros(X.shape)
        for j in range(0, b):
            newsnp = X[:, j]
            ii = newsnp != 9
            avg = np.mean(newsnp[ii])
            newsnp[np.logical_not(ii)] = avg
            denom = np.std(newsnp)
            if denom == 0:
                denom = 1

            if minorRef is not None and self.freq[self._currentSNP + j] > 0.5:
                denom = denom * -1

            Y[:, j] = (newsnp - avg) / denom

        self._currentSNP += b
        return Y


def _remove_dtype(x):
    """Removes dtype: float64 and dtype: int64 from pandas printouts"""
    x = str(x)
    x = x.replace("\ndtype: int64", "")
    x = x.replace("\ndtype: float64", "")
    return x


def _filter(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        x = ps.FilterFile(fname)
        logger.info(f"Read list of {len(x.IDList)} {noun} to {verb} {fname}")
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            logger.info(f"After merging, {len_merged_list} {noun}")
        else:
            logger.critical(f"No {noun} retained for analysis")
            raise ValueError(f"No {noun} retained for analysis")

        return merged_list


def annot_sort_key(s):
    """For use with --cts-bin. Fixes weird pandas crosstab column order."""
    if isinstance(s, tuple):
        s = [x.split("_")[0] for x in s]
        s = [float(x) if x != "min" else -float("inf") for x in s]
    else:  # type(s) = str:
        s = s.split("_")[0]
        if s == "min":
            s = float("-inf")
        else:
            s = float(s)

    return s


def ldscore(args):
    """
    Wrapper function for estimating l1, l1^2, l2 and l4 (+ optionally standard errors) from
    reference panel genotypes.

    Annot format is
    chr snp bp cm <annotations>

    """
    check_args(args)

    if args.bfile:
        snp_file, snp_obj = args.bfile + ".bim", ps.PlinkBIMFile
        ind_file, ind_obj = args.bfile + ".fam", ps.PlinkFAMFile
        array_file, array_obj = args.bfile + ".bed", PlinkBEDFile

    # read bim/snp
    array_snps = snp_obj(snp_file)
    m = len(array_snps.IDList)
    logger.info(f"Read list of {m} SNPs from {snp_file}")
    if args.annot is not None:  # read --annot
        try:
            if args.thin_annot:  # annot file has only annotations
                annot = ps.ThinAnnotFile(args.annot)
                n_annot, ma = len(annot.df.columns), len(annot.df)
                logger.info(
                    f"Read {n_annot} annotations for {ma} SNPs from {args.annot}"
                )
                annot_matrix = annot.df.values
                annot_colnames = annot.df.columns
                keep_snps = None
            else:
                annot = ps.AnnotFile(args.annot)
                n_annot, ma = len(annot.df.columns) - 4, len(annot.df)
                logger.info(
                    f"Read {n_annot} annotations for {ma} SNPs from {args.annot}"
                )
                annot_matrix = np.array(annot.df.iloc[:, 4:])
                annot_colnames = annot.df.columns[4:]
                keep_snps = None
                if np.any(annot.df.SNP.values != array_snps.df.SNP.values):
                    logger.critical(
                        "The .annot file must contain the same SNPs in the same"
                        + " order as the .bim file."
                    )
                    raise ValueError(
                        "The .annot file must contain the same SNPs in the same"
                        + " order as the .bim file."
                    )
        except Exception as e:
            logger.critical("Error parsing .annot file")
            raise e

    elif args.extract is not None:  # --extract
        keep_snps = _filter(args.extract, "SNPs", "include", array_snps)
        annot_matrix, annot_colnames, n_annot = None, None, 1

    elif args.cts_bin is not None and args.cts_breaks is not None:  # --cts-bin
        cts_fnames = sumstats._splitp(args.cts_bin)  # read filenames
        args.cts_breaks = args.cts_breaks.replace(
            "N", "-"
        )  # replace N with negative sign
        try:  # split on x
            breaks = [
                [float(x) for x in y.split(",")] for y in args.cts_breaks.split("x")
            ]
        except ValueError as e:
            raise ValueError(
                "--cts-breaks must be a comma-separated list of numbers: " + str(e.args)
            )

        if len(breaks) != len(cts_fnames):
            raise ValueError(
                "Need to specify one set of breaks for each file in --cts-bin."
            )

        if args.cts_names:
            cts_colnames = [str(x) for x in args.cts_names.split(",")]
            if len(cts_colnames) != len(cts_fnames):
                msg = "Must specify either no --cts-names or one value for each file in --cts-bin."
                raise ValueError(msg)

        else:
            cts_colnames = ["ANNOT" + str(i) for i in range(len(cts_fnames))]

        logger.info("Reading numbers with which to bin SNPs from {args.cts_bin}")

        cts_levs = []
        full_labs = []
        for i, fh in enumerate(cts_fnames):
            vec = ps.read_cts(cts_fnames[i], array_snps.df.SNP.values)

            max_cts = np.max(vec)
            min_cts = np.min(vec)
            cut_breaks = list(breaks[i])
            name_breaks = list(cut_breaks)
            if np.all(cut_breaks >= max_cts) or np.all(cut_breaks <= min_cts):
                raise ValueError(
                    "All breaks lie outside the range of the cts variable."
                )

            if np.all(cut_breaks <= max_cts):
                name_breaks.append(max_cts)
                cut_breaks.append(max_cts + 1)

            if np.all(cut_breaks >= min_cts):
                name_breaks.append(min_cts)
                cut_breaks.append(min_cts - 1)

            name_breaks.sort()
            cut_breaks.sort()
            n_breaks = len(cut_breaks)
            # so that col names are consistent across chromosomes with different max vals
            name_breaks[0] = "min"
            name_breaks[-1] = "max"
            name_breaks = [str(x) for x in name_breaks]
            labs = [
                name_breaks[i] + "_" + name_breaks[i + 1] for i in range(n_breaks - 1)
            ]
            cut_vec = pd.Series(pd.cut(vec, bins=cut_breaks, labels=labs))
            cts_levs.append(cut_vec)
            full_labs.append(labs)

        annot_matrix = pd.concat(cts_levs, axis=1)
        annot_matrix.columns = cts_colnames
        # crosstab -- for now we keep empty columns
        annot_matrix = pd.crosstab(
            annot_matrix.index,
            [annot_matrix[i] for i in annot_matrix.columns],
            dropna=False,
            colnames=annot_matrix.columns,
        )

        # add missing columns
        if len(cts_colnames) > 1:
            for x in product(*full_labs):
                if x not in annot_matrix.columns:
                    annot_matrix[x] = 0
        else:
            for x in full_labs[0]:
                if x not in annot_matrix.columns:
                    annot_matrix[x] = 0

        annot_matrix = annot_matrix[sorted(annot_matrix.columns, key=annot_sort_key)]
        if len(cts_colnames) > 1:
            # flatten multi-index
            annot_colnames = [
                "_".join([cts_colnames[i] + "_" + b for i, b in enumerate(c)])
                for c in annot_matrix.columns
            ]
        else:
            annot_colnames = [cts_colnames[0] + "_" + b for b in annot_matrix.columns]

        annot_matrix = np.matrix(annot_matrix)
        keep_snps = None
        n_annot = len(annot_colnames)
        if np.any(np.sum(annot_matrix, axis=1) == 0):
            # This exception should never be raised. For debugging only.
            raise ValueError(
                "Some SNPs have no annotation in --cts-bin. This is a bug!"
            )

    else:
        annot_matrix, annot_colnames, keep_snps = (
            None,
            None,
            None,
        )
        n_annot = 1

    # read fam
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    logger.info(f"Read list of {n} individuals from {ind_file}".format(n=n, f=ind_file))
    # read keep_indivs
    if args.keep:
        keep_indivs = _filter(args.keep, "individuals", "include", array_indivs)
    else:
        keep_indivs = None

    # read genotype array
    logger.info(f"Reading genotypes from {array_file}")
    geno_array = array_obj(
        array_file,
        n,
        array_snps,
        keep_snps=keep_snps,
        keep_indivs=keep_indivs,
        mafMin=args.maf,
    )

    # filter annot_matrix down to only SNPs passing MAF cutoffs
    if annot_matrix is not None:
        annot_keep = geno_array.kept_snps
        annot_matrix = annot_matrix[annot_keep, :]

    # determine block widths
    x = np.array((args.ld_wind_snps, args.ld_wind_kb, args.ld_wind_cm), dtype=bool)
    if np.sum(x) != 1:
        raise ValueError("Must specify exactly one --ld-wind option")

    if args.ld_wind_snps:
        max_dist = args.ld_wind_snps
        coords = np.array(range(geno_array.m))
    elif args.ld_wind_kb:
        max_dist = args.ld_wind_kb * 1000
        coords = np.array(array_snps.df["BP"])[geno_array.kept_snps]
    elif args.ld_wind_cm:
        max_dist = args.ld_wind_cm
        coords = np.array(array_snps.df["CM"])[geno_array.kept_snps]

    block_left = getBlockLefts(coords, max_dist)
    if block_left[len(block_left) - 1] == 0 and not args.yes_really:
        error_msg = (
            "Do you really want to compute whole-chomosome LD Score? If so, set the "
        )
        error_msg += "--yes-really flag (warning: it will use a lot of time / memory)"
        raise ValueError(error_msg)

    scale_suffix = ""
    if args.pq_exp is not None:
        logger.log("Computing LD with pq ^ {args.pq_exp}.")
        logger.info(
            "Note that LD Scores with pq raised to a nonzero power are not directly comparable to normal LD Scores."
        )
        scale_suffix = f"_S{args.pq_exp}"
        pq = np.matrix(geno_array.maf * (1 - geno_array.maf)).reshape((geno_array.m, 1))
        pq = np.power(pq, args.pq_exp)

        if annot_matrix is not None:
            annot_matrix = np.multiply(annot_matrix, pq)
        else:
            annot_matrix = pq

    logger.info("Estimating LD Score.")
    lN = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_matrix)
    col_prefix = "L2"
    file_suffix = "l2"

    if n_annot == 1:
        ldscore_colnames = [col_prefix + scale_suffix]
    else:
        ldscore_colnames = [y + col_prefix + scale_suffix for y in annot_colnames]

    # print .ldscore. Output columns: CHR, BP, RS, [LD Scores]
    out_fname = str(args.out) + "." + file_suffix + ".ldscore"
    new_colnames = geno_array.colnames + ldscore_colnames
    df = pd.DataFrame.from_records(np.c_[geno_array.df, lN])
    df.columns = new_colnames
    if args.print_snps:
        if args.print_snps.endswith("gz"):
            print_snps = pd.read_csv(args.print_snps, header=None, compression="gzip")
        elif args.print_snps.endswith("bz2"):
            print_snps = pd.read_csv(args.print_snps, header=None, compression="bz2")
        else:
            print_snps = pd.read_csv(args.print_snps, header=None)
        if len(print_snps.columns) > 1:
            raise ValueError(
                "--print-snps must refer to a file with a one column of SNP IDs."
            )
        logger.info(
            f"Reading list of {len(print_snps)} SNPs for which to print LD Scores from {args.print_snps}"
        )

        print_snps.columns = ["SNP"]
        df = df.loc[df.SNP.isin(print_snps.SNP), :]
        if len(df) == 0:
            raise ValueError("After merging with --print-snps, no SNPs remain.")
        else:
            logger.info(
                f"After merging with --print-snps, LD Scores for {len(df)} SNPs will be printed."
            )

    l2_suffix = ".gz"
    logger.info(f"Writing LD Scores for {len(df)} SNPs to {out_fname}.gz")
    df.drop(["CM", "MAF"], axis=1).to_csv(
        f"{out_fname}.gz",
        sep="\t",
        header=True,
        index=False,
        float_format="%.3f",
        compression="gzip"
    )
    if annot_matrix is not None:
        M = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix, axis=0))))
        ii = geno_array.maf > 0.05
        M_5_50 = np.atleast_1d(
            np.squeeze(np.asarray(np.sum(annot_matrix[ii, :], axis=0)))
        )
    else:
        M = [geno_array.m]
        M_5_50 = [np.sum(geno_array.maf > 0.05)]

    # print .M
    with open(f"{args.out}.{file_suffix}.M", "wb") as fout_M:
        fout_M.write("\t".join(map(str, M)).encode() + b'\n')

    # print .M_5_50
    with open(f"{args.out}.{file_suffix}.M_5_50", "wb") as fout_M_5_50:
        fout_M_5_50.write("\t".join(map(str, M_5_50)).encode() + b'\n')

    # print annot matrix
    if (args.cts_bin is not None) and not args.no_print_annot:
        out_fname_annot = args.out + ".annot"
        new_colnames = geno_array.colnames + ldscore_colnames
        annot_df = pd.DataFrame(np.c_[geno_array.df, annot_matrix])
        annot_df.columns = new_colnames
        del annot_df["MAF"]
        logger.info(
            f"Writing annot matrix produced by --cts-bin to {out_fname + '.gz'}"
        )
        annot_df.to_csv(
            out_fname_annot, sep="\t", header=True, index=False, compression="gzip"
        )

    # print LD Score summary
    pd.set_option("display.max_rows", 200)
    logger.info(f"\nSummary of LD Scores in {out_fname + l2_suffix}")
    t = df.iloc[:, 4:].describe()
    logger.info(t.iloc[1:, :])

    np.seterr(divide="ignore", invalid="ignore")  # print NaN instead of weird errors
    # print correlation matrix including all LD Scores and sample MAF
    logger.info("")
    logger.info("MAF/LD Score Correlation Matrix")
    logger.info(df.iloc[:, 4:].corr())

    # print condition number
    if (
        n_annot > 1
    ):  # condition number of a column vector w/ nonzero var is trivially one
        logger.info("\nLD Score Matrix Condition Number")
        cond_num = np.linalg.cond(df.iloc[:, 5:])
        logger.info(reg.remove_brackets(str(np.matrix(cond_num))))
        if cond_num > 10000:
            logger.warning("WARNING: ill-conditioned LD Score Matrix!")

    # summarize annot matrix if there is one
    if annot_matrix is not None:
        # covariance matrix
        x = pd.DataFrame(annot_matrix, columns=annot_colnames)
        logger.info("\nAnnotation Correlation Matrix")
        logger.info(x.corr())

        # column sums
        logger.info("\nAnnotation Matrix Column Sums")
        logger.info(_remove_dtype(x.sum(axis=0)))

        # row sums
        logger.info("\nSummary of Annotation Matrix Row Sums")
        row_sums = x.sum(axis=1).describe()
        logger.info(_remove_dtype(row_sums))

    np.seterr(divide="raise", invalid="raise")
