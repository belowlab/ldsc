from pathlib import Path
from ldscore import parse as ps
import unittest
import numpy as np
import pandas as pd
import os
from nose.tools import *
from numpy.testing import assert_array_equal

DIR = Path(os.path.dirname(__file__))


class TestModules(unittest.TestCase):
    def test_series_eq(self):
        x = pd.Series([1, 2, 3])
        y = pd.Series([1, 2])
        z = pd.Series([1, 2, 4])
        assert ps.series_eq(x, x)
        assert not ps.series_eq(x, y)
        assert not ps.series_eq(x, z)

    def test_get_compression(self):
        self.assertEqual(ps.get_compression(Path("test.gz")), "gzip")
        self.assertEqual(ps.get_compression(Path("test.bz2")), "bz2")
        self.assertEqual(ps.get_compression(Path("test.asdf")), None)

    def test_read_cts(self):
        match_snps = pd.Series(["rs1", "rs2", "rs3"])
        assert_array_equal(
            ps.read_cts(DIR / "parse_test/test.cts", match_snps), [1, 2, 3]
        )
        self.assertRaises(
            ValueError,
            ps.read_cts,
            DIR / "parse_test/test.cts",
            match_snps[0:2],
        )

    def test_read_sumstats(self):
        x = ps.sumstats(DIR / "parse_test/test.sumstats", dropna=True, alleles=True)
        self.assertEqual(len(x), 1)
        assert_array_equal(x.SNP, "rs1")
        self.assertRaises(
            ValueError, ps.sumstats, DIR / "parse_test/test.l2.ldscore.gz"
        )

    def test_frq_parser(self):
        x = ps.frq_parser(DIR / "parse_test/test1.frq", compression=None)
        assert_array_equal(x.columns, ["SNP", "FRQ"])
        assert_array_equal(x.SNP, ["rs_" + str(i) for i in range(8)])
        assert_array_equal(x.FRQ, [0.01, 0.1, 0.7, 0.2, 0.2, 0.2, 0.99, 0.03])
        x = ps.frq_parser(DIR / "parse_test/test2.frq.gz", compression="gzip")
        assert_array_equal(x.columns, ["SNP", "FRQ"])
        assert_array_equal(x.SNP, ["rs_" + str(i) for i in range(8)])
        assert_array_equal(x.FRQ, [0.01, 0.1, 0.3, 0.2, 0.2, 0.2, 0.01, 0.03])


class Test_ldscore(unittest.TestCase):

    def test_ldscore(self):
        x = ps.ldscore(DIR / "parse_test/test")
        self.assertEqual(list(x["SNP"]), ["rs" + str(i) for i in range(1, 23)])
        self.assertEqual(list(x["AL2"]), list(range(1, 23)))
        self.assertEqual(list(x["BL2"]), list(range(2, 46, 2)))

    def test_ldscore_loop(self):
        x = ps.ldscore(DIR / "parse_test/test", 2)
        self.assertEqual(list(x["SNP"]), ["rs" + str(i) for i in range(1, 3)])
        self.assertEqual(list(x["AL2"]), list(range(1, 3)))
        self.assertEqual(list(x["BL2"]), list(range(2, 6, 2)))

    def test_ldscore_fromlist(self):
        fh = DIR / "parse_test/test"
        x = ps.ldscore_fromlist([fh, fh])
        assert_array_equal(x.shape, (22, 5))
        y = ps.ldscore(DIR / "parse_test/test")
        assert_array_equal(x.iloc[:, 0:3], y)
        assert_array_equal(x.iloc[:, [0, 3, 4]], y)
        self.assertRaises(
            ValueError, ps.ldscore_fromlist, [fh, DIR / "parse_test/test2"]
        )


class Test_M(unittest.TestCase):

    def test_bad_M(self):
        self.assertRaises(ValueError, ps.M, DIR / "parse_test/test_bad")

    def test_M(self):
        x = ps.M(DIR / "parse_test/test")
        assert_array_equal(x.shape, (1, 3))
        assert_array_equal(x, [[1000, 2000, 3000]])

    def test_M_loop(self):
        x = ps.M(DIR / "parse_test/test", 2)
        assert_array_equal(x.shape, (1, 2))
        assert_array_equal(x, [[3, 6]])

    def test_M_fromlist(self):
        fh = DIR / "parse_test/test"
        x = ps.M_fromlist([fh, fh])
        assert_array_equal(x.shape, (1, 6))
        assert_array_equal(x, np.hstack((ps.M(fh), ps.M(fh))))


class Test_Fam(unittest.TestCase):

    def test_fam(self):
        fam = ps.PlinkFAMFile(DIR / "plink_test/plink.fam")
        self.assertEqual(fam.n, 5)
        correct = np.array(["per0", "per1", "per2", "per3", "per4"])
        assert_array_equal(fam.IDList.values.reshape((5,)), correct)

    def test_bad_filename(self):
        self.assertRaises(ValueError, ps.PlinkFAMFile, DIR / "plink_test/plink.bim")


class Test_Bim(unittest.TestCase):

    def test_bim(self):
        bim = ps.PlinkBIMFile(DIR / "plink_test/plink.bim")
        self.assertEqual(bim.n, 8)
        correct = np.array(
            ["rs_0", "rs_1", "rs_2", "rs_3", "rs_4", "rs_5", "rs_6", "rs_7"]
        )
        assert_array_equal(bim.IDList.values.reshape(8), correct)

    def test_bad_filename(self):
        self.assertRaises(ValueError, ps.PlinkBIMFile, DIR / "plink_test/plink.fam")
