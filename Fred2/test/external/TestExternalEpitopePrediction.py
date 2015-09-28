"""
Unittest for external epitope prediction methods
"""

import unittest
import os

from Fred2.Core import Allele
from Fred2.Core import Peptide
from Fred2.Core import Transcript

from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.EpitopePrediction import AExternalEpitopePrediction
from Fred2.EpitopePrediction import NetMHC_3_4


#only for internal testing (test are run with NetMHC 3.4
class NetMHC_3_0(NetMHC_3_4):

    __version = "3.0"

    @property
    def version(self):
        return self.__version


class TestExternalEpitopePredictionClass(unittest.TestCase):

    def setUp(self):
        self.peptides_mhcI = [Peptide("SYFPEITHI"),Peptide("IHTIEPFYS")]
        self.peptides_mhcII = [Peptide("AAAAAASYFPEITHI"),Peptide("IHTIEPFYSAAAAAA")]
        self.mhcI = [Allele("HLA-B*15:01"),Allele("HLA-A*02:01")]
        self.mhcII = [Allele("HLA-DRB1*07:01"), Allele("HLA-DRB1*15:01")]
        self.transcript = Transcript("")

    def test_multiple_inputs(self):
        for m in EpitopePredictorFactory.available_methods():
            mo = EpitopePredictorFactory(m)
            if isinstance(mo, AExternalEpitopePrediction):
                if any(a.name in mo.supportedAlleles for a in self.mhcII):
                    mo.predict(self.peptides_mhcII, alleles=self.mhcII)
                else:
                    mo.predict(self.peptides_mhcI, alleles=self.mhcI)

    def test_single_epitope_input(self):
        for m in EpitopePredictorFactory.available_methods():
            mo = EpitopePredictorFactory(m)
            if isinstance(mo, AExternalEpitopePrediction):
                if any(a.name in mo.supportedAlleles for a in self.mhcII):
                    mo.predict(self.peptides_mhcII[0], alleles=self.mhcII)
                else:
                    mo.predict(self.peptides_mhcI[0], alleles=self.mhcI)

    def test_single_allele_input(self):
        for m in EpitopePredictorFactory.available_methods():
            mo = EpitopePredictorFactory(m)
            if isinstance(mo, AExternalEpitopePrediction):
                if any(a.name in mo.supportedAlleles for a in self.mhcII):
                    mo.predict(self.peptides_mhcII, alleles=self.mhcII[0])
                else:
                    mo.predict(self.peptides_mhcI, alleles=self.mhcI[0])

    def test_wrong_epitope_input(self):
        with self.assertRaises(ValueError):
            EpitopePredictorFactory("NetMHC").predict(self.transcript, alleles=self.mhcI)

    def test_wrong_allele_input(self):
        with self.assertRaises(ValueError):
            EpitopePredictorFactory("NetMHC").predict(self.mhcI, alleles=self.transcript)

    def test_wrong_internal_to_external_version(self):
        with self.assertRaises(RuntimeError):
            EpitopePredictorFactory("NetMHC", version="3.0").predict(self.peptides_mhcI, alleles=self.mhcI)

    def test_path_option_and_optional_parameters(self):
        netmhc = EpitopePredictorFactory("NetMHC")
        exe = netmhc.command.split()[0]
        for try_path in os.environ["PATH"].split(os.pathsep):
            try_path = try_path.strip('"')
            exe_try = os.path.join(try_path, exe).strip()
            if os.path.isfile(exe_try) and os.access(exe_try, os.X_OK):
                netmhc.predict(self.peptides_mhcI, alleles=self.mhcI, path=exe_try, options="--sort")

    def test_path_and_optional_parameters_netctl(self):
        netctlpan = EpitopePredictorFactory("NetCTLpan")
        exe = netctlpan.command.split()[0]
        for try_path in os.environ["PATH"].split(os.pathsep):
            try_path = try_path.strip('"')
            exe_try = os.path.join(try_path, exe).strip()
            if os.path.isfile(exe_try) and os.access(exe_try, os.X_OK):
                print netctlpan.predict(self.peptides_mhcI, alleles=self.mhcI,
                                        path=exe_try,
                                        options="-wt 0.05 -wc 0.225 -ethr 0.5")

if __name__ == '__main__':
    unittest.main()