#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import sys
import warnings
import pandas
import tempfile
from pathlib import Path
from pypipegraph import FileGeneratingJob, MultiFileGeneratingJob
from mvariants import VarScan, OptionHandler, Mutect2, parse_vcf2

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


def test_parse_vcf2():
    vcf_file = "/project/code/mvariants/data/snp_calling_expected/raw.snps.indels.vcf"
    df = parse_vcf2(vcf_file)
    assert type(df) == pandas.DataFrame


def test_parse_vcf2_raise(tmpdir):
    f = Path(tmpdir) / "somefile.txt"
    with f.open("w") as outp:
        outp.write("...")
    with pytest.raises(ValueError):
        parse_vcf2(f)
