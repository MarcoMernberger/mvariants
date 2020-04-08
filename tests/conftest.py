#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for mvariants.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import sys
import subprocess
import pathlib
import pytest
import mbf_genomes
import mbf_align

from pypipegraph.testing.fixtures import (  # noqa:F401
    new_pipegraph,
    pytest_runtest_makereport,
)
from mbf_qualitycontrol.testing.fixtures import new_pipegraph_no_qc  # noqa:F401
from pypipegraph.testing import force_load

root = pathlib.Path(__file__).parent.parent
print("root", root)
sys.path.append(str(root / "src"))
print("the path is", sys.path)
subprocess.check_call(["python3", "setup.py", "build_ext", "-i"], cwd=root)

from mvariants import GATKPreprocess, Mutect2


@pytest.mark.usefixtures("new_pipegraph_no_qc")
@pytest.fixture
def test_lanes():
    genome_human = mbf_genomes.EnsemblGenome("Homo_sapiens", 96)
    input_samples = [
        [
            mbf_align.AlignedSample(
                "Test1",
                "/project/code/mvariants/data/base_raw_test_hg36_subread.bam",
                genome_human,
                is_paired=False,
                vid=None,
            )
        ],
        [
            mbf_align.AlignedSample(
                "Test2",
                "/project/code/mvariants/data/base_raw_test_hg3612.bam",
                genome_human,
                is_paired=False,
                vid=None,
            )
        ],
    ]
    return input_samples


@pytest.mark.usefixtures("new_pipegraph_no_qc")
@pytest.fixture
def gatk_test_lanes():
    genome_human = mbf_genomes.EnsemblGenome("Homo_sapiens", 96)
    input_samples = [
        [
            mbf_align.AlignedSample(
                "Test1GATK",
                "/project/code/mvariants/data/base_raw_test_hg36_Subread_gatk_rg.bam",
                genome_human,
                is_paired=False,
                vid=None,
            )
        ],
        [
            mbf_align.AlignedSample(
                "Test2GATK",
                "data/base_raw_test_hg3612_Subread_gatk_rg.bam",
                genome_human,
                is_paired=False,
                vid=None,
            )
        ],
    ]
    return input_samples


@pytest.mark.usefixtures("new_pipegraph_no_qc")
@pytest.fixture
def gatk_preprocessor(tmpdir):
    genome_human = mbf_genomes.EnsemblGenome("Homo_sapiens", 96)
    preprocessor = GATKPreprocess(
        genome=genome_human, name="TestGATK_Preprocessor", cache_dir=tmpdir
    )
    return preprocessor


@pytest.fixture
def mutect2_default(gatk_preprocessor):
    mutect = Mutect2(preprocessor=gatk_preprocessor)
    help_str, allowed_options, default_values = mutect._init_options()
    return mutect, help_str, allowed_options, default_values
