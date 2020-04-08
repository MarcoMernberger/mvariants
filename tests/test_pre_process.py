#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import sys
import warnings
import mbf_align
import mbf_genomes
import pypipegraph as ppg
from mvariants import SamtoolsmPileupSingleFile, GATKPreprocess
from pathlib import Path

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


@pytest.fixture
def samtools_preprocessor(tmpdir):
    preprocessor = SamtoolsmPileupSingleFile(
        name="Test",
        sampling_depth=100000,
        min_base_quality=15,
        chrm="2",
        start_stop=(72933703, 72933949),
        flag_filter=["UNMAP"],
        count_orphans=False,
        region_file=None,
        cache_dir=tmpdir,
    )
    return preprocessor


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_get_dependencies(samtools_preprocessor):
    deps = samtools_preprocessor.get_dependencies()
    assert len(deps) > 0
    assert all([isinstance(dep, ppg.Job) for dep in deps])


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_get_filename(samtools_preprocessor, test_lanes):
    ppath = samtools_preprocessor._get_filename(test_lanes[0])
    expected = samtools_preprocessor.cache_dir / "base_raw_test_hg36_subread.mpileup"
    assert ppath == expected
    ppath = samtools_preprocessor._get_filename(test_lanes[1])
    expected = samtools_preprocessor.cache_dir / "base_raw_test_hg3612.mpileup"
    assert ppath == expected


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_get_preprocessed_output(samtools_preprocessor, test_lanes):
    outputfiles = samtools_preprocessor.get_preprocessed_output(test_lanes)
    expected = [
        samtools_preprocessor.cache_dir / "base_raw_test_hg36_subread.mpileup",
        samtools_preprocessor.cache_dir / "base_raw_test_hg3612.mpileup",
    ]
    assert outputfiles == expected


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_run_modifier(samtools_preprocessor, test_lanes):
    def dummy_caller_function(input_lanes, *args):
        return [input_lanes] + list(args)

    modifier = samtools_preprocessor.run_modifier()
    modified_dummy = modifier(dummy_caller_function)
    expected = [
        [
            samtools_preprocessor.cache_dir / "base_raw_test_hg36_subread.mpileup",
            samtools_preprocessor.cache_dir / "base_raw_test_hg3612.mpileup",
        ],
        "somearg",
    ]
    print(modified_dummy(test_lanes, "somearg"))
    print(expected)
    assert modified_dummy(test_lanes, "somearg") == expected


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_preprocess(samtools_preprocessor, test_lanes, tmpdir):
    do_something = samtools_preprocessor.preprocess(test_lanes)
    do_something()
    for filepath in samtools_preprocessor.get_preprocessed_output(test_lanes):
        assert filepath.exists()


def test_gatk_get_dependencies(gatk_preprocessor):
    deps = gatk_preprocessor.get_dependencies()
    if len(deps) > 0:
        assert all([isinstance(dep, ppg.Job) for dep in deps])


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_gatk_prerequisite_jobs(gatk_preprocessor):
    jobs = gatk_preprocessor.prerequisite_jobs()
    assert len(jobs) == 3
    assert all([isinstance(job, ppg.FileGeneratingJob) for job in jobs])


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_gatk_get_preprocessed_output(gatk_preprocessor, gatk_test_lanes):
    expected = [
        Path(str(sample.bam_filename) + ".checked.txt")
        for samples in gatk_test_lanes
        for sample in samples
    ]
    print(expected)
    print(gatk_preprocessor.get_preprocessed_output(gatk_test_lanes))
    assert gatk_preprocessor.get_preprocessed_output(gatk_test_lanes) == expected


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_gatk_preprocess(gatk_preprocessor, gatk_test_lanes):
    assert callable(gatk_preprocessor.preprocess(gatk_test_lanes))


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_gatk_run_modifier(gatk_preprocessor, gatk_test_lanes):
    def dummy_caller_function(*args):
        return list(args)

    modifier = gatk_preprocessor.run_modifier()
    modified_dummy = modifier(dummy_caller_function)
    expected = [
        gatk_test_lanes,
        "some_output_filename",
        gatk_preprocessor.gatk_compliant_genome_file,
        "somearg",
    ]
    assert (
        modified_dummy(gatk_test_lanes, "some_output_filename", "somearg") == expected
    )


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_gatk_prerequisite_jobs_creation(gatk_preprocessor, tmpdir):
    gatk_preprocessor.prerequisite_jobs()
    ppg.run_pipegraph()
    gatk_compliant_genome_file = gatk_preprocessor.cache_dir / "genome.fasta"
    assert gatk_compliant_genome_file.exists()
    assert Path(str(gatk_compliant_genome_file) + ".fai").exists()
    assert (
        gatk_compliant_genome_file.parent / (gatk_compliant_genome_file.stem + ".dict")
    ).exists()
