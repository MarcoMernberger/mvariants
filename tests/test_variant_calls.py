import pytest
import sys
import warnings
import mbf_align
import mbf_genomes
import pandas
import pypipegraph as ppg
from pypipegraph import FileGeneratingJob, MultiFileGeneratingJob
from mvariants import (
    VarScan,
    VariantCall,
    SamtoolsmPileupSingleFile,
    parse_vcf2,
    GATKPreprocess,
    Mutect2,
)
from pathlib import Path


@pytest.fixture
def fix_samtools(tmpdir):
    preprocessor = SamtoolsmPileupSingleFile(
        name="AmpliconPileUp",
        sampling_depth=100000,
        min_base_quality=15,
        chrm="2",
        start_stop=(72933703, 72933949),
        flag_filter=["UNMAP"],
        count_orphans=False,
        region_file=None,
        cache_dir=str(tmpdir),
    )
    return preprocessor


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_vc_single(test_lanes, tmpdir, fix_samtools):
    input_samples = [test_lanes[0], []]
    expected_vcf = (
        "/project/code/mvariants/data/snp_calling_expected/raw.snps.indels.vcf"
    )
    preprocessor = fix_samtools
    options = {
        "--p-value": 0.000001,
        "--strand-filter": 0,
        "--min-var-freq": 0,
        "--min-coverage": 10,
        "--min-reads2": 5,
    }
    caller_snp = VarScan2("mpileup2snp", preprocessor=preprocessor, options=options)
    vc = VariantCall(input_samples, caller_snp, result_dir=str(tmpdir))
    vc.load()
    assert not Path(vc.output_file).exists()
    ppg.run_pipegraph()
    assert Path(vc.output_file).exists()
    assert Path(str(vc.output_file) + ".varscan.log").exists()
    reader = caller_snp.result_to_df_parser()
    df_produced = reader(vc.output_file)
    df_expected = reader(expected_vcf)
    pandas.testing.assert_frame_equal(df_expected, df_produced)

"""
@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_vc_somatic(test_lanes, tmpdir, fix_samtools):
    input_samples = test_lanes
    expected_indel = "/project/code/mvariants/data/snp_calling_expected/Varscan_somatic/base_raw_test_hg36_Subread_vs_base_raw_test_hg3612_Subread_by_Varscan_somatic.indel.vcf"
    expected_snp = "/project/code/mvariants/data/snp_calling_expected/Varscan_somatic/base_raw_test_hg36_Subread_vs_base_raw_test_hg3612_Subread_by_Varscan_somatic.snp.vcf"
    preprocessor = fix_samtools
    options = {
        "-A": "LikelihoodRankSumTest",
        "--strand-filter": 0,
        "--min-var-freq": 0,
        "--min-coverage": 10,
    }
    caller_soma = VarScan2("somatic", preprocessor=preprocessor, options=options)
    vc = VariantCall(input_samples, caller_soma, result_dir=str(tmpdir))
    print(caller_soma.options)
    assert not Path(vc.output_file).exists()
    vc.call()
    ppg.run_pipegraph()
    assert Path(vc.output_file).exists()
    for x in vc.output_file.parent.iterdir():
        print(x)
    res1 = vc.output_file.parent / (vc.output_file.stem + ".snp.vcf")
    assert res1.exists()
    res2 = vc.output_file.parent / (vc.output_file.stem + ".indel.vcf")
    assert res2.exists()
    assert Path(str(vc.output_file) + ".varscan.log").exists()
    df_expected = parse_vcf2(expected_snp, expected_indel)
    reader = caller_soma.result_to_df_parser()
    df_produced = reader(vc.output_file)
    pandas.testing.assert_frame_equal(df_expected, df_produced)


@pytest.mark.usefixtures("new_pipegraph_no_qc")
def test_gatk_mutect(gatk_test_lanes, tmpdir, gatk_preprocessor, mutect2_default):
    expected_vcf = "/project/code/mvariants/data/snp_calling_expected/mutect2/base_raw_test_hg36_Subread_gatk_rg_vs_base_raw_test_hg3612_Subread_gatk_rg_by_Mutect2.vcf"
    genome_human = mbf_genomes.EnsemblGenome("Homo_sapiens", 96)
    caller_mutect = Mutect2(
        name="Mutect2",
        preprocessor=GATKPreprocess(genome=genome_human, name="TestGATK_Preprocessor"),
        options={"-L": "2:72933703-72933949"},
    )
    vc = VariantCall(gatk_test_lanes, caller_mutect, result_dir=str(tmpdir))
    vc.load()
    assert not Path(vc.output_file).exists()
    ppg.run_pipegraph()
    assert Path(vc.output_file).exists()
    assert Path(str(vc.output_file) + ".mutect2.log").exists()
    reader = caller_mutect.result_to_df_parser()
    df_produced = reader(vc.output_file)
    df_expected = reader(expected_vcf)
    pandas.testing.assert_frame_equal(df_expected, df_produced)
"""