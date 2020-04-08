#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import sys
import warnings
import pandas as pd
import mbf_genomes

from pypipegraph import FileGeneratingJob, MultiFileGeneratingJob
from mvariants import VarScan, OptionHandler, Mutect2
from pathlib import Path

__author__ = "MarcoMernberger"
__copyright__ = "MarcoMernberger"
__license__ = "mit"


@pytest.fixture
def varscan_default():
    caller_snp = VarScan("mpileup2snp")
    default_arguments = {
        "--min-coverage": 8,
        "--min-reads2": 2,
        "--min-avg-qual": 15,
        "--min-var-freq": 0.01,
        "--min-freq-for-hom": 0.75,
        "--p-value": 99e-02,
        "--strand-filter": 1,
        "--output-vcf": 1,
        "--variants": 0,
    }
    help_str, allowed_options, default_values = caller_snp._init_options()
    return caller_snp, default_arguments, help_str, allowed_options, default_values


@pytest.fixture
def varscan_types():
    caller_snp = VarScan("mpileup2snp")
    caller_indel = VarScan("mpileup2indel")
    caller_somatic = VarScan("somatic")
    caller_cns = VarScan("mpileup2cns")
    caller_copynumber = VarScan("copynumber")

    return caller_snp, caller_indel, caller_somatic, caller_cns, caller_copynumber


def test_varscan_help_str(varscan_default):
    print(varscan_default[2])
    assert "USAGE" in varscan_default[2]


def test_varscan_allowed_options(varscan_default):
    assert type(varscan_default[3]) == dict


def test_varscan_default_param(varscan_default):
    assert type(varscan_default[4]) == dict


def test_varscan_allowed_options_keys(varscan_default):
    assert all([param in varscan_default[3] for param in varscan_default[1]])


def test_varscan_default_keys(varscan_default):
    assert all([param in varscan_default[4] for param in varscan_default[1]])


def test_option_handler_accepted_arguments(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    check = []
    for option in varscan_default[1]:
        check.append(option in vs_option_handler.accepted_arguments)
    assert all(check)


def test_option_handler_default_params(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    check = []
    for option in varscan_default[1]:
        check.append(option in vs_option_handler.default_parameter)
    assert all(check)


def test_option_handler_default_values(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    check = []
    for option in varscan_default[1]:
        check.append(
            varscan_default[1][option] == vs_option_handler.default_parameter[option]
        )
    assert all(check)


def test_option_handler_accepted_arguments_str(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    assert type(vs_option_handler.accepted_arguments_str()) == str


def test_option_handler_help_str(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    assert type(vs_option_handler._help_str) == str


def test_option_handler_option_list(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    options_list = OptionHandler.options_as_list(vs_option_handler.default_parameter)
    for key, value in vs_option_handler.default_parameter.items():
        options_list.remove(key)
        if value != "":
            options_list.remove(value)
    assert len(options_list) == 0


def test_option_handler_check(varscan_default):
    vs_option_handler = varscan_default[0].option_handler
    options = {"--variants": 0}
    ret = vs_option_handler.check_options(options)
    assert type(ret) == dict
    assert "--variants" in ret
    options = {"--variantsss": "", "--output-vcf": 1}
    with pytest.warns(UserWarning):
        checked = vs_option_handler.check_options(options)
        assert "--variantsss" not in checked
        assert "--output-vcf" in checked
        assert checked["--output-vcf"] == 1
    pytest.warns(UserWarning, lambda: vs_option_handler.check_options(options))
    vs_option_handler.warn = "error"
    with pytest.raises(UserWarning):
        vs_option_handler.check_options(options)
    options = {"--p-value": 0, "--output-vcf": 0}
    checked = vs_option_handler.check_options(options)
    assert checked["--p-value"] == 0
    assert checked["--output-vcf"] == 0


def test_varscan_init_wrong():
    with pytest.raises(ValueError):
        VarScan("this is an incorrect variant type")


def test_varscan_snp_somatic_options(varscan_types):
    options_snp = set(OptionHandler.options_as_list(varscan_types[0].options))
    options_somatic = set(OptionHandler.options_as_list(varscan_types[2].options))
    diff = options_somatic.difference(options_snp)
    somatic_specific = [
        "--min-coverage-normal",
        "--min-coverage-tumor",
        "--normal-purity",
        "--tumor-purity",
        "--somatic-p-value",
    ]
    assert all([x in diff for x in somatic_specific])


def test_mutect2_result_to_df_parser(mutect2_default):
    parser = mutect2_default[0].result_to_df_parser()
    df = parser(
        "/project/code/mvariants/data/snp_calling_expected/Varscan_somatic/base_raw_test_hg36_Subread_vs_base_raw_test_hg3612_Subread_by_Varscan_somatic.snp.vcf"
    )
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_mutect_help_str(mutect2_default):
    assert "USAGE" in mutect2_default[1]


def test_mutect_allowed_options(mutect2_default):
    assert type(mutect2_default[2]) == dict


def test_mutect_default_param(mutect2_default):
    assert type(mutect2_default[3]) == dict


def test_mutect_allowed_options_keys(mutect2_default):
    assert all([param in mutect2_default[1] for param in mutect2_default[2]])


def test_mutect_default_keys(mutect2_default):
    assert all([param in mutect2_default[1] for param in mutect2_default[3]])


def test_mutect_option_handler_accepted_arguments(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    check = []
    for option in mutect2_default[2]:
        check.append(option in option_handler.accepted_arguments)
    assert all(check)


def test_mutect_option_handler_default_params(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    check = []
    for option in mutect2_default[3]:
        check.append(option in option_handler.default_parameter)
    assert all(check)


def test_mutect_option_handler_default_values(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    check = []
    for option in mutect2_default[3]:
        check.append(
            mutect2_default[3][option] == option_handler.default_parameter[option]
        )
    assert all(check)


def test_mutect_option_handler_accepted_arguments_str(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    assert type(option_handler.accepted_arguments_str()) == str


def test_mutect_option_handler_help_str(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    assert type(option_handler._help_str) == str


def test_mutect_option_handler_option_list(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    options_list = OptionHandler.options_as_list(option_handler.default_parameter)
    for key, value in option_handler.default_parameter.items():
        if value != "":
            if isinstance(value, list):
                for item in value:
                    print(options_list)
                    print(key)
                    options_list.remove(key)
                    options_list.remove(item)
            else:
                options_list.remove(key)
                options_list.remove(value)
    assert len(options_list) == 0


def test_mutect_option_handler_check(mutect2_default):
    option_handler = mutect2_default[0].option_handler
    options = {"-add-output-vcf-command-line": True}
    ret = option_handler.check_options(options)
    assert type(ret) == dict
    assert "-add-output-vcf-command-line" in ret
    options = {"--variantsss": ""}
    with pytest.warns(UserWarning):
        checked = option_handler.check_options(options)
        assert "--variantsss" not in checked
    pytest.warns(UserWarning, lambda: option_handler.check_options(options))
    option_handler.warn = "error"
    with pytest.raises(UserWarning):
        option_handler.check_options(options)


def test_mutect_init():
    with pytest.raises(ValueError):
        Mutect2("this has no genome")
    genome_human = mbf_genomes.EnsemblGenome("Homo_sapiens", 96)
    Mutect2(genome=genome_human)


def test_mutect2_init_options(mutect2_default):
    help_str, allowed_options, default_values = mutect2_default[0]._init_options()
    assert isinstance(help_str, str)
    assert isinstance(allowed_options, dict)
    assert isinstance(default_values, dict)


def test_mutect2_filter(mutect2_default):
    cmd = mutect2_default[0].filter_mutect(Path("test_in"), Path("gatk"), Path("out"))
    assert isinstance(cmd, list)
    assert all(isinstance(x, str) for x in cmd)
    assert cmd[1:] == ["FilterMutectCalls", "-V", "test_in", "-R", "gatk", "-O", "out"]

