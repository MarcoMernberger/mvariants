#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""post_process.py: Contains different postprocessor for variant calling."""

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

from pathlib import Path
from mbf.align.post_process import _PostProcessor
import subprocess


class Mutect2Processor(mbf.align.post_process._PostProcessor):
    def __init__(
        self,
        result_folder="Mutect2",
        libname="mutect2",
        barcode="NB552003",
        platform="illumina",
    ):
        self.name = f"Mutectprocessor"
        self.platform = platform
        self.barcode = barcode
        self.libname = libname
        self.result_folder_name = result_folder
        # TODO: finish implementation and uncouple from the variant call

    def process(self, input_bam_name, output_bam_name, result_dir):
        cmd = [
            "python",
            self.code_path,
            "AddOrReplaceReadGroups",
            "-I",
            f"{str(input_bam_name)}",
            "-O",
            f"{str(output_bam_name)}",
            "--RGLB",
            f"{self.libname}",
            "--RGPL",
            f"{self.platform}",
            "--RGPU",
            f"{self.barcode}",
            "--RGSM",
            f"{input_bam_name.stem}",
        ]
        with Path(str(output_bam_name) + ".stderr").open("wb") as stderr:
            subprocess.check_call(cmd, stderr=stderr)
        subprocess.check_call(cmd)

    def register_qc(self, new_lane):
        pass  # pragma: no cover

    def get_version(self):
        return (
            subprocess.check_output(["python", self.code_path, "--version"])
            .decode("utf-8")
            .strip()
        )

    def get_parameters(self):
        return (self.get_version(),)

    def filter_mutect(self, infile, output_file, dependencies=[]):
        if isinstance(output_file, str):
            outfile = Path(output_file)
        outfile.parent.mkdir(parents=True, exist_ok=True)

        def __dump():
            cmd = [
                "python",
                "code/gatk/gatk-4.1.4.0/gatk",
                "FilterMutectCalls",
                "-V",
                str(infile),
                "-R",
                str(self.gatk_genome_file),
                "-O",
                str(outfile),
            ]
            print(" ".join(cmd))
            with Path(str(outfile) + ".stderr").open("wb") as stderr:
                subprocess.check_call(cmd, stderr=stderr)

        return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)


class DuplicateMarker:
    """
    Marks duplicate reads in alignments.

    [extended_summary]
    """

    def __init__(self):
        # TODO : implement
        pass


class GATKBaseScoreRecalibration:
    """
    Wraps the GATK base score recalibration.

    [extended_summary]
    """

    # TODO : implement
