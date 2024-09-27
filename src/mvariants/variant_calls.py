#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""variant_calling.py: This contains the result class for any variant caller."""

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

import pypipegraph2 as ppg
from typing import List, Union, Any, Callable
from mbf.align.lanes import AlignedSample
from pandas import DataFrame
from .caller import _Caller
from pypipegraph2 import CachedAttributeLoadingJob, FileGeneratingJob, Job
from .util import parse_vcf2
from pathlib import Path


class VariantCall:
    """
     Class to wrap the results of a variant caller.

    A VariantCall class is instantiated with a variant caller specifying the
    analysis and a list specifying the input samples. input samples are
    secified by a list of two lists which correspond to matched tumor-normal
    samples. Each of these lists may contain multiple samples to be analyzed
    together by the variant caller as a joint variant calling.

    Example
    -------
    # single sample call, no matched normal
    >>> VariantCall([[sample1], []], caller)
    # matched tumor normal sample call
    >>> VariantCall([[sample1_tumor], [sample1_normal]], caller)
    # joint variant calling, no matched normal
    >>> VariantCall([[sample1, sample2], []], caller)
    # joint variant calling, matched-normal
    >>> VariantCall([[sample1_tumor, sample2_tummor], [sample1_normal, sample2_normal]], caller)

    Parameters
    ----------
    input_samples : List[List[AlignedSample]]
        List of lists of matched tumor-normal samples. The second list may be empty.
    caller : Caller
        Variant caller to use.
    """

    def __init__(
        self,
        input_samples: List[List[AlignedSample]],
        caller: _Caller,
        dependencies: List[Job] = [],
        **wargs,
    ):
        """VariantCall constructor, see class documentation for details."""
        self.name = wargs.get("name", None)
        self.input_samples = input_samples
        if self.name is None:
            self.name = "-".join([sample.name for sample in input_samples[0]])
            if len(input_samples[1]) > 0:
                self.name += "_vs_" + "-".join([sample.name for sample in input_samples[1]])
            self.name += f"_by_{caller.name}"
        self.result_dir = Path(
            wargs.get("result_dir", input_samples[0][0].result_dir / "SNP_calling")
        )
        self.cache_dir = Path(wargs.get("cache_dir", Path("cache_dir") / self.name))
        self.filename = wargs.get("filename", f"{self.name}.vcf")
        ppg.util.assert_uniqueness_of_object(self)
        self.caller = caller
        self.result_dir = self.result_dir / self.caller.name
        # self.reference_file = input_samples[0][0].genome.find_file("genome.fasta")
        self.result_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.output_file = self.result_dir / self.filename
        self.vid = [i.vid for tup in input_samples for i in tup]
        self.dependencies = dependencies

    @classmethod
    def fromvcf(cls, vcf_file: Union[str, Path]) -> DataFrame:
        """
        fromvcf returns a DataFrame from vcf file.

        Convenience method to read a vcf file into a multi-indexed DataFrame.

        Parameters
        ----------
        vcf_file : Any[str, Path]
            vcf file to read.

        Returns
        -------
        pandas.DataFrame: DataFrame
            Multi-indexed DataFrame containing the result of a variant call.
        """
        reader_function = parse_vcf2
        df = reader_function(vcf_file)
        return df

    def get_vcf_reader_function(self) -> Callable:
        """
        Returns a parser for the result file to DataFrame.

        Returns the result file reader function from the corresponding _Caller
        class.

        Returns
        -------
        Callable
            Parser for the result file that returns a DataFrame.
        """
        if self.caller is not None:
            reader_function = self.caller.result_to_df_parser()
        else:
            reader_function = parse_vcf2
        return reader_function

    def load(self) -> CachedAttributeLoadingJob:
        """
        Ensures that a DataFrame from the result filer is loaded.

        Returns a load job that reads the result into a DataFrame which is then
        set as a class attribute df. Any downstream pypipegraph Jobs based
        on the variant call can depend on this.
        on this

        Returns
        -------
        pypipegraph.CachedAttributeLoadingJob
            Job dependency for downstream pipegraph jobs.
        """

        def do_load():
            reader_function = self.caller.result_to_df_parser()
            df = reader_function(str(self.output_file))
            return df

        return ppg.CachedAttributeLoadingJob(
            self.cache_dir / f"{self.name}_df_load", self, "df", do_load
        ).depends_on(self.call())

    def get_dependencies(self):
        return self.dependencies

    def call(self):
        """
        Creates the vcf generating job.

        Creates a pypipegraph Job that does the variant calling and takes care
        of preprocessing and dependencies.

        Returns
        -------
        pypipegraph.FileGeneratingJob
            The job that does the variant calling.
        """
        run_callable = self.caller.run()

        def run(output_file):
            run_callable(self.input_samples, output_file)

        job = ppg.FileGeneratingJob(self.output_file, run, empty_ok=False)
        # job can depend on preprocessor dependencies, caller dependencies and the preprocessor job
        job.depends_on(self.caller.get_dependencies())
        lanes_loaded = [
            input_sample.load()
            for sample_list in self.input_samples
            for input_sample in sample_list
        ]
        # If a snp caller needs some preparation, this is the place to do it.
        if self.caller.preprocessor is not None:
            job.depends_on(self.caller.preprocessor.get_dependencies())
            job.depends_on(self.caller.preprocessor.prerequisite_jobs())
            preprocessor_output = self.caller.preprocessor.get_preprocessed_output(
                self.input_samples
            )
            if len(preprocessor_output) > 0:
                preprocessor_job = ppg.MultiFileGeneratingJob(
                    preprocessor_output,
                    self.caller.preprocessor.preprocess(self.input_samples),
                )
                preprocessor_job.depends_on(lanes_loaded)
                preprocessor_job.depends_on(self.caller.preprocessor.get_dependencies()).depends_on(
                    self.caller.preprocessor.prerequisite_jobs()
                )
                job.depends_on(preprocessor_job)
        job.depends_on(self.get_dependencies())
        job.depends_on(lanes_loaded)
        job.depends_on(
            ppg.FunctionInvariant(
                f"{self.caller.__class__.__name__}.run",
                self.caller.__class__.run,
            )
        )
        for sample_list in self.input_samples:
            for input_sample in sample_list:
                job.depends_on(input_sample.load())
        job.cores_needed = self.caller.get_cores_needed()
        return job

    def write(self, output_filename: Union[str, Path]) -> FileGeneratingJob:
        """
        Returns a job that writes the result dataframe to a .tsv file.

        Returns
        -------
        pypipegraph.FileGeneratingJob
            Job that writes the result dataframe.
        """

        def __write(output_filename):
            self.df.to_csv(output_filename, sep="\t", index=False)

        return ppg.FileGeneratingJob(output_filename, __write).depends_on(self.load())
