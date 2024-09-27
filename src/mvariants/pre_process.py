#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""pre_process.py: Contains methods for preprocessing before actual mutation analysis."""

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Callable, Tuple, Union, Any
from pathlib import Path
from mbf.align import AlignedSample
from mbf.align.post_process import _PostProcessor
from mbf.genomes import GenomeBase
from pypipegraph import Job
from .base import OptionHandler, GATK
from mbf.externals import ExternalAlgorithmStore
import mbf.align
import subprocess
import pypipegraph as ppg
import shutil
import sys


class _PreProcessor(ABC):
    """
    Abstract class for preprocessors.

    Takes care of preprocessing steps and modifies the function call of the
    caller. Every preprocessing that is needed for mutation calling should be
    handled by a preprocessor.
    """

    @abstractmethod
    def get_preprocessed_output(
        self, input_samples: List[List[AlignedSample]]
    ) -> List[Path]:
        """
        Returns a list of Path objects to be created from input_samples.

        Returns a list of Path objects that correspond to the files created by
        the preprocessor from input_samples. That way, the preprocessor can tell
            which files it is going to create via ~self.preprocess.

        Parameters
        ----------
        input_samples : List[List[mbf.align.AlignedSample]]
            List of two lists containing the input samples to be analyzed.
            Multiple samples are given in the first list. Additional matched
            samples can be supplied in the second list.

        Returns
        -------
        List[Path]
            List of Path objects to be created.
        """
        pass  # pragma: no cover

    @abstractmethod
    def preprocess(
        self, input_samples: List[List[AlignedSample]]
    ) -> Optional[Callable]:
        """
        Returns a function that creates input files for the variant caller.

        Returns a function that creates the input files based on input samples.
        VariantCalls use this method for defining FileGeneratingJobs on which
        variant calling can depend.

        Parameters
        ----------
        input_samples : List[List[mbf.align.AlignedSample]]
            List of two lists containing the input samples to be analyzed.
            Multiple samples are given in the first list. Additional Matched
            samples can be supplied in the second list.

        Returns
        -------
        Callable
            The function that actually creates the files.
        """
        pass

    @abstractmethod
    def run_modifier(self) -> Callable:
        """
        Returns a wrapper that modifies the variant caller run function.

        Returns a wrapper that modifies the variant caller run function, e.g.
        to supply different arguments that only the preprocessor knows, such
        as pileup files instead of input samples.

        Returns
        -------
        Callable
            A function that wraps Caller.run().
        """
        pass  # pragma: no cover

    def prerequisite_jobs(self) -> List[Job]:
        """
        Returns a list of global prerequisite jobs.

        Returns a list of global prerequisite jobs necessary
        for all samples to be analyzed.

        Returns
        -------
        List[ppg.Job]
            List of jobs that have to be done before the actual mutation calling.
        """
        return []

    def get_dependencies(self) -> List[Job]:
        """
        Returns a list of dependencies.

        Returns a list of pypipegraph.Job instances that
        need to run before the actual mutation analysis.

        Returns
        -------
        typing.List[ppg.Job]
            List of pypipegraph.Job instances.
        """
        return []


class SamtoolsmPileupSingleFile(_PreProcessor):
    """
    Pileup preprocessor, turns input bams into pileup files.

    This preprocessor is a wrapper for samtools and creates a single
    pileup file for a list of samples using the mpileup command. 
    If matched tumor-normal samples are given, it will create two pileup
    files, one is tumor the other is normal.
    No additional preprocessing steps are taken.
    """

    def __init__(
        self,
        sampling_depth: int = 1000000,
        min_base_quality: int = 15,
        chrm: Optional[str] = None,
        start_stop: Optional[Tuple[int]] = None,
        flag_filter: List = ["UNMAP", "SECONDARY", "QCFAIL", "DUP"],
        count_orphans: bool = False,
        region_file: Optional[Union[str, Path]] = None,
        *args,
        **kwargs,
    ) -> None:
        """
        SamtoolsmPileupSingleFile constructor.

        The constructor takes a couple of parameters/options that are
        used by samtools for pileup creation.

        Parameters
        ----------
        sampling_depth : int, optional
            Sampling depth of pileup per position, by default 1000000.
        min_base_quality : int, optional
            Minimum base quality to accept for pileup, by default 15.
        chrm : Optional[str], optional
            Chromosome of interest, by default None.
        start_stop : Optional[Tuple[int]], optional
            Genomic region of interest, by default None.
        flag_filter : List, optional
            List of samtools-specific bam flag filter indicating which reads
            to filter out, by default ["UNMAP", "SECONDARY", "QCFAIL", "DUP"].
        count_orphans : bool, optional
            Ignore orphan reads or count them, by default False.
        region_file : Optional[Union[str, Path]], optional
            File of genomic regions of interest, by default None.
        """
        self.instance_name = kwargs.get(
            "instance_name", "_".join(["Samtools", "mpileup"])
        )
        self.options = [
            "-t",
            "DP",
            "-Q",
            str(min_base_quality),
        ]
        if len(flag_filter) > 0:
            self.options.extend(["--ff", ",".join(flag_filter)])
        if count_orphans:
            self.options.append("-A")
        if chrm is not None:
            reg = str(chrm)
            if start_stop is not None:
                reg += f":{start_stop[0]}-{start_stop[1]}"
                self.options.extend(["-r", reg])
        if region_file is not None:
            self.options.extend(["-l", str(region_file)])
        if sampling_depth is not None:
            # samtools sampling depth defaults to 8000
            self.options.extend(["-d", str(sampling_depth)])
        self.options.extend(args)
        self.cache_dir = Path(
            kwargs.get("cache_dir", Path("cache") / self.instance_name)
        )
        if isinstance(self.cache_dir, str):
            self.cache_dir = Path(self.cache_dir)
        self.cache_dir.mkdir(exist_ok=True, parents=True)

    def get_dependencies(self) -> List[ppg.Job]:
        """
        Returns a list of dependencies. Overrides superclass method.

        Returns a list of pypipegraph.Job instances that need to run before the
        actual mutation analysis.

        Returns
        -------
        typing.List[ppg.Job]
            List of pypipegraph.Job instances.
        """
        return [
            ppg.ParameterInvariant(f"{self.instance_name}", list(self.options)),
            ppg.FunctionInvariant(self.instance_name + "_pre_process", self.preprocess),
            ppg.FunctionInvariant(
                self.instance_name + "_run_modifier", self.run_modifier
            ),
        ]

    def _get_filename(self, input_samples: List[AlignedSample]) -> Path:
        """
        Returns the filename of the pileup file to be created from input_samples.

        Returns the filename of the pileup file to be created from input_samples
        by self.preprocess.
                
        Parameters
        ----------
        input_samples : List[mbf.align.AlignedSample]
            List of one or multiple samples as input for the pileup file.
        
        Returns
        -------
        str
            File name of mpileup file to be created from input_samples.
        """
        input_name = f"{'.'.join([input_sample.bam_filename.stem for input_sample in input_samples])}.mpileup"
        if len(input_name) >= 250:
            input_name = input_name[:250]
        return self.cache_dir / input_name

    def get_preprocessed_output(
        self, input_samples: List[List[AlignedSample]]
    ) -> List[Path]:
        """
        Returns a list of Path objects to be created from input_samples.

        Returns a list of Path objects that correspond to the files created by
        the preprocessor from input_samples. That way, the preprocessor can tell
        which files it is going to create via self.preprocess(). Overrides
        abstract superclass method ~Preprocessor.get_preprocessed_output.

        Parameters
        ----------
        input_samples : List[List[mbf.align.AlignedSample]]
            List of two lists containing the input samples to be analyzed. 
            Multiple samples are given in the first list. Additional matched
            samples can be supplied in the second list.

        Returns
        -------
        List[Path]
            List of Path objects to be created.
        """

        samples_name = self._get_filename(input_samples[0])
        if len(input_samples[1]) > 0:
            return [samples_name, self._get_filename(input_samples[1])]
        return [samples_name]

    def run_modifier(self) -> Callable:
        """
        Returns a wrapper that modifies the variant caller run function.

        Returns a wrapper that modifies the variant caller run function by
        replacing the default input_samples argument with mpileup files.
        Overrides abstract superclass method ~PreProcessor.run_modifier.

        Returns
        -------
        Callable
            A function that wraps Caller.run().
        """

        def wrapper(caller_func):
            def inner(*args, **kwargs):
                new_args = [self.get_preprocessed_output(args[0])] + list(args[1:])
                return caller_func(*new_args, **kwargs)

            return inner

        return wrapper

    def __do_pileup_single_file(
        self, input_samples: List[AlignedSample], reference_path: Path
    ) -> None:
        """
        Creates a single pileup file from input_Samples.

        Invokes samtools to create a sinlge potentially multi-sample pileup
        file.

        Parameters
        ----------
        input_samples : List[mbf.align.AlignedSample]
            List of one or multiple samples as input for the pileup file.
        reference_path : Path
            Path to reference file (e.g. genome.fasta).

        Raises
        ------
        subprocess.CalledProcessError
            If subprocess call failed.
        """
        input_bams = [input_sample.get_bam_names()[0] for input_sample in input_samples]
        mpileup = self._get_filename(input_samples)
        # for each list of samples in input_samples we must pile up
        with Path(str(mpileup) + ".stderr").open("w") as stderr:
            cmd_generate_pileup = [
                "samtools",
                "mpileup",
                "-f",
                str(reference_path),
                "-o",
                str(mpileup),
            ]
            cmd_generate_pileup.extend(self.options)
            for input_bam in input_bams:
                cmd_generate_pileup.extend(["-I", input_bam])
            stderr.write(" ".join(cmd_generate_pileup) + "\n")
            stderr.flush()
            try:
                subprocess.check_call(cmd_generate_pileup, stderr=stderr)
            except subprocess.CalledProcessError:
                print(
                    f"Samtools pileup didn't work, Command was:{' '.join(cmd_generate_pileup)}"
                )
                raise

    def preprocess(self, input_samples: List[List[AlignedSample]]) -> Callable:
        """
        Returns a function that creates input files for the variant caller.

        Returns a function that creates a single pileup file for the first
        and optionally second list of input_samples. Overrides the abstract
        superclass method ~PreProcessor.preprocess.

        Parameters
        ----------
        input_samples : List[List[mbf.align.AlignedSample]]
            List of two lists containing the input samples to be analyzed.
            Multiple samples are given in the first list. Additional Matched
            samples can be supplied in the second list.

        Returns
        -------
        Callable
            A function that wraps Caller.run().
        """

        def do_pileup():
            reference_path = input_samples[0][0].genome.find_file("genome.fasta")
            for sample_list in input_samples:
                if len(sample_list) > 0:
                    self.__do_pileup_single_file(sample_list, reference_path)

        return do_pileup


class GATKPreprocessor(GATK, _PreProcessor):
    """
    GATK preprocessor, ensures that reference files are compliant to GATK
    best practices.

    This preprocessor is a wrapper for the GATK toolbox and creates a number
    of input files that GATK requires to work.
    """

    def __init__(
        self,
        genome: GenomeBase,
        instance_name: str = None,
        version: str = "_last_used",
        options: Dict[str, Any] = {},
        store=None,
        **kwargs,
    ):
        """
        GATKPreprocess constructor.

        Accepts a genome that needs to be preprocessed to ensure GATK compliance.

        Parameters
        ----------
        genome : GenomeBase
            The reference genome of the aligned samples.
        name : str, optional
            [description], by default None
        """
        super().__init__(
            tool="ValidateSamFile", options=options, version=version, store=store,
        )
        self.instance_name = (
            f"GATKPreprocessor_{genome.name}"
            if instance_name is None
            else instance_name
        )
        self.genome = genome
        self.cache_dir = Path(
            kwargs.get("cache_dir", Path("cache") / self.instance_name)
        )
        self.cache_dir.mkdir(exist_ok=True, parents=True)
        self.gatk_compliant_genome_file = self.cache_dir / "genome.fasta"

    def lane_postprocessor(self) -> _PostProcessor:
        return GATKLanePostProcessor()

    def create_post_processed_lane(
        self, input_sample: AlignedSample, **kwargs
    ) -> AlignedSample:
        return input_sample.post_process(
            self.lane_postprocessor(),
            new_name=kwargs.get("name", f"{input_sample.name}_{self.instance_name}"),
            result_dir=kwargs.get("result_dir", input_sample.result_dir),
        )

    def prerequisite_jobs(self) -> List[Job]:
        """
        Returns a list of global prerequisite jobs.

        Returns a list of global prerequisite jobs necessary for all samples:
        copying and reindexing the reference fasta, creating a dictionary file.

        Returns
        -------
        List[ppg.Job]
            List of jobs that have to be done before the actual mutation calling.
        """

        def copy_genome(gatk_compliant_genome_file):
            # GATK needs a dictionary in the same directory as the genome so we copy
            shutil.copy(
                self.genome.find_file("genome.fasta"), str(gatk_compliant_genome_file)
            )

        def create_index(gatk_compliant_genome_file):
            # create the index file
            cmd = ["samtools", "faidx", self.gatk_compliant_genome_file]
            with Path(str(gatk_compliant_genome_file) + ".fai.stderr").open(
                "wb"
            ) as stderr:
                subprocess.check_call(cmd, stdout=stderr)

        def create_dict():
            # create the dictionary
            arguments = [
                "CreateSequenceDictionary",
                "-R",
                str(self.gatk_compliant_genome_file),
            ]
            cmd = self.build_cmd(self.gatk_compliant_genome_file.parent, 1, arguments)
            with Path(str(self.gatk_compliant_genome_file) + ".stderr").open(
                "wb"
            ) as stderr:
                subprocess.check_call(cmd, stderr=stderr)

        job1 = ppg.FileGeneratingJob(self.gatk_compliant_genome_file, copy_genome)
        job2 = ppg.FileGeneratingJob(
            str(self.gatk_compliant_genome_file) + ".fai", create_index
        ).depends_on(job1)
        job3 = ppg.FileGeneratingJob(
            self.gatk_compliant_genome_file.parent
            / (self.gatk_compliant_genome_file.stem + ".dict"),
            create_dict,
        ).depends_on([job1, job2])
        return [job1, job2, job3]

    def _check_bam(self, bam_filename: Path) -> None:
        """
        Checks a bam file for compliance with GATK standards.

        This will complain if GATK cannot deal with your bam files.        

        Parameters
        ----------
        bam_filename : Path
            Bam file to be tested.
        """
        arguments = [
            "-I",
            str(bam_filename),
            "-R",
            str(self.gatk_compliant_genome_file),
            "-M",
            "SUMMARY",
            "-O",
            str(bam_filename) + ".checked.txt",
        ]
        cmd = self.build_cmd(None, 1, arguments)
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            print(subprocess.STDOUT)
            print(
                "Check of bam files failed. You may want to postprocess your lanes using GATKLanePostProcessor."
            )
            print("-".join(cmd))
            raise

    def preprocess(self, input_samples: List[List[AlignedSample]]) -> Callable:
        """
        Returns a function that creates input files for the variant caller.

        Creating input files for GATK includes adding read groups and checking
        the resulting files afterwards.

        Parameters
        ----------
        input_samples : List[List[mbf.align.AlignedSample]]
            List of two lists containing the input samples to be analyzed.
            Multiple samples are given in the first list. Additional Matched
            samples can be supplied in the second list.

        Returns
        -------
        Callable
            The function that actually creates the files.
        """

        def do_process():
            for samples in input_samples:
                for sample in samples:
                    self._check_bam(sample.bam_filename)

        return do_process

    def get_preprocessed_output(self, input_samples: List[List[AlignedSample]]) -> List:
        """
        Overrides abstract method in superclass.

        Returns
        -------
        List
            Empty list, since no files are created.
        """
        ret = [
            Path(str(sample.bam_filename) + ".checked.txt")
            for samples in input_samples
            for sample in samples
        ]
        return ret

    def run_modifier(self) -> Callable:
        """
        Returns a wrapper that modifies the variant caller run function.

        Overrides the abstract superclass method.

        Returns
        -------
        Callable
            A function that wraps Caller.run().
        """

        def wrapper(caller_func):
            def inner(*args, **kwargs):
                new_args = [
                    args[0],
                    args[1],
                    self.gatk_compliant_genome_file,
                    *args[2:],
                ]
                return caller_func(*new_args, **kwargs)

            return inner

        return wrapper

    def get_dependencies(self) -> List[ppg.Job]:
        """
        Returns a list of dependencies.

        Returns a list of pypipegraph.Job instances that
        need to run before the actual mutation analysis. Overrides the 
        superclass method.

        Returns
        -------
        typing.List[ppg.Job]
            List o  f pypipegraph.Job instances.
        """
        return [
            ppg.FunctionInvariant(self.instance_name + "_pre_process", self.preprocess),
            ppg.FunctionInvariant(
                self.instance_name + "_run_modifier", self.run_modifier
            ),
        ]


class GATKLanePostProcessor(GATK, _PostProcessor):
    """
    PostProcessor creating GATK-compliant bam files.

    PostProcessor for mbf.align.AlignedSample classes that creates GATK
    compliant bam files. Since this is a preprocessing step for variant calling
    with GATK it is placed here.
    Subclasses GATK and mbf.align.post_process._PostProcessor.

    Parameters
    ----------
    options : Optional[Dict[str, str]], optional
        Dictionary of parameter-values to be supplied to the GATK call.
    version : str, optional
        Version of the tool to be used, by default "_last_used".
    store : ExterrnalAlgorithmStore, optional
        Store that handles downloads and provides thh external tool, by default
        None.
    """

    def __init__(
        self,
        options: Dict[str, str] = {"-PL": "ILLUMINA", "-PU": "na", "-LB": "na"},
        version: str = "_last_used",
        store: Optional[ExternalAlgorithmStore] = None,
    ):
        """GATKLanePostProcessor constructor, see class documentation for details."""
        super().__init__(
            tool="AddOrReplaceReadGroups",
            options=options,
            version=version,
            store=store,
        )

    def process(
        self, input_bam_name: Path, output_bam_name: Path, result_dir: Path
    ) -> None:
        """
        Adds a read group to bam files and outputs them to a new bam file in 
        cache directory.

        Parameters
        ----------
        inputbam : Path
            Input bam file.
        outputbam : Path
            The resulting output bam file.
        sample_name : str
            Name of the sample.
        """
        cmd_arguments = [
            "-I",
            str(input_bam_name),
            "-O",
            str(output_bam_name),
            "-SM",
            str(output_bam_name.stem),
        ]
        cmd_arguments.extend(self.optionhandler.options_as_list_str(self.options))
        cmd = self.build_cmd(result_dir, 1, cmd_arguments)
        print(" ".join(cmd))
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError:
            print(" ".join(cmd))
            raise

    def register_qc(self, new_lane):
        pass  # pragma: no cover
