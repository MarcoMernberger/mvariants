#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""effect_prediction.py: Contains methods vor predicting the result of a mutation."""

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

from mbf.externals import ExternalAlgorithm, ExternalAlgorithmStore
from mbf.externals.util import download_file
from typing import Optional, Dict, List
from pathlib import Path
from pypipegraph import Job
from .variant_calls import VariantCall
import pypipegraph as ppg
import os
import time
import random
import docker


class VEP(ExternalAlgorithm):
    """
    VEP Wrapper for the ensembl variant effect predictor docker.

    This subclasses from ExternalAlgorithm to take care of the cache download.
    In addition, it retrieves the ensembl-vep docker image and attaches the
    docker container to the project container.

    Parameters
    ----------
    species : str, optional
        The species name, by default "homo_sapiens". Should be the same as
        the genome species.
    revision : int, optional
        Ensembl revision to be used, by default 99.
    grch : str, optional
        The assembly used, by default "GRCh38".
    store : ExterrnalAlgorithmStore, optional
        Store that handles downloads and provides thh external tool, by default
        None.
    """

    def __init__(
        self,
        species: Optional[str] = "homo_sapiens",
        revision: Optional[int] = 99,
        grch: Optional[str] = "GRCh38",
        store: Optional[ExternalAlgorithm] = None,
    ) -> None:
        """VEP constructor, see class documentation for details."""
        self.grch = grch
        self.revision = revision
        self.species = species.lower()
        self.full_version_cache = f"{self.species}_vep_{self.revision}_{self.grch}"
        version = f"{self.species}_{self.revision}"
        super().__init__(version, store)
        self.image = "ensemblorg/ensembl-vep"
        self.wdir = "/project"
        self.volumes = {
            os.environ["ANYSNAKE_PROJECT_PATH"]: {"bind": "/project", "mode": "rw"},
        }
        self.client = docker.from_env()

    latest_version = "homo_sapiens_99"

    def fetch_version(self, version: str, target_filename: Path) -> None:
        """
        Takes care of the tool download.

        Overrides the ExternalAlgorithm methood. Downloads the VEP cache
        to the prebuild location specified by the corresponding
        ExternalAlgorithmStore and packs it into a tar.gz file.

        Parameters
        ----------
        version : str
            The tool version to be used.
        target_filename : Path
            The path to the local tar.gz file.
        """
        url = f"ftp://ftp.ensembl.org/pub/release-{self.revision}/variation/indexed_vep_cache/{self.species}_vep_{self.revision}_{self.grch}.tar.gz"
        download_file(url, target_filename.open("wb"))

    @property
    def name(self) -> str:
        """
        Getter for the name of the external method for version handling.

        Overrides the ExternalAlgorithm method.

        Returns
        -------
        str
            Name of the external method.
        """

        return "VEP"

    @property
    def multi_core(self) -> bool:
        """
        Returns wehther to use multiple cores.

        Overrides the ExternalAlgorithm method.

        Returns
        -------
        bool
            Whether the external method can use multiple cores.
        """
        return False

    def get_latest_version(self):
        """Getter for the latest_version attribute."""
        return self.latest_version

    def get_version(self) -> str:
        """
        Returns the version of the tool.

        Retrieves the version from the help string.

        Returns
        -------
        str
            The version.
        """
        s = self.get_help()
        version = s[s.find("ensembl-vep") : s.find("Help")].split(":")[1].strip()
        return version

    def get_help(self) -> str:
        """
        Returns a help string from ther VEP tool.

        Runs the docker container with a help command and returns the output.

        Returns
        -------
        str
            Help message string.
        """
        command_help = ["vep", "--help"]
        container = self.client.containers.run(
            image=self.image,
            volumes=self.volumes,
            working_dir=self.wdir,
            command=command_help,
            detach=True,
        )
        ret = ""
        for line in container.logs(stream=True):
            ret += line.decode("utf-8")

        return ret

    def print_help(self) -> None:
        """Prints the help string."""
        print(self.get_help())

    def vep_run(
        self,
        variant_call: VariantCall,
        dependencies: List[Job] = [],
        options: Dict[str, str] = {},
        **kwargs,
    ):
        """
        Returns a Job that runs the variant effect prediction analysis.

        Returns a job that runs the VEP docker for a given variant call and
        creates an output file next to the variant call output file.

        Parameters
        ----------
        variant_call : VariantCall
            Variant call for which the prediction should be done.
        dependencies : list, optional
            List of pypipegraph.Job containing dependencies, by default [].
        options : dict, optional
            Additional options to be passed to VEP, by default {}.

        Returns
        -------
        pypipegraph.FileGeneratingJob
            Job that runs the analysis.
        """
        outputfolder = kwargs.get("result_dir", variant_call.result_dir)
        if isinstance(outputfolder, str):
            outputfolder = Path(outputfolder)
        output_file = f"{variant_call.output_file.name}.vep.txt"
        outfile = outputfolder / output_file
        outfile.parent.mkdir(parents=True, exist_ok=True)
        sentinel = Path(str(outfile) + ".log")
        deps = dependencies
        deps.append(
            ppg.ParameterInvariant(
                f"PI_{output_file}",
                [
                    self.volumes,
                    self.wdir,
                    str(options),
                ],
            )
        )
        deps.append(variant_call.load())

        def __dump():
            command = [
                "vep",
                "--database",
                "--dir",
                str(self.path),
                "-i",
                f"{str(variant_call.output_file)}",
                "-o",
                str(outfile),
                "--hgvs",
                "--hgvsg",
                "--symbol",
                "--everything",
                "--force_overwrite",
            ]
            for k in options:
                command.extend([k, str(options[k])])
            time.sleep(random.random() * 2)
            print(" ".join(command))
            container = self.client.containers.run(
                self.image,
                volumes=self.volumes,
                working_dir=self.wdir,
                command=command,
                detach=True,
            )
            with sentinel.open("w") as outp:
                outp.write(" ".join(command))
                for line in container.logs(stream=True):
                    line = line.strip().decode("utf-8")
                    outp.write(line)

        job = ppg.FileGeneratingJob(sentinel, __dump).depends_on(deps)
        return job

    def build_cmd(
        self, output_directory: Optional[Path], ncores: int, arguments: List[str]
    ):
        return []
