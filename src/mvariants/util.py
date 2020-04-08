#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains some utility functions."""

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

from pandas import DataFrame
from pathlib import Path
from typing import Union
import pandas as pd
import tempfile
import tarfile


def download_file_and_targz(url: str, filename: str, tar_filename: str) -> None:
    """
    Downloads a file and puts it into a tar.gz archive.

    Downloads a file then puts it into a tar.gz file, which may be needed by 
    ExternalAlgorithm.

    Parameters
    ----------
    url : str
        Download location.
    filename : str
        Name of the file inside the tar.gz.
    tar_filename : str
        Path to the tar.gz archive.

    Raises
    ------
    ValueError
        Raises if tar file does not have the correct suffix.
    """
    tar_filename = str(tar_filename)
    if not tar_filename.endswith(".tar.gz"):  # pragma: no cover
        raise ValueError(
            f"output filename did not end with .tar.gz, was {tar_filename}."
        )

    with tarfile.open(tar_filename, "w:gz") as tarout:
        with tempfile.NamedTemporaryFile(mode="wb") as tf:
            download_file(url, tf)
            tf.seek(0)
            tarout.add(tf.name, arcname=filename)


def parse_vcf2(input_vcf: Union[str, Path], *args: str) -> DataFrame:
    """
    Parses a vcf file and returns a DataFrame containing the variant calls.

    This returns a pd.DataFrame containing all default columns as well as
    optional info field values split into single columns and converted to a
    multilevel index.

    Parameters
    ----------
    input_vcf : str
        Full path to the vcf inpout file.

    Returns
    -------
    pd.DataFrame
        A pandas dataframe with mutlilevel index.
    """

    def join_split_info(df):
        df_sub = df.INFO.str.split(";", expand=True)
        df_sub.columns = [x[: x.find("=")] for x in df_sub.loc[0]]
        df_sub = df_sub.apply(lambda series: [x[x.find("=") + 1 :] for x in series])
        index = pd.MultiIndex.from_product(
            [["INFO"], df_sub.columns], names=["first", "second"]
        )
        df_sub.columns = index
        del df[("INFO", (""))]
        df = df.join(df_sub)
        return df

    def split_samples(df):
        fixed = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT", "INFO"]
        samples = [x for x in df.columns.get_level_values("first") if x not in fixed]
        if len(samples) > 0:
            # optional genotype fields present
            second = df.FORMAT[0].split(":")
            for sample in samples:
                df_sub = df[sample].str.split(":", expand=True)
                index = pd.MultiIndex.from_product(
                    [[sample], second], names=["first", "second"]
                )
                df_sub.columns = index
                del df[(sample, (""))]
                df = df.join(df_sub)
        return df

    dfs = []
    input_vcfs = [Path(input_vcf)]
    if len(args) > 0:
        input_vcfs += [Path(x) for x in args]
    for input_vcf in input_vcfs:
        comments = []
        to_df = []
        print(input_vcf)
        with input_vcf.open("r") as inp:
            # unfortunately, pandas does not take multiple char comment indicators
            for line in inp.readlines():
                if line.startswith("#"):
                    comments.append(line)
                else:
                    to_df.append(line[:-1].split("\t"))
        if len(comments) == 0:
            raise ValueError(f"{input_vcf} is not in valid vcf format.")

        columns = comments[-1][1:-1].split("\t")
        mcolumns = pd.MultiIndex.from_product(
            [columns, [""]], names=["first", "second"]
        )
        df = pd.DataFrame(to_df, columns=mcolumns)
        # split INFO
        df = join_split_info(df)
        # split samples
        df = split_samples(df)
        dfs.append(df)
    df = pd.concat(dfs)
    return df


def download_file(url, file_object):
    """Download an url"""
    if isinstance(file_object, (str, Path)):
        raise ValueError("download_file needs a file-object not a name")

    try:
        if url.startswith("ftp"):
            return download_ftp(url, file_object)
        else:
            return download_http(url, file_object)
    except Exception as e:
        raise ValueError("Could not download %s, exception: %s" % (repr(url), e))


def download_http(url, file_object):
    """Download a file from http"""
    import requests
    import shutil

    r = requests.get(url, stream=True)
    if r.status_code != 200:
        raise ValueError("HTTP Error return: %i fetching %s" % (r.status_code, url))
    r.raw.decode_content = True
    shutil.copyfileobj(r.raw, file_object)


def download_ftp(url, file_object):
    """Download a file from ftp"""
    import ftplib
    import urllib

    schema, host, path, parameters, query, fragment = urllib.parse.urlparse(url)
    with ftplib.FTP(host) as ftp:
        try:
            ftp.login("anonymous", "")
            if "\n" in path:  # pragma: no cover
                raise ValueError("New line in path: %s" % (repr(path),))
            if path.endswith("/"):
                ftp.retrbinary("LIST " + path, file_object.write)
            else:
                ftp.retrbinary("RETR " + path, file_object.write)
        except ftplib.Error as e:
            raise ValueError("Error retrieving urls %s: %s" % (url, e))
