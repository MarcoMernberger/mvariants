import pypipegraph2 as ppg
import rpy2.robjects as ro
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Union, Any, Callable
from pypipegraph2 import Job
from .caller import _Caller
from mbf.r import convert_dataframe_to_r
from mbf.genomes.ensembl import _EnsemblGenome


def summary__varscan_matched_dump_exons(
    output_file: Path,
    genome: _EnsemblGenome,
    dependencies: List[Job] = [],
) -> Job:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    def __dump(outfile):
        genome.df_exons.to_csv(output_file, sep="\t")

    return ppg.FileGeneratingJob(output_file, __dump).depends_on(dependencies)


def summary__varscan_matched(
    output_file: Path,
    sample_to_result_files: Dict[str, Path],
    caller: _Caller,
    genome: _EnsemblGenome,
    filter_func: Optional[Callable] = None,
    dependencies: List[Job] = [],
) -> Job:
    """
    summary__varscan_matched creates a job that summarizes varscan outputs.

    _extended_summary_

    Parameters
    ----------
    output_file : Path
        The output summary file to be generated.
    sample_to_result_files : Dict[str, Path]
        a dictrionary mapping sample names to the varscan result files
        (.vcf.snp and .vcf.indel).
    caller : _Caller
        The Varscan instance used.
    dependencies : List[Job], optional
        List of Jobs, should include the pileups, by default [].

    Returns
    -------
    Job
        _description_

    Raises
    ------
    ValueError
        _description_
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    def __dump(outfile):
        dfs_to_concat = []
        rename_columns = {
            "ALT": "alt",
            "REF": "ref",
            "#CHROM": "chrom",
            "POS": "pos",
        }
        for sample_name in sample_to_result_files:
            for filename in sample_to_result_files[sample_name]:
                # df = caller.read_vcf(filename)
                df = pd.read_csv(filename, sep="\t", skiprows=17)

                df_normal = df["NORMAL"].str.split(":", expand=True)
                df_normal.columns = [f"{s} NORMAL" for s in df["FORMAT"].values[0].split(":")]
                df_tumor = df["TUMOR"].str.split(":", expand=True)
                df_tumor.columns = [f"{s} TUMOR" for s in df["FORMAT"].values[0].split(":")]
                df = df.rename(columns=rename_columns)

                df_expand = pd.concat([df, df_normal, df_tumor], axis=1)
                df_expand["AD TUMOR"] = df_expand["AD TUMOR"].astype(int)
                df_expand["RD TUMOR"] = df_expand["RD TUMOR"].astype(int)
                df_expand["DP TUMOR"] = df_expand["DP TUMOR"].astype(int)
                df_expand["AD NORMAL"] = df_expand["AD NORMAL"].astype(int)
                df_expand["RD NORMAL"] = df_expand["RD NORMAL"].astype(int)
                df_expand["DP NORMAL"] = df_expand["DP NORMAL"].astype(int)

                df_expand["AF TUMOR"] = (
                    df_expand["AD TUMOR"].astype(float).div(df_expand["DP TUMOR"].astype(float))
                )
                df_expand["RF TUMOR"] = (
                    df_expand["RD TUMOR"].astype(float).div(df_expand["DP TUMOR"].astype(float))
                )
                df_expand["AF NORMAL"] = (
                    df_expand["AD NORMAL"].astype(float).div(df_expand["DP NORMAL"].astype(float))
                )
                df_expand["RF NORMAL"] = (
                    df_expand["RD NORMAL"].astype(float).div(df_expand["DP NORMAL"].astype(float))
                )
                df_expand["FREQ TUMOR"] = df_expand["FREQ TUMOR"].str[:-1].astype(float)
                df_expand["FREQ NORMAL"] = df_expand["FREQ NORMAL"].str[:-1].astype(float)

                def get_filter(df_exons):
                    maps = dict([(chrm, set()) for chrm in df_exons["chr"].unique()])
                    for _, row in df_exons.iterrows():
                        maps[row["chr"]].update(list(np.arange(row["start"], row["stop"] + 1)))

                    def __filter(row):
                        try:
                            return row["pos"] in maps[row["chrom"]]
                        except KeyError:
                            return False

                    return __filter

                df_expand["Exon"] = df_expand.apply(get_filter(genome.df_exons), axis=1)

                # df = df[(df["somatic_status"] == "Somatic")]  # check only somatic mutations
                # df["tumor_var_freq"] = [float(x[:-1]) / 100 for x in df["tumor_var_freq"].values]
                # df["reference_frequency"] = np.ones(len(df)) - df["tumor_var_freq"].values
                # df["normal_var_freq"] = [float(x[:-1]) / 100 for x in df["normal_var_freq"].values]
                df_expand["sample"] = [sample_name] * len(df)
                if ".indel." in str(filename):
                    df_expand["caller"] = [f"{caller.name}_indel"] * len(df)
                elif ".snp." in str(filename):
                    df_expand["caller"] = [f"{caller.name}_snp"] * len(df)
                else:
                    raise ValueError(f"Don't know how to handle the file {filename}.")

                dfs_to_concat.append(df_expand)
        df_ret = pd.concat(dfs_to_concat)
        df_ret = df_ret.rename(columns=rename_columns)
        df_ret.to_csv(str(outfile), sep="\t", index=False)

    return ppg.FileGeneratingJob(output_file, __dump).depends_on(dependencies)


def signature_analysis(
    output_files: Dict[str, Path],
    summary_file_input: Path,
    dependencies: List[Job] = [],
    used_callers: Optional[List[str]] = None,
    filterfunc: Optional[Callable] = None,
    dbfile: Optional[Path] = None,
) -> Job:
    """
    signature_analysis performs a signature analysis using deconstructSigs.

    _extended_summary_

    Parameters
    ----------
    output_files : Dict[str, List[str]]
        the
    summary_file_input : Path
        _description_
    dependencies : List[Job], optional
        _description_, by default []
    used_callers : Optional[List[str]], optional
        _description_, by default None
    filterfunc : Optional[Callable], optional
        a filterfuction to be applied on the summary dataframe to filter out
        unwanted snps (e.g. intron mutations), by default None

    Returns
    -------
    MultiFileGeneratingJob
        The job that creates the output files.
    """
    # for sample_name in output_files:
    #    output_files[sample_name].parent.mkdir(parents=True, exist_ok=True)
    first_output = list(output_files.values())[0]
    parent_folder = first_output.parent
    parent_folder.mkdir(parents=True, exist_ok=True)
    mutout_file = str(parent_folder / f"{first_output.stem}.signature_mutations.tsv")
    logfile = str(parent_folder / f"{first_output.stem}.signature_mutations.log")
    weight_file = parent_folder / f"{first_output.stem}.weights.tsv"

    def __dump(*outfiles):
        ro.r("library(deconstructSigs)")
        ro.r("library(RColorBrewer)")
        ro.r("mm39 = BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39")
        df_mutations = pd.read_csv(summary_file_input, sep="\t")
        if filterfunc is not None:
            df_mutations = filterfunc(df_mutations)
        df_mutations.to_csv(mutout_file, sep="\t", index=False)
        canonical = [str(i) for i in range(1, 20)] + ["X", "Y", "MT"]
        df_mutations = df_mutations[[c in canonical for c in df_mutations["chrom"].values]]
        suffixes = [".signature.png", ".pie.png", ".tumor.png", ".weights.csv"]
        samples = output_files.keys()
        concat_weights = []
        with Path(logfile).open("w") as outl:
            for sample_name in samples:
                df = df_mutations[df_mutations["sample"] == sample_name].copy()
                if len(df) > 50:
                    df["chrom"] = [f"chr{x}" if x != "MT" else "chrM" for x in df["chrom"].values]
                    # we might want to exclude indels here ...
                    if used_callers is not None:
                        df = df[df["caller"].isin(used_callers)]
                    df = df[["sample", "chrom", "pos", "ref", "alt"]]
                    # this should reproduce the R error
                    # df["pos"] = df["pos"]-1
                    ro.r(
                        """
                        f = function(df_mut, sample_id, filename_signature, filename_pie, filename_tumor, filename_weights, dbfile)
                        {
                            tryCatch(
                                {
                                    print(dbfile)
                                    if (dbfile == "None"){
                                        signature.db = deconstructSigs::signatures.cosmic
                                    }
                                    else{
                                        signature.db = read.csv(dbfile, sep="\t", check.names = FALSE, header = TRUE, row.names="Type")
                                    }
                                    si = deconstructSigs::mut.to.sigs.input(df_mut, sample.id = "sample", chr = "chrom", pos = "pos", ref = "ref", alt = "alt", bsg = mm39)
                                    sigs.output = deconstructSigs::whichSignatures(tumor.ref = si, signatures.ref = signature.db, sample.id = sample_id, contexts.needed = TRUE, tri.counts.method = deconstructSigs::tri.counts.exome)
                                    write.table(sigs.output$weights, file=filename_weights)
                                    png(filename_signature)
                                    deconstructSigs::plotSignatures(sigs.output)
                                    dev.off()
                                    png(filename_pie)
                                    deconstructSigs::makePie(sigs.output, add.color=c(brewer.pal(n = sum(sigs.output$weights > 0), name="Spectral")))
                                    dev.off()
                                    png(filename_tumor)
                                    deconstructSigs::plotTumor(sigs.output[['tumor']])
                                    dev.off()
                                },
                                warning = function(w)
                                {
                                    print(w)
                                    stop(w)
                                },
                                error = function(e)
                                {
                                print(e)
                                stop(e)
                                }
                            )
                        }
                        """
                    )(
                        convert_dataframe_to_r(df),
                        sample_name,
                        str(output_files[sample_name].with_suffix(suffixes[0])),
                        str(output_files[sample_name].with_suffix(suffixes[1])),
                        str(output_files[sample_name].with_suffix(suffixes[2])),
                        str(output_files[sample_name].with_suffix(suffixes[3])),
                        str(dbfile),
                    )
                    outl.write(f"For {sample_name} {len(df)} SNPs were used for decomposition.\n")
                    concat_weights.append(
                        pd.read_csv(output_files[sample_name].with_suffix(suffixes[3]), sep=" ")
                    )
                else:
                    outl.write(
                        f"For {sample_name} {len(df)} SNPs remain after filtering. Samples with less than 50 SNPs are discarded.\n"
                    )
            df_weights = pd.concat(concat_weights)
            df_weights.to_csv(weight_file, sep="\t")

    return ppg.MultiFileGeneratingJob([logfile, mutout_file], __dump).depends_on(dependencies)
