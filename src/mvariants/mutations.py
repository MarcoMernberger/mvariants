import pandas as pd
import mbf_genomes
import numpy as np
_mutation_ignored_transcripts = [
    "ENST00000361390"  # since it is not a complete coding sequence
]

def get_cds_sequence(genome, protein_id, protein_info=None):
    """Get the coding sequence (rna) of a protein"""
    if protein_info is None:
        protein_info = genome.df_proteins.loc[protein_id]
    elif protein_info.name != protein_id:
        raise ValueError("protein_id != protein_info['protein_id']")
    cdna = ""
    chr = protein_info["chr"]
    for start, stop in protein_info["cds"]:
        cdna += genome.get_genome_sequence(chr, start, stop)
    if protein_info["strand"] not in (1, -1):  # pragma: no cover
        raise ValueError(f'{protein_info["strand"]} was not 1/-1')
    if protein_info["strand"] == -1:
        cdna = reverse_complement(cdna)
    return cdna


def get_amino_acid_exchange(chr, pos, reference_base, actual_bases, genome):
    """Check whether a mutation is silent (False) or leads to an amino acid sequence in any transcript"""
    ret = {}
    altered = {}
    subs = []
    transcript_id = ""
    gene_name = ""
    gene_id = ""
    for mb in actual_bases.split(","):
        snp_df = _to_snp_df(chr, pos, reference_base, mb)
        altered_transcripts = snp_induced_protein_changes_full_info(snp_df, genome)
        if altered_transcripts:
            altered[mb] = altered_transcripts
            ret[mb] = {}
            for transcript_id in altered_transcripts:
                exchange = _get_amino_acid_exchange(
                    altered_transcripts[transcript_id]["original_aa"],
                    altered_transcripts[transcript_id]["mutated_aa"],
                    reference_base,
                    mb,
                )
                ret[mb][transcript_id] = exchange
            if len(altered_transcripts.keys()) > 0:
                for tid in altered_transcripts:
                    if transcript_id == "":
                        transcript_id = tid
                    subs.append(ret[mb][tid])
    gene_ids = {}
    to_df = {'gene_stable_id' : [], 'gene_name' : [], 'transcript_stable_id' : [], 'mb' : [], 'change' : []}
    for mb in ret:
        for transcript_id in altered_transcripts:
            gene_id = genome.df_transcripts.loc[transcript_id]['gene_stable_id']
            gene_name = genome.df_genes.loc[gene_id]['name']
            gene_ids[gene_id] = gene_name
            to_df['gene_stable_id'].append(gene_id)        
            to_df['gene_name'].append(gene_name)
            to_df['mb'].append(mb)
            to_df['transcript_stable_id'].append(transcript_id)
            to_df['change'].append(ret[mb][transcript_id])
    df = pd.DataFrame(to_df)
    return df, altered


def _get_amino_acid_exchange(
    original_AA_sequence, altered_AA_sequences, ref_nucleotides, altered_nucleotides
):
    """FF: Not sure what this does exactly.
    I think it turns the differences between orignal_AA and altered_aa it into a 'edit' representation
    """
    ret = ""
    for altered_AA_sequence in altered_AA_sequences.split(","):
        if len(ref_nucleotides) == 1 and len(altered_nucleotides) == 1:
            # SNP
            for num, aminoacid in enumerate(original_AA_sequence):
                if aminoacid != altered_AA_sequence[num]:
                    ret = aminoacid + str(num + 1) + altered_AA_sequence[num]
                    break
        elif len(ref_nucleotides) < len(altered_nucleotides):
            # insertion
            if (len(altered_nucleotides) - len(ref_nucleotides)) % 3 == 0:
                # inframe
                # inserted_len = (len(altered_nucleotides) - len(ref_nucleotides)) / 3
                for num, aminoacid in enumerate(original_AA_sequence):
                    if aminoacid != altered_AA_sequence[num]:
                        ret = (
                            aminoacid
                            + str(num + 1)
                            + "ins"
                            + original_AA_sequence[num + 1]
                        )
                        break
            else:
                # frameshift
                ret = "fs"
        elif len(ref_nucleotides) > len(altered_nucleotides):
            # deletion
            if (
                len(altered_AA_sequences) != 0
                and (len(ref_nucleotides) - len(altered_nucleotides)) % 3 == 0
            ):
                # inframe
                deleted_aa = (len(ref_nucleotides) - len(altered_nucleotides)) / 3
                for num, aminoacid in enumerate(altered_AA_sequence):
                    if aminoacid != original_AA_sequence[num]:
                        ret = (
                            original_AA_sequence[num]
                            + str(num + 1)
                            + "del"
                            + original_AA_sequence[num + 1 + deleted_aa]
                        )
                        break
            else:
                ret = "fs"
    return ret


def snp_induced_cdna_changes(snp_df, genome):
    """Given a df of {chr, pos, ref, alt} SNP information,
    this extracts all ensembl transcripts where the snps induce a protein change.
    Result is a dictionary: transcript_stable_id -> alterned dna sequence"""
    potentially_affected_transcripts = snp_hit_transcripts(snp_df, genome)
    result = {}
    for transcript_info in potentially_affected_transcripts.iter_rows():
        original_cdna, mutated_cdna = snp_changed_cds(
            snp_df, genome, transcript_info)
        result[transcript_info['transcript_stable_id']] = mutated_cdna
    return result


def snp_changed_cds(snp_df, genome, transcript_info_or_stable_id):
    """Retrieve original coding sequence+, altered coding sequence
    for a single transcript_info (=genome.get_transcript_info(stable_id)) or stable_id.
    Please note that while the start is defined, the end is simply the end of the cdan,
    not the end is not - you can either use get_translation_stop(), or look for the stop codon...
    """
    if isinstance(transcript_info_or_stable_id, str):
        transcript_info = genome.get_transcript_info(
            transcript_info_or_stable_id)
    else:
        transcript_info = transcript_info_or_stable_id
    reference_sequence = ""
    mutated_sequence = ""
    start_correction = 0
    ii = 0
    proteins = genome.df_proteins[genome.df_proteins['transcript_stable_id'] == transcript_info_or_stable_id[0]]
    if len(proteins) == 1:
        cds_positions = []
        for protein_start, protein_stop in proteins['cds'].values[0]:
            cds_positions.extend([protein_stop, protein_start])
        if transcript_info[1]['strand'] == -1:
            translation_start = np.max(cds_positions)
        else:
            translation_start = np.min(cds_positions)
    else:
        if transcript_info[1]['strand'] == -1:
            translation_start = transcript_info[1]['stop']
        else:
            translation_start = transcript_info[1]['start']
    for exon_start, exon_stop in transcript_info[1]['exons']: # no need to to anything for exons before the start, but do till the end, since we hight find something interesting
        exon_seq = genome.get_genome_sequence(
            transcript_info[1]['chr'], exon_start, exon_stop).lower()
        mutated_exon = exon_seq
        for i, row in snp_df.iterrows():
            pos = int(row['pos'])
            if row['chr'] == transcript_info[1]['chr']:
                if exon_start <= pos < exon_stop:
                    # snps always match the +1 strand
#                    if exon_seq[row['pos'] - exon_start: row['pos'] - exon_start + len(row['ref'])].upper() != row['ref'].upper():
                    try:
                        if row['ref'].upper().find(exon_seq[pos - exon_start: pos - exon_start + len(row['ref'])].upper()) != 0:
                            raise ValueError("Reference base did not match %s Was: %s" % (
                                row, exon_seq[pos - exon_start: pos - exon_start + len(row['ref'])]))
                    except TypeError:
                        print(pos, exon_start, len(row['ref']))
                        raise
                    mutated_exon = mutated_exon[
                        :pos - exon_start] + row['alt'].upper() + mutated_exon[pos - exon_start + len(row['ref']):]
                    if transcript_info[1]['strand'] == 1 and len(mutated_sequence) + (pos - exon_start) <=  translation_start: # we are before the ATG
                        start_correction +=  len(row['alt']) - len(row['ref'])
                    elif transcript_info[1]['strand'] == -1 and len(mutated_sequence) + exon_stop - pos >=  translation_start: # we are before the ATG
                        start_correction +=  len(row['alt']) - len(row['ref'])
        reference_sequence += exon_seq
        mutated_sequence += mutated_exon
        ii += 1
    start_offset = transcript_info[1]['start'] - translation_start
    if transcript_info[1]['strand'] == -1:
        start_offset = translation_start - transcript_info[1]['stop']
        reference_sequence = mbf_genomes.common.reverse_complement(
            reference_sequence
            )
        mutated_sequence = mbf_genomes.common.reverse_complement(
            mutated_sequence)

    # ensembl is 1 based..-1
    reference_sequence = reference_sequence[start_offset:]
    # ensembl is 1 based..
    mutated_sequence = mutated_sequence[start_offset + start_correction:]
    return reference_sequence, mutated_sequence

def _to_snp_df(chr, pos, reference_base, actual_bases):
    return pd.DataFrame(
        {"chr": [chr], "pos": [pos], "ref": [reference_base], "alt": [actual_bases]}
    )

def snp_induced_protein_changes(snp_df, genome):
    """Given a df of {chr, pos, ref, alt} SNP information,
    this extracts all ensembl transcripts where the snps induce a protein change.
    Result is a dictionary: transcript_stable_id -> altered amino acid sequence"""
    potentially_affected_transcripts = snp_hit_transcripts(snp_df, genome)
    result = {}
    for transcript_info in potentially_affected_transcripts.iterrows():
        original_cdna, mutated_cdna = snp_changed_cds(
            snp_df, genome, transcript_info)
        original_aa = genome.genetic_code.translate_dna_till_stop(original_cdna)
        mutated_aa = genome.genetic_code.translate_dna_till_stop(mutated_cdna)
        if original_aa != mutated_aa:
            result[transcript_info[1]['transcript_stable_id']] = mutated_aa
    return result


def snp_induced_protein_changes_full_info(snp_df, genome):
    """Given a df of {chr, pos, ref, alt} SNP information,
    this extracts all ensembl transcripts where the snps induce a protein change.
    Result is a dictionary: transcript_stable_id -> 
    {'original_cds': '', 'mutated_cds: '', 'original_aa': '', 'mutated_aa': ''}
    
    altered amino acid sequence"""
    potentially_affected_transcripts = snp_hit_transcripts(snp_df, genome)
    result = {}
    for transcript_info in potentially_affected_transcripts.iterrows():
        original_cdna, mutated_cdna = snp_changed_cds(
            snp_df, genome, transcript_info)
        original_aa = genome.genetic_code.translate_dna_till_stop(original_cdna)
        try:
            mutated_aa = genome.genetic_code.translate_dna_till_stop(mutated_cdna)
        except ValueError:
            mutated_aa = "NNNNNN"
        if original_aa != mutated_aa:
            result[transcript_info[0]] = {
                    'original_cds': original_cdna,
                    'mutated_cds': mutated_cdna,
                    'original_aa': original_aa,
                    'mutated_aa': mutated_aa,
                    }
    return result


def snp_hit_transcripts(snp_df, genome):
    """Given a df of {chr, pos, ref, alt} SNP information,
    this returns a set of transcript_stable_ids of transcripts where
    one of the snps is in one of the exons.
    result is a df with transcript_info
    python3 ready
    """
    hit_transcripts = set()
    hit_transcript_info = []
    for i, row in snp_df.iterrows():
        genes_overlapping = genome.get_genes_overlapping(
            row["chr"], row["pos"], row["pos"] + 1
        )  # returns None if there are no overlapping
        if genes_overlapping is not None:
            for stable_id, row_genes in genes_overlapping.iterrows():
                for transcript_stable_id in row_genes['transcript_stable_ids']:
                    # don't bother if it is already in the list
                    if not transcript_stable_id in hit_transcripts:
                        transcript_info = genome.df_transcripts.loc[transcript_stable_id]
                        
                        for exon_start, exon_stop in transcript_info["exons"]:
                            if (
                                row["chr"] == transcript_info["chr"]
                                and exon_start <= row["pos"] < exon_stop
                            ):
                                hit_transcripts.add(transcript_stable_id)
                                hit_transcript_info.append(transcript_info)
                                break
                        else:  # no exon hit...
                            continue
    if hit_transcript_info:
        return pd.DataFrame(hit_transcript_info)
    else:
        return pd.DataFrame(
            {"transcript_stable_id": [], "exons": [], "gene_stable_id": []}
        )

'''
def get_translation_stop(transcript_info):
    """given a transcript_info, how far into the concatenated (reversed) exon sequence is the stop codon?"""
    # there is a common ensembl bug where the translation is
    # annotated to be a bp shorter than it actually is
    for exon_start, exon_stop, exon_phase in transcript_info['exons']:
        stop_offset += exon_stop - exon_start
        if ii == transcript_info['translation_stop_exon']:
            break

    if ((stop_offset + transcript_info['translation_stop']) - (start_offset + transcript_info['translation_start'] - 1)) % 3 == 2:
        fix = 1
    # and some have only the first bp of a triplet, but it ain't a
    # stop, so we cut it of - it's irrelevant for the question
    # anyhow
    elif ((stop_offset + transcript_info['translation_stop']) - (start_offset + transcript_info['translation_start'] - 1)) % 3 == 1:
        fix = -1
    else:
        fix = 0
    return stop_offset + transcript_info['translation_stop'] + fix





def aa_short_changes(original_protein, mutated_protein):
    changes = []
    last_hit = -1
    run_start = -1
    for ii in xrange(min(len(original_protein), len(mutated_protein))):
        if original_protein[ii] != mutated_protein[ii]:
            changes.append(
                (original_protein[ii], (ii + 1), mutated_protein[ii]))
    ii = 0
    for a, b, c in zip(changes, changes[1:], changes[2:]):
        if a[1] == b[1] - 1 == c[1] - 2:  # we have a run...
            changes = changes[:ii + 1]
            changes[ii] = (changes[ii][0], changes[ii][1], 'fs')
            break
        ii += 1
    return "".join([(x[0] + str(x[1]) + x[2]) for x in changes])
'''