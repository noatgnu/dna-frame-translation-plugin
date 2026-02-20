#!/usr/bin/env python3
"""
DNA/RNA Frame Translation Plugin Runner.

Translates DNA/RNA sequences to amino acid sequences using 3-frame or 6-frame
translation with customizable start/stop codons and genetic code tables.
"""

import os
from pathlib import Path
from typing import List, Optional

import click
import pandas as pd
from Bio import SeqIO

from translator import (
    TranslationResult,
    three_frame_translation,
    six_frame_translation
)


def detect_format(file_path: str) -> str:
    """
    Detect file format based on extension.

    Args:
        file_path: Path to input file.

    Returns:
        Format string: 'fasta', 'csv', or 'tsv'.
    """
    ext = Path(file_path).suffix.lower()
    if ext in [".fasta", ".fa", ".fna", ".faa"]:
        return "fasta"
    elif ext == ".csv":
        return "csv"
    elif ext in [".tsv", ".txt"]:
        return "tsv"
    else:
        return "tsv"


def parse_codon_list(codon_string: str) -> List[str]:
    """
    Parse comma-separated codon string into list.

    Args:
        codon_string: Comma-separated codons (e.g., "ATG,GTG,CTG").

    Returns:
        List of uppercase codon strings.
    """
    if not codon_string or not codon_string.strip():
        return []
    return [c.strip().upper().replace("U", "T") for c in codon_string.split(",") if c.strip()]


def load_sequences_fasta(file_path: str) -> List[tuple]:
    """
    Load sequences from FASTA file.

    Args:
        file_path: Path to FASTA file.

    Returns:
        List of (id, sequence) tuples.
    """
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences


def load_sequences_tabular(
    file_path: str,
    sep: str,
    id_column: str,
    sequence_column: str
) -> List[tuple]:
    """
    Load sequences from tabular file (CSV/TSV).

    Args:
        file_path: Path to tabular file.
        sep: Column separator.
        id_column: Name of ID column.
        sequence_column: Name of sequence column.

    Returns:
        List of (id, sequence) tuples.
    """
    df = pd.read_csv(file_path, sep=sep)

    if id_column not in df.columns:
        available = ", ".join(df.columns.tolist()[:10])
        raise ValueError(f"ID column '{id_column}' not found. Available: {available}")

    if sequence_column not in df.columns:
        available = ", ".join(df.columns.tolist()[:10])
        raise ValueError(f"Sequence column '{sequence_column}' not found. Available: {available}")

    sequences = []
    for _, row in df.iterrows():
        seq_id = str(row[id_column])
        seq = row[sequence_column]
        if pd.notna(seq) and str(seq).strip():
            sequences.append((seq_id, str(seq).strip()))

    return sequences


def write_fasta(results: List[TranslationResult], output_path: str) -> int:
    """
    Write translation results to FASTA file.

    Args:
        results: List of TranslationResult objects.
        output_path: Output file path.

    Returns:
        Number of sequences written.
    """
    count = 0
    with open(output_path, "w") as f:
        for result in results:
            header = f">{result.sequence_id}|frame_{result.frame}|{result.strand}|pos_{result.start_position}-{result.end_position}"
            f.write(f"{header}\n{result.peptide_sequence}\n")
            count += 1
    return count


def write_summary(results: List[TranslationResult], output_path: str) -> None:
    """
    Write translation summary to TSV file.

    Args:
        results: List of TranslationResult objects.
        output_path: Output file path.
    """
    data = []
    for result in results:
        data.append({
            "sequence_id": result.sequence_id,
            "frame": result.frame,
            "strand": result.strand,
            "start_position": result.start_position,
            "end_position": result.end_position,
            "nucleotide_length": result.nucleotide_length,
            "peptide_length": result.peptide_length,
            "peptide_sequence": result.peptide_sequence
        })

    df = pd.DataFrame(data)
    df.to_csv(output_path, sep="\t", index=False)


def write_report(
    results: List[TranslationResult],
    total_sequences: int,
    parameters: dict,
    output_path: str
) -> None:
    """
    Write translation report.

    Args:
        results: List of TranslationResult objects.
        total_sequences: Total input sequences processed.
        parameters: Dictionary of run parameters.
        output_path: Output file path.
    """
    with open(output_path, "w") as f:
        f.write("DNA/RNA Frame Translation Report\n")
        f.write("=" * 50 + "\n\n")

        f.write("Parameters:\n")
        f.write("-" * 30 + "\n")
        for key, value in parameters.items():
            f.write(f"  {key}: {value}\n")
        f.write("\n")

        f.write("Results Summary:\n")
        f.write("-" * 30 + "\n")
        f.write(f"  Input sequences: {total_sequences}\n")
        f.write(f"  Output peptides: {len(results)}\n")

        if results:
            lengths = [r.peptide_length for r in results]
            f.write(f"  Min peptide length: {min(lengths)} aa\n")
            f.write(f"  Max peptide length: {max(lengths)} aa\n")
            f.write(f"  Mean peptide length: {sum(lengths) / len(lengths):.1f} aa\n")

            frame_counts = {}
            for r in results:
                frame_counts[r.frame] = frame_counts.get(r.frame, 0) + 1

            f.write("\n")
            f.write("Peptides per frame:\n")
            for frame in sorted(frame_counts.keys()):
                f.write(f"  Frame {frame:+d}: {frame_counts[frame]}\n")

            strand_counts = {"+": 0, "-": 0}
            for r in results:
                strand_counts[r.strand] += 1

            f.write("\n")
            f.write("Peptides per strand:\n")
            f.write(f"  Forward (+): {strand_counts['+']}\n")
            f.write(f"  Reverse (-): {strand_counts['-']}\n")

        f.write("\n")
        f.write("=" * 50 + "\n")
        f.write("Translation completed successfully.\n")


@click.command()
@click.option("--input_file", required=True, help="Input sequence file (FASTA/CSV/TSV)")
@click.option("--input_format", default="auto", help="Input format (auto/fasta/csv/tsv)")
@click.option("--sequence_column", default="sequence", help="Sequence column name for tabular input")
@click.option("--id_column", default="gene", help="ID column name for tabular input")
@click.option("--translation_mode", default="six_frame", help="Translation mode (three_frame/six_frame)")
@click.option("--use_start_codons", is_flag=True, help="Require start codons")
@click.option("--start_codons", default="ATG", help="Comma-separated start codons")
@click.option("--use_stop_codons", is_flag=True, help="Stop at stop codons")
@click.option("--stop_codons", default="TAA,TAG,TGA", help="Comma-separated stop codons")
@click.option("--min_length", default=0, type=int, help="Minimum peptide length")
@click.option("--max_length", default=0, type=int, help="Maximum peptide length (0=no limit)")
@click.option("--genetic_code", default="1", help="NCBI genetic code table ID")
@click.option("--output_dir", required=True, help="Output directory")
def main(
    input_file: str,
    input_format: str,
    sequence_column: str,
    id_column: str,
    translation_mode: str,
    use_start_codons: bool,
    start_codons: str,
    use_stop_codons: bool,
    stop_codons: str,
    min_length: int,
    max_length: int,
    genetic_code: str,
    output_dir: str
) -> None:
    """
    Translate DNA/RNA sequences to amino acids using frame translation.
    """
    print("[1/5] Parsing arguments...")

    os.makedirs(output_dir, exist_ok=True)

    if input_format == "auto":
        input_format = detect_format(input_file)

    start_codon_list: Optional[List[str]] = None
    if use_start_codons:
        start_codon_list = parse_codon_list(start_codons)
        if not start_codon_list:
            start_codon_list = ["ATG"]

    stop_codon_list: Optional[List[str]] = None
    if use_stop_codons:
        stop_codon_list = parse_codon_list(stop_codons)
        if not stop_codon_list:
            stop_codon_list = ["TAA", "TAG", "TGA"]

    genetic_code_int = int(genetic_code)

    parameters = {
        "input_file": os.path.basename(input_file),
        "input_format": input_format,
        "translation_mode": translation_mode,
        "use_start_codons": use_start_codons,
        "start_codons": start_codon_list if use_start_codons else "N/A (from position 0)",
        "use_stop_codons": use_stop_codons,
        "stop_codons": stop_codon_list if use_stop_codons else "N/A (translate entire frame)",
        "min_length": min_length if min_length > 0 else "No filter",
        "max_length": max_length if max_length > 0 else "No limit",
        "genetic_code": genetic_code_int
    }

    print("[2/5] Loading sequences...")

    if input_format == "fasta":
        sequences = load_sequences_fasta(input_file)
    elif input_format == "csv":
        sequences = load_sequences_tabular(input_file, ",", id_column, sequence_column)
    else:
        sequences = load_sequences_tabular(input_file, "\t", id_column, sequence_column)

    total_sequences = len(sequences)
    print(f"       Loaded {total_sequences} sequences")

    print("[3/5] Translating sequences...")

    all_results: List[TranslationResult] = []
    processed = 0

    translate_func = six_frame_translation if translation_mode == "six_frame" else three_frame_translation

    for seq_id, seq in sequences:
        if translation_mode == "six_frame":
            results = six_frame_translation(
                seq=seq,
                sequence_id=seq_id,
                start_codons=start_codon_list,
                stop_codons=stop_codon_list,
                use_start_codons=use_start_codons,
                use_stop_codons=use_stop_codons,
                genetic_code=genetic_code_int,
                min_length=min_length,
                max_length=max_length
            )
        else:
            results = three_frame_translation(
                seq=seq,
                sequence_id=seq_id,
                start_codons=start_codon_list,
                stop_codons=stop_codon_list,
                use_start_codons=use_start_codons,
                use_stop_codons=use_stop_codons,
                genetic_code=genetic_code_int,
                min_length=min_length,
                max_length=max_length,
                reverse=False
            )

        all_results.extend(results)
        processed += 1

        if processed % 100 == 0:
            print(f"       Processed {processed}/{total_sequences} sequences...", end="\r")

    print(f"       Processed {processed}/{total_sequences} sequences")
    print(f"       Generated {len(all_results)} peptides")

    print("[4/5] Writing output files...")

    fasta_path = os.path.join(output_dir, "translated_sequences.fasta")
    write_fasta(all_results, fasta_path)
    print(f"       Written: {fasta_path}")

    summary_path = os.path.join(output_dir, "translation_summary.tsv")
    write_summary(all_results, summary_path)
    print(f"       Written: {summary_path}")

    report_path = os.path.join(output_dir, "translation_report.txt")
    write_report(all_results, total_sequences, parameters, report_path)
    print(f"       Written: {report_path}")

    print("[5/5] Translation completed successfully!")
    print(f"       Input: {total_sequences} sequences")
    print(f"       Output: {len(all_results)} peptides")


if __name__ == "__main__":
    main()
