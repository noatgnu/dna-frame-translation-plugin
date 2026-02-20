"""
Core translation module for DNA/RNA to amino acid frame translation.

Provides functions for 3-frame and 6-frame translation with customizable
start codons, stop codons, and genetic code tables.
"""

from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass
from Bio import Seq
from Bio.Data import CodonTable


@dataclass
class TranslationResult:
    """Represents a single translated peptide."""

    sequence_id: str
    frame: int
    strand: str
    start_position: int
    end_position: int
    nucleotide_length: int
    peptide_sequence: str
    peptide_length: int


def get_codon_table(table_id: int = 1) -> CodonTable.CodonTable:
    """
    Get NCBI codon table by ID.

    Args:
        table_id: NCBI genetic code table number (1=Standard, 2=Vertebrate Mito, etc.)

    Returns:
        CodonTable object for translation.
    """
    return CodonTable.unambiguous_dna_by_id[table_id]


def find_start_positions(
    seq: str,
    start_codons: Optional[List[str]],
    frame_offset: int = 0
) -> List[int]:
    """
    Find all positions of start codons in a sequence for a given frame.

    Args:
        seq: DNA/RNA sequence string.
        start_codons: List of start codons to search for. None means start from position 0.
        frame_offset: Reading frame offset (0, 1, or 2).

    Returns:
        List of start positions (nucleotide indices).
    """
    if start_codons is None:
        return [frame_offset]

    positions = []
    seq_upper = seq.upper().replace("U", "T")
    seq_length = len(seq)

    for pos in range(frame_offset, seq_length - 2, 3):
        codon = seq_upper[pos:pos + 3]
        if codon in start_codons:
            positions.append(pos)

    return positions


def translate_frame(
    seq: str,
    start_position: int,
    stop_codons: Optional[List[str]],
    genetic_code: int = 1,
    use_stop_codons: bool = True
) -> Tuple[str, int]:
    """
    Translate a sequence from a given start position.

    Args:
        seq: DNA/RNA sequence string.
        start_position: Position to start translation.
        stop_codons: List of stop codons. None uses standard stops from codon table.
        genetic_code: NCBI genetic code table number.
        use_stop_codons: Whether to stop at stop codons.

    Returns:
        Tuple of (translated peptide sequence, end position in nucleotides).
    """
    seq_obj = Seq.Seq(seq.upper().replace("U", "T"))
    seq_length = len(seq_obj)

    fragment_end = start_position + 3 * ((seq_length - start_position) // 3)
    fragment = seq_obj[start_position:fragment_end]

    if use_stop_codons:
        if stop_codons is None:
            peptide = str(fragment.translate(table=genetic_code, to_stop=True))
        else:
            peptide_chars = []
            codon_table = get_codon_table(genetic_code)

            for i in range(0, len(fragment), 3):
                codon = str(fragment[i:i + 3])
                if len(codon) < 3:
                    break
                if codon in stop_codons:
                    break
                try:
                    aa = codon_table.forward_table.get(codon, "X")
                    peptide_chars.append(aa)
                except KeyError:
                    peptide_chars.append("X")

            peptide = "".join(peptide_chars)
    else:
        peptide = str(fragment.translate(table=genetic_code, to_stop=False))

    end_position = start_position + len(peptide) * 3
    return peptide, end_position


def three_frame_translation(
    seq: str,
    sequence_id: str,
    start_codons: Optional[List[str]] = None,
    stop_codons: Optional[List[str]] = None,
    use_start_codons: bool = False,
    use_stop_codons: bool = True,
    genetic_code: int = 1,
    min_length: int = 0,
    max_length: int = 0,
    reverse: bool = False
) -> List[TranslationResult]:
    """
    Perform 3-frame translation on a sequence.

    Args:
        seq: DNA/RNA sequence string.
        sequence_id: Identifier for the sequence.
        start_codons: List of start codons (e.g., ["ATG", "GTG"]).
        stop_codons: List of stop codons (e.g., ["TAA", "TAG", "TGA"]).
        use_start_codons: Whether to require start codons.
        use_stop_codons: Whether to stop at stop codons.
        genetic_code: NCBI genetic code table number.
        min_length: Minimum peptide length (0 = no filter).
        max_length: Maximum peptide length (0 = no limit).
        reverse: Whether to use reverse complement.

    Returns:
        List of TranslationResult objects.
    """
    results = []
    strand = "-" if reverse else "+"

    if reverse:
        seq_obj = Seq.Seq(seq.upper().replace("U", "T"))
        seq = str(seq_obj.reverse_complement())

    seq_length = len(seq)

    for frame in range(3):
        frame_num = frame + 1
        if reverse:
            frame_num = -(frame + 1)

        if use_start_codons and start_codons:
            start_positions = find_start_positions(seq, start_codons, frame)
            if not start_positions:
                continue
            start_pos = start_positions[0]
        else:
            start_pos = frame

        peptide, end_pos = translate_frame(
            seq,
            start_pos,
            stop_codons if use_stop_codons else None,
            genetic_code,
            use_stop_codons
        )

        if not peptide:
            continue

        peptide_len = len(peptide)

        if min_length > 0 and peptide_len < min_length:
            continue

        if max_length > 0 and peptide_len > max_length:
            continue

        results.append(TranslationResult(
            sequence_id=sequence_id,
            frame=frame_num,
            strand=strand,
            start_position=start_pos,
            end_position=end_pos,
            nucleotide_length=end_pos - start_pos,
            peptide_sequence=peptide,
            peptide_length=peptide_len
        ))

    return results


def six_frame_translation(
    seq: str,
    sequence_id: str,
    start_codons: Optional[List[str]] = None,
    stop_codons: Optional[List[str]] = None,
    use_start_codons: bool = False,
    use_stop_codons: bool = True,
    genetic_code: int = 1,
    min_length: int = 0,
    max_length: int = 0
) -> List[TranslationResult]:
    """
    Perform 6-frame translation on a sequence (3 forward + 3 reverse complement).

    Args:
        seq: DNA/RNA sequence string.
        sequence_id: Identifier for the sequence.
        start_codons: List of start codons.
        stop_codons: List of stop codons.
        use_start_codons: Whether to require start codons.
        use_stop_codons: Whether to stop at stop codons.
        genetic_code: NCBI genetic code table number.
        min_length: Minimum peptide length.
        max_length: Maximum peptide length.

    Returns:
        List of TranslationResult objects.
    """
    forward_results = three_frame_translation(
        seq=seq,
        sequence_id=sequence_id,
        start_codons=start_codons,
        stop_codons=stop_codons,
        use_start_codons=use_start_codons,
        use_stop_codons=use_stop_codons,
        genetic_code=genetic_code,
        min_length=min_length,
        max_length=max_length,
        reverse=False
    )

    reverse_results = three_frame_translation(
        seq=seq,
        sequence_id=sequence_id,
        start_codons=start_codons,
        stop_codons=stop_codons,
        use_start_codons=use_start_codons,
        use_stop_codons=use_stop_codons,
        genetic_code=genetic_code,
        min_length=min_length,
        max_length=max_length,
        reverse=True
    )

    return forward_results + reverse_results
