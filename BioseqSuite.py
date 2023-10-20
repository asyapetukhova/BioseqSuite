import sys
from collections import Counter

from modules.protein_functions import (
    calculate_amino_acid_percentages,
    calculate_pI,
    classify_amino_acid,
    counting_point_mutations,
    counting_molecular_weight,
    count_variant_rna,
    find_amino_acid_indices,
    get_occurrences,
    is_protein,
    determine_total_protein_charge
)
from modules.dna_rna_functions import (
    reverse,
    complement,
    reverse_complement,
    transcribe
)


#Constants
DNA_BASES = "ACGTacgt"
RNA_BASES = "ACGUacgu"


def protein_tool(*args: str) -> str:
    """
    Main function that is used to get sequence(s) and command. It performs a given action with the entered sequence

    Arguments:
    - args (str): amino acid sequence(s) and command.
    The input must use the single letter amino acid code
    The last element of the string must be the command

    Returns:
    - result (str): the result of a given sequence processing
    """
    *sequences, action = args
    sequences = [seq.upper() for seq in sequences]
    commands = {'calculate_amino_acid_percentages': calculate_amino_acid_percentages,
                'classify_amino_acid': classify_amino_acid,
                'find_amino_acid_indices': find_amino_acid_indices,
                'counting_point_mutations': counting_point_mutations,
                'counting_molecular_weight': counting_molecular_weight,
                'get_occurrences': get_occurrences,
                'count_variant_rna': count_variant_rna,
                'determine_total_protein_charge': determine_total_protein_charge,
                'calculate_pI': calculate_pI}
    command = commands[action]
    for seq in sequences:
        if not is_protein(seq):
            print('Sequence is not protein', file=sys.stderr)
            sys.exit(1)

    return str(command(sequences[0], sequences[1])) if len(sequences) == 2 else str(command(sequences[0]))


def run_dna_rna_tools(*args: str) -> str:
    """
    Performs the specified operation on the given sequences

    Arguments:
    - args (str): RNA/DNA sequence(s) and command.
    The last element of the string must be the command

    Returns:
    - result (str): the result of a given sequence processing
    """
    seqs, operation = args[:-1], args[-1]
    results = []  # Checks the operation
    commands = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    command = commands[operation]
    for seq in seqs:
        if not set(seq).issubset(DNA_BASES) and not set(seq).issubset(RNA_BASES):
            raise ValueError("Not a DNA or RNA! Check your entered sequence.")
        results.append(command(seq))
    # Return the result
    if len(results) == 1:
        return results[0]
    else:
        return results


def fastq_tool(
    seqs: dict,
    gc_bounds: int | float | tuple = (0, 100),
    length_bounds: int | tuple = (0, 2**32),
    quality_threshold: float | int = 0,
) -> dict:
    """
    Filter fastq sequences by specified parameters

    Arguments:
    - seqs (dict): fastq sequences. Key - the name of the sequence (str). Value - tuple of type (sequence, quality).
    - gc_bounts (int|float|tuple): GC composition interval (in percent) for filtering (default is (0, 100)), i.e all reads passed.
    - length_bounds (int|tuple): filtering length interval (default is (0, 2**32)).
    - quality_threshold (float|int): threshold value of average read quality for filtering, default is 0 (phred33 scale).

    Returns:
    - seqs_filtered(dict): sequences that passed the filter.
    """
    seqs_filtered = {}
    for name, sequence in seqs.items():
        gc_content = (
            (sequence[0].count("G") + sequence[0].count("C")) / len(sequence[0]) * 100
        )
        if isinstance(gc_bounds, tuple):
            if gc_content < gc_bounds[0] or gc_content > gc_bounds[1]:
                continue
        elif isinstance(gc_bounds, (int, float)):
            if gc_content > gc_bounds:
                continue
        else:
            raise ValueError(
                "gc_bounds must be of type tuple, int or float! Change your input values."
            )
        if isinstance(length_bounds, tuple):
            if (
                len(sequence[0]) < length_bounds[0]
                or len(sequence[0]) > length_bounds[1]
            ):
                continue
        elif isinstance(length_bounds, int):
            if len(sequence[0]) > length_bounds:
                continue
        else:
            raise ValueError(
                "length_bounds must be of type tuple or int! Change your input values."
            )
        mean_quality = sum([ord(q) - 33 for q in sequence[1]]) / len(sequence[1])
        if isinstance(quality_threshold, (int, float)):
            if mean_quality < quality_threshold:
                continue
        else:
            raise ValueError("quality_threshold must be of type int or float.")
        seqs_filtered[name] = sequence
    if seqs_filtered == {}:
        print("No sequences found! Change your input parameters.")
    else:
        return seqs_filtered
