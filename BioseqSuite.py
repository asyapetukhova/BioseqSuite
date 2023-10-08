from collections import Counter
from modules.calculate_protein import calculate_amino_acid_percentages, calculate_pI
import sys

# Constants
NUMBER_CODONS = {'F': 2, 'L': 6, 'I': 3, 'M': 1, 'V': 4, 'S': 6, 'P': 4, 'T': 4,
                 'A': 4, 'Y': 2, 'H': 2, 'Q': 2, 'N': 2, 'K': 2, 'D': 2, 'E': 2,
                 'C': 2, 'W': 1, 'R': 6, 'G': 4
                 }
NEG_CHARGED = ['D', 'E']
POS_CHARGED = ['H', 'K', 'R']
PK1 = {'F': 2.2, 'L': 2.36, 'I': 2.36, 'M': 2.28,
       'V': 2.32, 'S': 2.21, 'P': 1.99, 'T': 2.71,
       'A': 2.34, 'Y': 2.2, 'H': 1.82, 'Q': 2.17,
       'N': 2.02, 'K': 2.18, 'D': 1.88, 'E': 2.19,
       'C': 1.71, 'W': 2.38, 'R': 2.17, 'G': 2.34
       }
PK2 = {'F': 9.09, 'L': 9.6, 'I': 9.68, 'M': 9.21,
       'V': 9.62, 'S': 9.15, 'P': 10.96, 'T': 9.62,
       'A': 9.69, 'Y': 9.11, 'H': 9.17, 'Q': 9.13,
       'N': 9.8, 'K': 8.95, 'D': 9.6, 'E': 9.67,
       'C': 8.33, 'W': 9.39, 'R': 9.04, 'G': 9.6
       }
PK3 = {'Y': 10.07, 'H': 6.0, 'K': 10.53,
       'C': 10.78, 'D': 3.65, 'E': 4.25, 'R': 12.48
       }
AA_SET = set('FLIMVSPTAYHQNKDECWRG')
DICT_MOLECULAR_MASS = {
    'G': 75, 'A': 89, 'V': 117, 'L': 131, 'I': 131, 'P': 115,
    'F': 165, 'Y': 181, 'W': 204, 'S': 105, 'T': 119, 'C': 121,
    'M': 149, 'N': 132, 'Q': 146, 'D': 133, 'E': 147, 'K': 146,
    'R': 174, 'H': 155
}
DNA_DICT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
}
RNA_DICT = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
    "a": "u",
    "u": "a",
    "c": "g",
    "g": "c",
}
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
