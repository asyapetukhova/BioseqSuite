import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
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


def filter_fastq(
    input_path: str, 
    output_path: str, 
    gc_bounds: int | float | tuple = (0, 100), 
    length_bounds: int | tuple = (0, 2**32), 
    quality_threshold: float | int = 0
) -> dict:
    
    def filter_record(record: SeqRecord) -> bool:
        seq_len: int = len(record.seq)
        mean_quality: float = sum(record.letter_annotations["phred_quality"]) / seq_len
        gc_content: float = gc_fraction(record.seq)
        
        return (length_bounds[0] <= seq_len <= length_bounds[1] and
                mean_quality >= quality_threshold and
                gc_bounds[0] <= gc_content <= gc_bounds[1])
    
    with open(input_path, "r") as input_handle, open(output_path, "w") as output_handle:
        records = SeqIO.parse(input_handle, "fastq")
        filtered_records = (record for record in records if filter_record(record))
        SeqIO.write(filtered_records, output_handle, "fastq")
