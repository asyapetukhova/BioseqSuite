def counting_point_mutations(seq1: str, seq2: str) -> int:
    """
    Counts the number of mutations - amino acid substitutions in the sequence seq2 relative to seq1.
    Input sequences must have the same length.

    Arguments:
    - seq1 (str): sequence to compare with
    - seq2 (str): sequence to compare to

    Return:
    - output (int): number of amino acid substitutions
    """
    output = 0
    for number_amino_acid in range(len(seq1)):
        if seq1[number_amino_acid] != seq2[number_amino_acid]:
            output += 1
    return output


def counting_molecular_weight(seq: str) -> int:
    """
    Counts the molecular mass of a protein sequence seq

    Arguments:
    - seq (str): sequence to count the molecular weight

    Return:
    - output (int): molecular weight value
    """
    output = 0
    for amino_acid in seq:
        output += DICT_MOLECULAR_MASS[amino_acid]
    return output - 18 * (len(seq) - 1)


def count_variant_rna(seq: str) -> int:
    """
    Сounts the number of RNAs that can be a template for the synthesis of the entered sequence

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code

    Returns:
    - output (int): number of RNAs that can be a template for the synthesis of the entered sequence
    """
    output = 1
    for i in seq:
        output = output * NUMBER_CODONS[i]
    return output
