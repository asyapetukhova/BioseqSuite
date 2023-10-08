def classify_amino_acid(seq: str) -> str:
    """
    Determine the percentage of acidic, basic and neutral amino acids in protein.

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code.

    Returns
    - output (str): percentage of neutral, acidic and basic amino acids in protein.
    """
    amino_acid_counts = {'acidic': 0, 'neutral': 0, 'basic': 0}
    aa_classification = {'acidic': ['D', 'E'], 'neutral': ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
                         'basic': ['H','R', 'K']}
    for i in seq:
        key = [key for key, value in aa_classification.items() if i in value]
        amino_acid_counts[key[0]] += 1
    acidic_percentage = round(amino_acid_counts['acidic'] / len(seq) * 100, 2)
    neutral_percentage = round(amino_acid_counts['neutral'] / len(seq) * 100, 2)
    basic_percentage = round(amino_acid_counts['basic'] / len(seq) * 100, 2)
    output = f'neutral: {neutral_percentage}, acidic: {acidic_percentage}, basic: {basic_percentage}'
    return output


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


def get_occurrences(seq1: str, seq2: str) -> str:
    """
   Counting the number of occurrences of string seq2 in string seq1.
   Getting indexes of occurrences of string seq2 in string seq1.

    Arguments:
    - seq1 (str): sequence in which search
    - seq2 (str): sequence to search in

    Return:
    - output (str): str, first element is the number of occurrences (int).
      All subsequent elements are indexes of occurrences of seq2 in seq1.
    """
    output = [seq1.count(seq2)]
    for i in range(len(seq1)):
        if seq1.startswith(seq2, i):
            output.append(i + 1)
    return (f'Number of occurrences: {output[0]}; '
            f'indexes: {", ".join(str(element) for element in output[1:])}')


def find_amino_acid_indices(seq: str, amino_acid: str) -> str:
    """
    Finds the amino acid indices specified in the input.

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code.
    - amino_acid (str): amino acid for which indices need to be found.

    Returns:
    - output (str): all found indices in the protein for which the entered amino acid corresponds to.
    """
    indices = []
    if amino_acid not in seq:
        raise ValueError('Amino acid not found')
    for index, aa in enumerate(seq):
        if aa == amino_acid:
            indices.append(index + 1)
    output = ', '.join(str(i) for i in indices)
    return output


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


def determine_total_protein_charge(seq) -> str:
    """
    Determine whether the protein has positive, negative or neutral charge in neutral pH

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code

    Returns:
    - output (str): positive, negative or neutral charge of protein in neutral pH
    """
    seq_list = list(seq.strip())
    aa_cnt = Counter(seq_list)
    number_of_pos = sum([aa_cnt[aa] for aa in POS_CHARGED])
    number_of_neg = sum([aa_cnt[aa] for aa in NEG_CHARGED])
    if number_of_pos == number_of_neg:
        return 'neutral'
    return 'positive' if number_of_pos > number_of_neg else 'negative'
