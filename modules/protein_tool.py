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
