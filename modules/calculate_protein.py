def calculate_amino_acid_percentages(seq: str) -> str:
    """
    Calculating the percentage of amino acids in protein.

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code.

    Returns:
    - output (str): percentage of amino acids in protein in descending order.
    """
    aa_count = {} #counting amino acids in sequence
    protein_length = len(seq)
    for amino_acid in seq:
        if amino_acid in aa_count:
            aa_count[amino_acid] += 1
        else:
            aa_count[amino_acid] = 1
    composition_rates = {}
    for aa in aa_count:
        composition_rates[aa] = aa_count[aa] / protein_length * 100
    output = ', '.join([ f'{key}: {round(value, 2)}' for key, value in sorted(composition_rates.items(),
                                               key=lambda item: -item[1])])
    return output


def calculate_pI(seq: str) -> float:
    """
    Calculation pI of the protein in neutral pH

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code

    Returns:
    - output (float): approximate value of the pI of the protein in neutral pH
    """
    seq_list = list(seq.strip())
    first_aa = seq_list[0]
    last_aa = seq_list[-1]
    aa_cnt = Counter(seq_list)

    summ_charge = [PK2[first_aa], PK1[last_aa]]

    for key, value in aa_cnt.items():
        try:
            summ_charge.extend([PK3[key] for _ in range(value)])
        except:
            pass

    return sum(summ_charge) / len(summ_charge)
