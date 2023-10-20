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
