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
