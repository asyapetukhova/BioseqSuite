def is_protein(seq: str, aa_set = set('FLIMVSPTAYHQNKDECWRG')) -> bool:
    """
    Check whether the transmitted sequence consists of amino acids

    Arguments:
    - seq (str): amino acid sequence. The input must be uppercased and use the single letter amino acid code
    - aa_set(set): all amino acid that we use. May be replaced with extra amino acids or it's modifications. But other functions aren't not intended to work with unusual amino acids

    Returns:
    - bool: True if the sequence is an amino acid sequence, False if it contains other symbols
    """
    return True if set(seq) <= AA_SET else False
