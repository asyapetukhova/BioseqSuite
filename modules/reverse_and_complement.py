def reverse(seq: str) -> str:
    """
    Reverses the sequence

    Arguments:
    - seq (str): DNA or RNA sequence

    Returns:
    - output (str): reversed sequence
    """
    return seq[::-1]


def complement(seq: str) -> str:
    """
    Returns the complementary sequence

    Arguments:
    - seq (str): DNA or RNA sequence

    Returns:
    - output (str): complementary sequence
    """
    if "U" in seq or "u" in seq:
        return "".join([RNA_DICT[base] for base in seq])
    else:
        return "".join([DNA_DICT[base] for base in seq])


def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complementary sequence

    Arguments:
    - seq (str): DNA or RNA sequence

    Returns:
    - output (str): Reverse complement sequence
    """
    return complement(reverse(seq))
