# Constants
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
