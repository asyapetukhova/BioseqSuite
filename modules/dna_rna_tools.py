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


def transcribe(seq: str) -> str:
    """
    Transcribes DNA sequence into RNA

    Arguments:
    - seq (str): DNA sequence

    Returns:
    - output (str): RNA sequence
    """
    if "U" in seq or "u" in seq:
        raise ValueError("Can't transcribe RNA!")
    else:
        return seq.replace("T", "U").replace("t", "u")


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
