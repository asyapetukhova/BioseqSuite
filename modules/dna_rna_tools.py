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


def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complementary sequence

    Arguments:
    - seq (str): DNA or RNA sequence

    Returns:
    - output (str): Reverse complement sequence
    """
    return complement(reverse(seq))


def run_dna_rna_tools(*args: str) -> str:
    """
    Performs the specified operation on the given sequences

    Arguments:
    - args (str): RNA/DNA sequence(s) and command.
    The last element of the string must be the command

    Returns:
    - result (str): the result of a given sequence processing
    """
    seqs, operation = args[:-1], args[-1]
    results = []  # Checks the operation
    commands = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    command = commands[operation]
    for seq in seqs:
        if not set(seq).issubset(DNA_BASES) and not set(seq).issubset(RNA_BASES):
            raise ValueError("Not a DNA or RNA! Check your entered sequence.")
        results.append(command(seq))
    # Return the result
    if len(results) == 1:
        return results[0]
    else:
        return results
