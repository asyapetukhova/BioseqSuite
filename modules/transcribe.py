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
