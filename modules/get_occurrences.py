def get_occurrences(seq1: str, seq2: str) -> str:
    """
   Counting the number of occurrences of string seq2 in string seq1.
   Getting indexes of occurrences of string seq2 in string seq1.

    Arguments:
    - seq1 (str): sequence in which search
    - seq2 (str): sequence to search in

    Return:
    - output (str): str, first element is the number of occurrences (int).
      All subsequent elements are indexes of occurrences of seq2 in seq1.
    """
    output = [seq1.count(seq2)]
    for i in range(len(seq1)):
        if seq1.startswith(seq2, i):
            output.append(i + 1)
    return (f'Number of occurrences: {output[0]}; '
            f'indexes: {", ".join(str(element) for element in output[1:])}')
