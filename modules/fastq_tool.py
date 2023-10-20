def fastq_tool(
    seqs: dict,
    gc_bounds: int | float | tuple = (0, 100),
    length_bounds: int | tuple = (0, 2**32),
    quality_threshold: float | int = 0,
) -> dict:
    """
    Filter fastq sequences by specified parameters

    Arguments:
    - seqs (dict): fastq sequences. Key - the name of the sequence (str). Value - tuple of type (sequence, quality).
    - gc_bounts (int|float|tuple): GC composition interval (in percent) for filtering (default is (0, 100)), i.e all reads passed.
    - length_bounds (int|tuple): filtering length interval (default is (0, 2**32)).
    - quality_threshold (float|int): threshold value of average read quality for filtering, default is 0 (phred33 scale).

    Returns:
    - seqs_filtered(dict): sequences that passed the filter.
    """
        seqs_filtered = {}
    for name, sequence in seqs.items():
        gc_content = (
            (sequence[0].count("G") + sequence[0].count("C")) / len(sequence[0]) * 100
        )  # check gc-content
        if isinstance(gc_bounds, tuple):
            if gc_content < gc_bounds[0] or gc_content > gc_bounds[1]:
                continue
        elif isinstance(gc_bounds, (int, float)):
            if gc_content > gc_bounds:
                continue
        else:
            raise ValueError(
                "gc_bounds must be of type tuple, int or float! Change your input values."
            )
