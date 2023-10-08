def protein_tool(*args: str) -> str:
    """
    Main function that is used to get sequence(s) and command. It performs a given action with the entered sequence

    Arguments:
    - args (str): amino acid sequence(s) and command.
    The input must use the single letter amino acid code
    The last element of the string must be the command

    Returns:
    - result (str): the result of a given sequence processing
    """
    *sequences, action = args
    sequences = [seq.upper() for seq in sequences]
    commands = {'calculate_amino_acid_percentages': calculate_amino_acid_percentages,
                'classify_amino_acid': classify_amino_acid,
                'find_amino_acid_indices': find_amino_acid_indices,
                'counting_point_mutations': counting_point_mutations,
                'counting_molecular_weight': counting_molecular_weight,
                'get_occurrences': get_occurrences,
                'count_variant_rna': count_variant_rna,
                'determine_total_protein_charge': determine_total_protein_charge,
                'calculate_pI': calculate_pI}
    command = commands[action]
    for seq in sequences:
        if not is_protein(seq):
            print('Sequence is not protein', file=sys.stderr)
            sys.exit(1)

    return str(command(sequences[0], sequences[1])) if len(sequences) == 2 else str(command(sequences[0]))
