"""
Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.

Return: The adjacency list corresponding to O3. You may return edges in any order.
"""

def overlap_graphs(sequences, ids, l = 3):
    """
    Generate an overlap graph from a list of DNA sequences.
    
    Parameters:
        sequences (list of str): list of DNA sequences.
        ids (list of str): list of identifiers corresponding to the sequences.
        l (int): length of the overlap to consider (default is 3).
        
    Returns:
        list of tuples: each tuple contains two identifiers (id1, id2) where the
                        sequence with id1 overlaps with the sequence with id2.
    """
    result = []
    
    for idx1, seq1 in enumerate(sequences):
        for idx2, seq2 in enumerate(sequences):
            # Skip if the same sequence is compared
            if idx1 == idx2:
                continue
            # Check if seq1 ends with the first l characters of seq2
            if str(seq1).endswith(seq2[:l]):
                result.append((ids[idx1], ids[idx2]))
    
    return result

if __name__ == "__main__":
    dna_sequences = """"""
    
    with open('./data/rosalind_grph.txt', 'r') as file:
        dna_sequences = file.read()
        
    # Create a lists with all the sequences and their ids
    ids, sequences, seq = [], [], ""

    for line in dna_sequences.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if seq:
                sequences.append(seq)
                seq = ""
            ids.append(line[1:])
        else:
            seq += line

    if seq:
        sequences.append(seq)
        
    # Input validation
    for idx, sequence_list in enumerate(sequences):
        seq = ''.join(sequence_list)
        # Check length of each sequence
        if len(seq) > 10000:
            raise ValueError(f"Sequence {ids[idx]} is too long, must be 10000 characters or less.")
        # Check nucleotides
        if not all(nt in 'ACGT' for nt in seq.upper()):
            raise ValueError(f"Sequence {ids[idx]} contains invalid characters, only A, C, G, T are allowed.")
    
    edges = overlap_graphs(sequences, ids)
    # Output the edges to a file
    with open("./output/output_grph.txt", "w") as f:
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]}\n")