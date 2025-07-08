"""
Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several possible consensus 
    strings exist, then you may return any one of them.)
"""

import numpy as np

NUCLEOTIDE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
INDEX_MAP = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

def get_consensus(profile, cols):
    """
    Generate a consensus string from the profile matrix.
    
    Parameters:
        profile (numpy.ndarray): profile matrix where each row corresponds to a nucleotide (A, C, G, T) 
            and each column corresponds to a position in the sequences.
        cols (int): number of columns in the profile matrix.
        
    Return:
        consensus (str): consensus string derived from the profile matrix.
    """
    consensus = []
    for c in range(cols):
        # Find the index of the maximum value in the column c
        # This index corresponds to the nucleotide with the highest count at that position
        idx_max = profile[:,c].argmax(axis=0)
        # Append the corresponding nucleotide to the consensus string
        # Using INDEX_MAP to convert index to nucleotide
        consensus.append(INDEX_MAP[idx_max])
        
    return ''.join(consensus)

def get_profile(sequences):
    """
    Create a profile matrix and consensus string from a list of DNA sequences.
    
    Parameters:
        sequences (list of str): list of DNA sequences.
        
    Return:
        profile (numpy.ndarray): profile matrix where each row corresponds to a nucleotide (A, C, G, T) 
            and each column corresponds to a position in the sequences.
        consensus (str): consensus string derived from the profile matrix.
    """
    
    cols = len(sequences[0])
    rows = 4
    profile = np.zeros((rows, cols), dtype=int)
    
    for c in range(cols):
        for seq in sequences:
            # For each nucleotide in the sequence, update the profile matrix
            nucleotide = seq[c]
            # Get the index corresponding to the nucleotide
            r = NUCLEOTIDE_MAP[nucleotide.upper()]
            # Increment the count for the corresponding nucleotide at position r,c
            profile[r][c] += 1
            
    consensus = get_consensus(profile, cols)

    return profile, consensus
        

if __name__ == "__main__":
    dna_sequence = """
    """
    
    with open('./data/rosalind_cons.txt', 'r') as file:
        dna_sequence = file.read()
    
    # Create a list with all the sequences and their ids
    ids, sequences, seq = [], [], ""

    for line in dna_sequence.split('\n'):
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
    # Max number of sequences = 10
    if len(ids) > 10:
        raise ValueError("The number of sequences must be 10 or less.")
    
    sequences_len = []
    for idx, sequence_list in enumerate(sequences):
        seq = ''.join(sequence_list)
        # Check length of each sequence
        if len(seq) > 1000:
            raise ValueError(f"Sequence {ids[idx]} is too long, must be 1000 characters or less.")
        # Check nucleotides
        if not all(nt in 'ACGT' for nt in seq.upper()):
            raise ValueError(f"Sequence {ids[idx]} contains invalid characters, only A, C, G, T are allowed.")
        sequences_len.append(len(seq))
    # Check if the length of sequences are the same
    if not all(x == sequences_len[0] for x in sequences_len):
        raise ValueError("Sequences must have the same length.")

    profile, consensus = get_profile(sequences)
    
    # Print the consensus string and profile matrix
    print(consensus)
    labels = list(NUCLEOTIDE_MAP.keys())
    for i, label in enumerate(labels):
        counts = ' '.join(map(str, profile[i]))
        print(f"{label}: {counts}")
    