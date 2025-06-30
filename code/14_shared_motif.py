"""
Given: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.

Return: A longest common substring of the collection. 

(If multiple solutions exist, you may return any single solution.)
"""

def find_shared_motif(sequences):
    """
    Finds the longest common substring (motif) shared among all DNA sequences.

    Parameters:
        sequences (list of str): a list of DNA sequences (strings composed of A, C, G, T).

    Returns:
        str: the longest common substring present in all sequences.
            If multiple substrings of the same maximum length exist, returns any one of them.
            Returns an empty string if no common substring exists.
    """
    
    # Select the shortest sequence as a reference.
    # No common substring can be longer than this sequence.
    shortest_seq = min(sequences, key=len)
    len_shortest = len(shortest_seq)

    # Check substrings from longest to shortest
    for length in range(len_shortest, 0, -1):
        # Slide a window of the current length across the shortest sequence
        for start in range(len_shortest - length + 1):
            motif = shortest_seq[start:(start+length)]

            # Check if this motif exists in all sequences
            if all(motif in seq for seq in sequences):
                return motif  # First match found will be the longest

    # If no common substring is found
    return ""

if __name__ == "__main__":
    dna_sequences = """"""
    
    with open('./data/rosalind_lcsm.txt', 'r') as file:
        dna_sequences = file.read()
        
    # Create a lists with all the sequences
    sequences, seq = [], ""

    for line in dna_sequences.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if seq:
                sequences.append(seq)
                seq = ""
        else:
            seq += line

    if seq:
        sequences.append(seq)
        
    # Input validation
    if len(sequences) > 100:
        raise ValueError("The number of sequences must be less than or equal to 100.")
    
    for sequence in sequences:
        if len(sequence) > 1000:
            raise ValueError("The length of the sequences must be less than or equal to 1000.")
        
    motif = find_shared_motif(sequences)
    print(motif)
    # Save motif to file
    with open('./output/output_lcsm.txt', 'w') as file:
        file.write(str(motif))
        
    