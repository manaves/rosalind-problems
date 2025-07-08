"""
For a collection of strings, a larger string containing every one of the smaller strings as a substring is called a superstring.
By the assumption of parsimony, a shortest possible superstring over a collection of reads serves as a candidate chromosome.

Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format 
    (which represent reads deriving from the same strand of a single linear chromosome).

The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct 
    the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).
"""

def find_overlap(seq_a, seq_b):
    """
    Finds the length of the maximum overlap between the end of seq_a and the start of seq_b.
    The overlap is defined as the longest suffix of seq_a that matches a prefix of seq_b
    
    Parameters:
        seq_a (str): the first sequence.
        seq_b (str): the second sequence.
        
    Returns:
        int: the length of the maximum overlap.
    """
    # Find the maximum overlap length between the end of seq_a and the start of seq_b
    max_len = min(len(seq_a), len(seq_b))
    # Check for overlaps from the maximum possible length down to 1
    for l in range(max_len, 0, -1):
        if seq_a[-l:] == seq_b[:l]:
            return l
    return 0

def genome_assembly(sequences):
    """
    Assembles a genome from a collection of sequences by merging them based on overlaps.
    
    Parameters:
        sequences (list): a list of DNA sequences.
    
    Returns:
        str: the assembled genome as a single string.
    """
    # While there are more than one sequence, keep merging them
    while len(sequences) > 1:
        # Initialize variables to track the maximum overlap and the best pair of sequences
        max_len = -1
        best_pair = (0, 0)
        best_merged = ""

        # Try every pair of sequences to find the pair with the maximum overlap
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if i != j:
                    # Find the overlap length between sequences[i] and sequences[j]
                    olen = find_overlap(sequences[i], sequences[j])
                    # If the overlap is greater than the current maximum, update the best pair and the best merged string
                    if olen > max_len:
                        max_len = olen
                        best_pair = (i, j)
                        # Create the merged string by combining the two sequences at the overlap
                        best_merged = sequences[i] + sequences[j][olen:]

        # Merge the best pair
        i, j = best_pair
        # Remove the merged sequences
        sequences.pop(max(j, i))
        sequences.pop(min(j, i))
        # Add the merged string
        sequences.append(best_merged)

    return sequences[0]

if __name__ == "__main__":
    inf = """"""
    
    with open('./data/rosalind_long.txt', 'r') as file:
        inf = file.read()
        
    # Create a list with all the sequences and their ids
    sequences, seq = [], ""

    for line in inf.split('\n'):
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
    if len(sequences) > 50:
        raise ValueError("The number of sequences must be less than or equal to 50.")
    
    for sequence in sequences:
        if len(sequence) > 1000:
            raise ValueError("The length of the sequences must be less than or equal to 1 kpb.")
    
    final_seq = genome_assembly(sequences)
    
    with open('./output/output_long.txt', 'w') as output_file:
        print(final_seq)
        output_file.write(final_seq)