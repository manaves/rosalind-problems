"""
A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string
    is the number of symbols that it contains.

An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s
"""

def count_nt(sequence):
    """
        Count the occurrences of each nucleotide in a DNA sequence.

        Parameters:
        sequence (str): a string representing a DNA sequence.

        Returns:
        dict: a dictionary with nucleotides as keys and their counts as values.
    """
    # Initialize a dictionary to hold the counts
    nt_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    # Iterate through each character in the sequence
    for nt in sequence.upper():
        if nt in nt_count:
            nt_count[nt] += 1

    return nt_count

# Example usage
if __name__ == "__main__":
    # Read DNA sequences from txt file
    dna_sequence = ""
    with open('./data/rosalind_dna.txt', 'r') as file:
        dna_sequence = file.read().strip()
    
    for sequence in dna_sequence.splitlines():
        # Validate the DNA sequence   
        if len(sequence) > 1000:
            raise ValueError("DNA sequence is too long, must be 1000 characters or less.")
        
        if not all(nt in 'ACGT' for nt in sequence.upper()):
            raise ValueError("DNA sequence contains invalid characters, only A, C, G, T are allowed.")
        
        # Count the nucleotides in the DNA sequence
        counts = count_nt(sequence)
        print(f"{counts['A']} {counts['C']} {counts['G']} {counts['T']} ")