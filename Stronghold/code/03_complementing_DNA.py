"""
In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s, 
    then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s of length at most 1000 bp.

Return: The reverse complement sc of s
"""

def reverse_complement(dna_sequence):
    """
    Return the reverse complement of a DNA sequence.
    
    Parameters:
        dna_sequence (str): a string representing a DNA sequence.
        
    Returns:
        str: the reverse complement of the DNA sequence.
    """
    
    complement_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    
    # Reverse the DNA sequence and replace each nucleotide with its complement
    reversed_dna = reversed(dna_sequence)
    reverse_complemented_sequence = ''.join(complement_dict[nt] for nt in reversed_dna)
    
    return reverse_complemented_sequence

if __name__ == "__main__":
    dna_sequence = ""
    with open('./data/rosalind_revc.txt', 'r') as file:
        dna_sequence = file.read().strip()
        
    # Validate the DNA sequence
    if len(dna_sequence) > 1000:
        raise ValueError("DNA sequence is too long, must be 1000 characters or less.")
    
    if not all(nt in 'ACGT' for nt in dna_sequence.upper()):
        raise ValueError("DNA sequence contains invalid characters, only A, C, G, T are allowed.")
    
    reverse_complemented_dna = reverse_complement(dna_sequence)
    print(reverse_complemented_dna)