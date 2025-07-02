"""
Given: A DNA string of length at most 1 kbp in FASTA format.

Return: The position and length of every reverse palindrome 
    in the string having length between 4 and 12. 
    You may return these pairs in any order.
"""

import re

NUCLEOTIDES = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

def reverse_complement(sequence):
    """
    Returns the reverse complement of a given DNA sequence.
    
    Parameters:
        sequence (str): a DNA sequence consisting of the characters A, T, C, and G.
        
    Returns:
        str: the reverse complement of the input sequence.
    """
    c_seq = ''.join(NUCLEOTIDES[base] for base in reversed(sequence))
    return c_seq

def find_reverse_palindromes(seq, min_len=4, max_len=12):
    """
    Finds all reverse palindromes in a given DNA sequence within specified length range.
    
    Parameters:
        seq (str): a DNA sequence consisting of the characters A, T, C, and G.
        min_len (int): minimum length of the palindrome to find (default is 4).
        max_len (int): maximum length of the palindrome to find (default is 12).
        
    Returns:
        list: a list of tuples, each containing the starting position (1-based) and length of the reverse palindrome.
    """
    results = []
    # Iterate through lengths from min_len to max_len
    for l in range(min_len, max_len + 1):
        # Use a sliding window to check each fragment of length l
        for i in range(0, len(seq) - l + 1):
            # Extract the fragment of length l
            fragment = seq[i:i+l]
            # Check if the fragment is a reverse palindrome
            if fragment == reverse_complement(fragment):
                results.append((i + 1, l))
                
        results.sort(key=lambda x: x[0])
    return results

if __name__ == "__main__":
    sequence = """"""
    
    with open('./data/rosalind_revp.txt', 'r') as file:
        sequence = file.read()
        
    # Get the sequence
    seq = ""

    for line in sequence.split('\n'):
        line = line.strip()
        if not line.startswith('>'):
            seq += line
        
    # Input validation
    if len(seq) > 1000:
        raise ValueError("The DNA sequence length must be less than or equal to 1000 bp.")
    
    restriction_sites = find_reverse_palindromes(seq)
    
    # Output the result
    with open('./output/output_revp.txt', 'w') as file:
        for start, length in restriction_sites:
            file.write(f"{start} {length}\n")
        print('\n'.join(f"{start} {length}" for start, length in restriction_sites))