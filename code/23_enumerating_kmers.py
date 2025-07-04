"""
Given: A collection of at most 10 symbols defining an ordered alphabet, 
    and a positive integer n (n≤10).

Return: All strings of length n that can be formed from the alphabet, 
    ordered lexicographically (use the standard order of symbols in the English alphabet).
"""

from itertools import product

def enumerate_kmers(alphabet, n):
    """
    Enumerates all possible kmers of length n from the given alphabet.
    
    Parameters:
        alphabet (list): a list of symbols representing the alphabet.
        n (int): the length of the kmers to generate.
        
    Returns:
        list: a sorted list of all kmers of length n formed from the alphabet.
    """
    # Generate all combinations with replacement of the alphabet of length n        
    kmers = product(alphabet, repeat=n)
    
    # Join each tuple to form strings
    kmers_strings = [''.join(kmer) for kmer in kmers]
    
    # Sort the kmers lexicographically
    kmers_strings.sort()
    
    return kmers_strings

if __name__ == "__main__":
    alphabet = ""
    n = 0
    
    with open('./data/rosalind_lexf.txt', 'r') as file:
        alphabet, n = file.read().strip().split('\n')
        n = int(n)
        
    # Get the alphabet
    alphabet = alphabet.split(' ')
    
    # Input validation
    if len(alphabet) > 10:
        raise ValueError("The length of the input string must be less than or equal to 10.")
    
    if n > 10 and n > 0:
        raise ValueError("n must be between 1 and 10 (both included).")
    # Enumerate kmers
    kmers = enumerate_kmers(alphabet, n)
    
    # Output the kmers
    with open('./output/output_lexf.txt', 'w') as output_file:
        for kmer in kmers:
            output_file.write(kmer + '\n')
            print(kmer)