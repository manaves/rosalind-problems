"""
Given: A DNA string s of length at most 1000 bp.

Return: Four integers (separated by spaces) representing the respective number 
    of times that the symbols 'A', 'C', 'G', and 'T' occur in s. Note: You must 
    provide your answer in the format shown in the sample output below.
"""

from Bio.Seq import Seq

NUCLEOTIDES = ['A', 'C', 'G', 'T']

def count_nucleotides(sequence):
    """
    Counts the occurrences of each nucleotide in the given DNA sequence.
    
    Args:
        sequence (Seq): a DNA sequence.
        
    Returns:
        list: a list of four integers representing the counts of 'A', 'C', 'G', and 'T' in that order.
    """
    results = []
    for nt in NUCLEOTIDES:
        count = sequence.count(nt)
        results.append(count)
    
    return results

if __name__ == "__main__":
    sequence = None
    
    with open('./Armory/data/rosalind_ini.txt', 'r') as file:
        sequence = Seq(file.read().strip()).upper()
    
    # Input validation
    if len(sequence) > 1000:
        raise ValueError("The length of the sequence must be less than or equal to 1000 pb.")
    
    if not all(nt in NUCLEOTIDES for nt in sequence):
        raise ValueError("The sequence must only contain the nucleotides A, C, G, and T.")
    
    counts = count_nucleotides(sequence=sequence)
    
    with open('./Armory/output/output_ini.txt', 'w') as output_file:
        result = ' '.join(map(str, counts))
        print(result)
        output_file.write(result)