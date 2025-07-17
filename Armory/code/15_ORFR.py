"""
Given: A DNA string s of length at most 1 kbp.

Return: The longest protein string that can be translated from an ORF of s. 
    If more than one protein string of maximum length exists, then you may output any solution.
"""

from Bio.Seq import Seq, translate
from re import finditer

def trim_to_codon(seq):
    """
    Trims the sequence to the nearest complete codon.
    
    Parameters:
        seq (Seq): the sequence to be trimmed.
        
    Returns:
        Seq: the trimmed sequence.
    """
    return seq[:(len(seq) - (len(seq) % 3))]

def find_orfs(orfs, s):
    """
    Finds open reading frames (ORFs) in a given sequence.
    
    Parameters:
        orfs (set): a set to store unique ORFs.
        s (Seq): the sequence in which to find ORFs.

    Returns:
        None: modifies the orfs set in place.
    """
    for i in finditer('ATG', str(s)):
        subseq = s[i.start():]
        trimmed = trim_to_codon(subseq)
        prot = translate(trimmed, table=1, stop_symbol='', to_stop=True)
        orfs.add(str(prot))

if __name__ == "__main__":
    with open('./Armory/data/rosalind_orfr.txt', 'r') as f:
        dna = Seq(f.read().strip())
        
    # Input validation
    if len(str(dna)) > 1000:
        raise ValueError("The length of the DNA sequence must be less than or equal to 1000.")
        
    # Calculate the reverse complement
    dna_rc = dna.reverse_complement()
    # Set to store unique ORFs
    orfs = set() 
    
    # Find ORFs in both strands
    # Direct strand
    find_orfs(orfs, dna)
    
    # Reverse complement strand
    find_orfs(orfs, dna_rc)
    
    # Find the longest ORF
    longest_orf = max(orfs, key=len)
    
    with open('./Armory/output/output_orfr.txt', 'w') as output_file:
        print(longest_orf)
        output_file.write(longest_orf)