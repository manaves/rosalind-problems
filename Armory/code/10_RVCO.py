"""
Recall that in a DNA string s, 'A' and 'T' are complements of each other, as are 'C' and 'G'. 
    Furthermore, the reverse complement of s is the string sc formed by reversing the symbols of s
    and then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

The Reverse Complement program from the SMS 2 package can be run online.

Given: A collection of n (n≤10) DNA strings.

Return: The number of given strings that match their reverse complements.
"""

from Bio import SeqIO

def check_palindromes(records):
    """
    Check if each DNA string in the records matches its reverse complement.
    
    Parameters:
        records (iterable): an iterable of SeqIO records (DNA strings in FASTA format
        
    Returns:
        result (int): the count of DNA strings that match their reverse complements.
    """
    
    result = 0
    # Calculate the reverse complement for each record
    for record in records:
        s = record.seq
        sc = s.reverse_complement()
        # Check if are equal
        if s == sc:
            result += 1
    
    return result

if __name__ == "__main__":
    records = list(SeqIO.parse('./Armory/data/rosalind_rvco.txt', 'fasta'))
    
    # Input validation
    if len(records) > 10:
        raise ValueError("The number of sequences must be less than or equal to 10.")
    
    result = check_palindromes(records)
    
    with open('./Armory/output/output_rvco.txt', 'w') as output_file:
        print(result)
        output_file.write(str(result))