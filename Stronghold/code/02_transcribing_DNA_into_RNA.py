"""
An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by
    replacing all occurrences of 'T' in t with 'U' in u.

Given: A DNA string t having length at most 1000 nt.

Return: The transcribed RNA string of t.
"""

def transcribe_dna_to_rna(dna_sequence):
    """
    Transcribe a DNA sequence into an RNA sequence by replacing 'T' with 'U'.
    
    Parameters:
        dna_sequence (str): a string representing a DNA sequence.
        
    Returns:
        str: the transcribed RNA sequence.
    """
    return dna_sequence.replace('T', 'U')

if __name__ == "__main__":
    dna_sequence = ""
    with open('./data/rosalind_rna.txt', 'r') as file:
        dna_sequence = file.read().strip()
    
    # Validate the DNA sequence
    if len(dna_sequence) > 1000:
        raise ValueError("DNA sequence is too long, must be 1000 characters or less.")
    
    if not all(nt in 'ACGT' for nt in dna_sequence.upper()):
        raise ValueError("DNA sequence contains invalid characters, only A, C, G, T are allowed.")
    
    transcribed_dna = transcribe_dna_to_rna(dna_sequence)
    print(transcribed_dna)