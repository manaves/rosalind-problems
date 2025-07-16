"""
Given: A DNA string s of length at most 10 kbp, and a protein string translated by s.

Return: The index of the genetic code variant that was used for translation. (If multiple solutions exist, you may return any one.)
"""

from Bio.Seq import translate
# Excluded genetic code variants (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
EXCLUDED = [7, 8, 17, 18, 19, 20,]

def found_genetic_code_variant(dna, protein):
    """
    Finds the index of the genetic code variant used to translate a given DNA sequence into a protein.
    
    Parameters:
        dna (str): the DNA sequence to be translated.
        protein (str): the protein sequence that results from the translation of the DNA.
        
    Returns:
        int: the index of the genetic code variant used for translation.
    """
    
    for n in range(1,33):
        if n in EXCLUDED:
            continue
        # Translate the DNA sequence using the specified genetic code variant
        translation = translate(dna, table=n, stop_symbol='')
        
        # Check if the length of the translated sequence matches the protein sequence
        if len(translation) != len(protein):
            continue
        # Check if the translated sequence matches the protein sequence
        # If it matches, return the index of the genetic code variant
        if translation == protein:
            return n

if __name__ == "__main__":
    dna, protein = None, None
    
    with open('./Armory/data/rosalind_ptra.txt', 'r') as file:
        dna, protein = file.read().strip().split('\n')
        
    # Input validation
    if len(dna) > 10000:
        raise ValueError("The length of the DNA sequence must be less than or equal to 10 kpb.")
    if not all(nt in 'ACGT' for nt in dna.upper()):
        raise ValueError("DNA sequence contains invalid characters, only A, C, G, T are allowed.")
    
    gcv = str(found_genetic_code_variant(dna, protein))
    
    with open('./Armory/output/output_ptra.txt', 'w') as output_file:
        print(gcv)
        output_file.write(gcv)