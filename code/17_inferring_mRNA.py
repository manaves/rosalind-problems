"""
Given: A protein string of length at most 1000 aa.

Return: The total number of different RNA strings from which the protein could 
    have been translated, modulo 1,000,000. 
    (Don't neglect the importance of the stop codon in protein translation.)
"""

AA = { 'G': 4, 'A': 4, 'V': 4, 'L': 6, 'I': 3, 'P': 4, 
      'F': 2, 'Y': 2, 'W': 1, 'S': 6, 'T': 4, 'C': 2, 
      'M': 1, 'N': 2, 'Q': 2, 'K': 2, 'R': 6, 'H': 2, 
      'D': 2, 'E': 2, 'B': 4, 'Z': 4 }

MODULO_VALUE = 1000000

def inferring_mRNA(protein):
    """
    Calculate the number of different mRNA sequences that can produce a given protein sequence.
    
    Parameters:
        protein (str): the protein sequence.
        
    Returns:
        int: the total number of different mRNA sequences.
    """
    if not protein:
        return 0
    
    result = 1
    for aa in protein:
        if aa in AA.keys():
            result = (result * AA[aa]) % MODULO_VALUE
        else:
            raise ValueError(f"Invalid amino acid: {aa}")
    
    # Multiply by 3 for the stop codon
    result = (result * 3) % MODULO_VALUE
    
    return result

if __name__ == "__main__":
    protein = ""
    with open('./data/rosalind_mrna.txt', 'r') as file:
        protein = file.read().strip()
    
    # Validate input
    if not protein or any(aa not in AA for aa in protein):
        raise ValueError("Invalid input: protein sequence must contain valid amino acids.")
    
    if len(protein) > 1000:
        raise ValueError("Error: the length of the protein sequence must be less than or equal to 1000.")
    
    result = inferring_mRNA(protein)
    print(result)
    
    # Save result into txt
    with open('./output/rosalind_mrna.txt', 'w') as file:
        file.write(str(result))