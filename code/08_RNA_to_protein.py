"""
Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The protein string encoded by s.
"""

CODON_TABLE = {
    'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V',
    'UUC':'F', 'CUC':'L', 'AUC':'I', 'GUC':'V',
    'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V',
    'UUG':'L', 'CUG':'L', 'AUG':'M', 'GUG':'V',
    'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A',
    'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A',
    'UCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A',
    'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A',
    'UAU':'Y', 'CAU':'H', 'AAU':'N', 'GAU':'D',
    'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D',
    'UAA':'Stop','CAA':'Q', 'AAA':'K', 'GAA':'E',
    'UAG':'Stop','CAG':'Q', 'AAG':'K', 'GAG':'E',
    'UGU':'C', 'CGU':'R', 'AGU':'S', 'GGU':'G',
    'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
    'UGA':'Stop','CGA':'R', 'AGA':'R', 'GGA':'G',
    'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G'
}

def translation(s):
    """
    Function that translates an RNA string into a protein.
    
    Parameters:
        s (string): RNA string.
    
    Return:
        final_protein (string): protein string.
    """
    protein = []
    for p in range(0,len(s)-2,3):
        codon = s[p:(p+3)]
        aa = CODON_TABLE[codon]
        
        if aa == 'Stop':
            break
        else:
            protein.append(aa)
    
    final_protein = ''.join(protein)
    
    return final_protein        
        
if __name__ == "__main__":
    s = ""
    
    with open('./data/rosalind_prot.txt', 'r') as file:
        s = file.read().strip()
    
    # Input validation
    if len(s) > 10000:
        raise ValueError("Sequence is too long, must be 10000 characters or less.")
    if not all(nt in 'ACGU' for nt in s.upper()):
        raise ValueError("Sequence contains invalid characters, only A, C, G, U are allowed.")
    
    protein = translation(s)
    print(protein)