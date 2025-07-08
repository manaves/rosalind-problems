"""
Given: A DNA string s of length at most 1 kbp in FASTA format.

Return: Every distinct candidate protein string that can be translated from ORFs of s. 
    Strings can be returned in any order.
"""

CODON_TABLE = {
    'TTT':'F', 'CTT':'L', 'ATT':'I', 'GTT':'V',
    'TTC':'F', 'CTC':'L', 'ATC':'I', 'GTC':'V',
    'TTA':'L', 'CTA':'L', 'ATA':'I', 'GTA':'V',
    'TTG':'L', 'CTG':'L', 'ATG':'M', 'GTG':'V',
    'TCT':'S', 'CCT':'P', 'ACT':'T', 'GCT':'A',
    'TCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A',
    'TCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A',
    'TCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A',
    'TAT':'Y', 'CAT':'H', 'AAT':'N', 'GAT':'D',
    'TAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D',
    'TAA':'Stop','CAA':'Q', 'AAA':'K', 'GAA':'E',
    'TAG':'Stop','CAG':'Q', 'AAG':'K', 'GAG':'E',
    'TGT':'C', 'CGT':'R', 'AGT':'S', 'GGT':'G',
    'TGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
    'TGA':'Stop','CGA':'R', 'AGA':'R', 'GGA':'G',
    'TGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G'
}

NUCLEOTIDES = {
    'A': 'T', 'T': 'A',
    'G': 'C', 'C': 'G'
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
    
def find_orfs(sequence, l = 3):
    """
    Finds all open reading frames (ORFs) in a given DNA sequence.
    
    Parameters:
        sequence (str): a DNA sequence consisting of the characters A, T, C, and G.
        l (int): the length of the codon, default is 3.
        
    Returns:
        list: a list of distinct protein strings that can be translated from the ORFs.
    """
    
    proteins = set()
    # Iterate through the three possible reading frames
    for frame in range(l):
        s = sequence[frame:]
        start_codons_indices = []
        # Iterate through the sequence in codon-sized chunks
        for i in range(0, (len(s) - l + 1), l):
            # Extract the current codon
            codon = s[i:i+l]
            # Check if the codon is a start or stop codon
            if codon == 'ATG':
                start_codons_indices.append(i)
            elif codon in ['TAA', 'TAG', 'TGA']:
                # If it's a stop codon, translate all proteins from the start codons found so far
                for start_index in start_codons_indices:
                    protein = ""
                    for j in range(start_index, i, l):
                        current_codon = s[j:j+l]
                        aa = CODON_TABLE[current_codon]
                        protein += aa
                    # Ensure protein is not empty
                    if protein:
                        proteins.add(protein)
                # Clear start codons as we've processed them
                start_codons_indices = []

    return list(proteins)

if __name__ == "__main__":
    sequence = """"""
    
    with open('./data/rosalind_orf.txt', 'r') as file:
        sequence = file.read()
        
    # Get the sequence
    seq = ""

    for line in sequence.split('\n'):
        line = line.strip()
        if not line.startswith('>'):
            seq += line
            
    # Input validation
    if len(seq) > 1000:
        raise ValueError("The sequence length must be equal to or less than 1000 pb.")

    rev_com = reverse_complement(seq)
    
    proteins = find_orfs(seq)
    proteins_rev = find_orfs(rev_com)
    total_proteins = set(proteins+proteins_rev)
    with open('./output/output_orf.txt', 'w') as f:
        for protein in total_proteins:
            print(protein)
            f.write(f"{protein}\n")