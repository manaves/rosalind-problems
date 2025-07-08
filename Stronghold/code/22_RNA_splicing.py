"""
Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s 
    acting as introns. All strings are given in FASTA format.

Return: A protein string resulting from transcribing and translating the exons of s.
    (Note: Only one solution will exist for the dataset provided.)
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

def remove_introns(sequence, introns):
    """
    Removes introns from the given DNA sequence.
    
    Parameters:
        sequence (str): the original DNA sequence from which introns will be removed.
        introns (list): a list of intron sequences to be removed from the original sequence
    
    Returns:
        str: the cleaned DNA sequence with all introns removed.
    """
    clean_sequence = sequence
    
    for intron in introns:
        # Replace each intron in the sequence with an empty string
        clean_sequence = clean_sequence.replace(intron, '')
    
    return clean_sequence

def traduction(clean_sequence, l = 3):
    """
    Translates a cleaned DNA sequence into a protein string.
    
    Parameters:
        clean_sequence (str): the cleaned DNA sequence without introns.
        l (int): the length of each codon (default is 3).
        
    Returns:
        str: the resulting protein string from the translation of the cleaned sequence.
    """
    protein = ""
    
    for n in range(0, (len(clean_sequence) - l + 1), l):
        # Extract the codon of length l starting from position n
        curr_codon = clean_sequence[n:(n+l)]
        # If the codon is a stop codon, break the loop
        if CODON_TABLE[curr_codon] == 'Stop':
            break
        # Append the corresponding amino acid to the protein string
        protein += CODON_TABLE[curr_codon]
        
    return protein

if __name__ == "__main__":
    dna_sequence = """
    """
    
    with open('./data/rosalind_splc.txt', 'r') as file:
        dna_sequence = file.read()
    
    # Create a list with all the sequences and their ids
    sequences, seq = [], ""

    for line in dna_sequence.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if seq:
                sequences.append(seq)
                seq = ""
        else:
            seq += line

    if seq:
        sequences.append(seq)
        
    sequence, introns = sequences[0], sequences[1:]
    
    # Input validation
    if len(sequence) > 1000:
        raise ValueError("The length of the input string must be less than or equal to 1000.")
    
    clean_sequence = remove_introns(sequence, introns)
    protein = traduction(clean_sequence)
    
    with open('./output/output_splc.txt', 'w') as output_file:
        print(protein)
        output_file.write(str(protein))