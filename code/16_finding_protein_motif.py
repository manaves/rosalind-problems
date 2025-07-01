"""
Given: At most 15 UniProt Protein Database access IDs.

Return: For each protein possessing the N-glycosylation motif, output its given access ID 
    followed by a list of locations in the protein string where the motif can be found.
"""

import urllib.request 
import re

def retrieve_sequences(ids, output='./data/rosalind_mprt_seqs.txt'):
    """
    Retrieve protein sequences from UniProt based on given IDs and save them to a file.
    
    Parameters:
        ids (list): a list of UniProt IDs.
        output (str): the path to the output file where sequences will be saved.
        
    Returns:
        list: a list of sequences corresponding to the given IDs.
    """
        
    sequences = []
    with open(output, 'w') as f:
        for id in ids:
            id = id.split('_')[0]
            actual_url = f"http://www.uniprot.org/uniprot/{id}.fasta"
            try:
                with urllib.request.urlopen(actual_url) as response:
                    # Read and write sequences
                    fasta_data = response.read().decode('utf-8')
                    f.write(fasta_data + '\n')
                    
                    # Extract raw sequence (remove header lines)
                    lines = fasta_data.strip().split('\n')
                    seq = ''.join(line for line in lines if not line.startswith('>'))
                    sequences.append(seq)
            except Exception as e:
                print(f"Failed to retrieve {id}: {e}")
                sequences.append(None)
    
    return sequences
def finding_motif(sequence, pattern = r'(?=(N[^P][ST][^P]))'):
    """
    Find all occurrences of the N-glycosylation motif in a protein sequence.
    
    Parameters:
        sequence (str): the protein sequence to search.
        pattern (str): the regex pattern for the N-glycosylation motif.
        
    Returns:
        list: a list of positions (1-based index) where the motif is found.
    """
    if not sequence:
        return []    
    
    # Find all matches of the pattern in the sequence
    pattern = re.compile(pattern)

    matches = [m.start() + 1 for m in pattern.finditer(sequence)]
    return matches

if __name__ == "__main__":
    ids = []
    with open('./data/rosalind_mprt.txt', 'r') as file:
        ids = file.read().strip().split('\n')
        
    # Input validation
    if len(ids) > 15:
        raise ValueError("Error: the number of UniProt IDs must be less than or equal to 15.")
    
    # Get sequences
    sequences = retrieve_sequences(ids)
    
    # Output file
    with open('./output/output_mprt.txt', 'w') as f:
        for idx, id in enumerate(ids):
            sequence = sequences[idx]
            if sequence is None:
                continue
            positions = finding_motif(sequence)
            if positions:
                print(id)
                print(' '.join(map(str, positions)))
                f.write(f"{id}\n")
                f.write(' '.join(map(str, positions)) + '\n')