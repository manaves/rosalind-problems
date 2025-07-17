"""
Given: Set of nucleotide strings in FASTA format.

Return: ID of the string most different from the others.
"""

import subprocess
from Bio import AlignIO
from itertools import combinations

def hamm(s1, s2):
    """
    Computes the Hamming distance between two strings of equal length.
    
    Parameters:
        s1 (str): the first string.
        s2 (str): the second string.
        
    Returns:
        int: the Hamming distance between the two strings.
    """
    
    count = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
            
    return count

def run_clustal():
    """
    Run ClustalW to align the sequences in the FASTA file.
    
    Returns:
        alignment (MultipleSeqAlignment): the aligned sequences.
    """
    # Run ClustalW command
    command = "clustalw -INFILE=./Armory/data/rosalind_clus.txt -OUTFILE=./Armory/output/output_clus_aligned.fasta -OUTPUT=FASTA"
    try:
        subprocess.run(command, shell=True, check=True)
        print("ClustalW alignment completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running ClustalW: {e}")
    
    # Read the alignment
    alignment = AlignIO.read("./Armory/output/output_clus_aligned.fasta", "fasta")
    
    return alignment

if __name__ == "__main__":
    # Run ClustalW to get the alignment
    alignment = run_clustal()
    
    # Calculate pairwise distances
    distances = {}
    for seq1, seq2 in combinations(alignment, 2):
        dist = hamm(seq1.seq, seq2.seq)
        # Store the distance in both directions
        distances[(seq1.id, seq2.id)] = dist
    
    # Find the sequence with the maximum distance from all others
    max_distance = -1
    most_different_id = None
    
    for seq in alignment:
        total_distance = sum(distances.get((seq.id, other.id), 0) + distances.get((other.id, seq.id), 0) for other in alignment if other.id != seq.id)
        if total_distance > max_distance:
            max_distance = total_distance
            most_different_id = seq.id
    
    # Output the ID of the most different sequence
    with open('./Armory/output/output_clus.txt', 'w') as output_file:
        print(most_different_id)
        output_file.write(most_different_id)