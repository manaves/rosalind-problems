"""
Given: Two DNA strings s and t in FASTA format that share some short inexact repeat r
    of 32-40 bp. By "inexact" we mean that r may appear with slight modifications 
    (each repeat differ by ≤3 changes/indels).

Return: The total number of occurrences of r as a substring of s, followed by the 
    total number of occurrences of r as a substring of t.
    
IMPORTANT: it's necessary to use the LALIGN algorithm. You can download ir from:
    https://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml
"""

from Bio import SeqIO
import subprocess

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


def count_hamming(pattern, seq, dist=3):
    """
    Counts the number of occurrences of a pattern in a sequence with a given Hamming distance.
    
    Parameters:
        pattern (str): the pattern to search for.
        seq (str): the sequence in which to search for the pattern.
        dist (int): the maximum allowed Hamming distance.
        
    Returns:
        int: the count of occurrences of the pattern in the sequence within the specified Hamming distance.
    """
    
    count = 0
    np = len(pattern)
    ns = len(seq)
    # Iterate through the sequence to find occurrences of the pattern
    for i in range(ns - np + 1):
        if hamm(seq[i : i + np], pattern) <= dist:
            count += 1
            
    return count

def run_lalign():
    """
    Runs the LALIGN algorithm to find the longest common subsequence between two sequences.
    
    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        
    Returns:
        str: The longest common subsequence.
    """
    
    command = '../../fasta-36.3.8i/bin/lalign36 -m 9i ./Armory/data/seq1.fasta ./Armory/data/seq2.fasta -n -r +4/-8 -f -8 -g -8 -E 0.000000001 -m 0 -O "./Armory/output/output_subo_laling.lav"'
    with open('./Armory/output/output_subo_laling.lav', 'w') as output_file:
        subprocess.run(command, shell=True, stdout=output_file, stderr=subprocess.PIPE)    
    
import re

def extract_target_alignments(file_path, min_len=32, max_len=40):
    """
    Extracts target alignments from a LALIGN output file based on specified length criteria.
    
    Parameters:
        file_path (str): the path to the LALIGN output file.
        min_len (int): the minimum length of the alignment to consider.
        max_len (int): the maximum length of the alignment to consider.
        
    Returns:
        list: a list of target alignments that meet the length criteria.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    collecting = False
    alignment_block = set()
    results = set()

    for i, line in enumerate(lines):
        line = line.strip()

        # Check for the start of a new alignment block
        if "identity" in line and "overlap" in line:
            match = re.search(r"in (\d+) nt overlap", line)
            if match:
                # Extract the alignment length from the line
                aln_len = int(match.group(1))
                if min_len <= aln_len <= max_len:
                    for j in range(i + 1, len(lines)):
                        # Collect lines until the end of the alignment block
                        actual_line = lines[j].strip()
                        # Check for the end of the alignment block
                        if lines[j].startswith(">--") or lines[j].startswith(">>"):
                            break
                        else:
                            # Collect lines that start with "Rosali" to form the alignment block
                            if actual_line.startswith("Rosali"):
                                collecting = True
                                alignment_block.add(actual_line)
                else:
                    collecting = False
                    alignment_block = set()
        elif collecting:
            alignment_block_act = set(alignment_block)
            for alignment in alignment_block_act:
                # Check if the alignment starts with "Rosali" and add it to results
                if alignment.startswith("Rosali"):
                    ali = alignment.split(" ")[1]
                    results.add(ali)

    return list(results)[0]
    

if __name__ == "__main__":
    records = SeqIO.index('./Armory/data/rosalind_subo.txt', 'fasta')
    # Get the ids
    ids = list(records.keys())
    # Save sequences fo FASTA files
    for idx,record in enumerate(records.values()):
        SeqIO.write(record, f'./Armory/data/seq{idx+1}.fasta', 'fasta')
    
    # Input validation
    if len(ids) != 2:
        raise ValueError("The FASTA file must contain exactly two sequences.")

    sequences = [records[ids[0]].seq, records[ids[1]].seq]
    
    run_lalign()
    target_alignments = extract_target_alignments('./Armory/output/output_subo_laling.lav')
    
    with open('./Armory/output/output_subo.txt', 'w') as output_file:
        for idx, sequence in enumerate(sequences):
            count = count_hamming(target_alignments, sequence, dist=3)
            print(count)
            if idx == 0:
                output_file.write(f"{str(count)} ")
            else:
                output_file.write(f"{str(count)}")
        
    print(target_alignments)