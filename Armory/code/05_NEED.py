"""
Given: Two GenBank IDs.

Return: The maximum global alignment score between the DNA strings associated with these IDs.
"""

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

def retrieve_records(ids):
    """
    Retrieves the GenBank records for the given IDs and returns them as a list of SeqRecord objects.
    
    Parameters:
        ids (str): a comma-separated string of GenBank IDs.
    
    Returns:
        list: a list of SeqRecord objects containing the sequences.
    """
    
    Entrez.email = 'manaves@uoc.edu'
    handle = Entrez.efetch(db="nucleotide", id=[ids], rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    
    return records

def score_alignment(records, match=5, mismatch=-4, op_gap=-10, ex_gap=-1):
    """
    Calculates the maximum global alignment score between two sequences.
    
    Parameters:
        records (list): a list of SeqRecord objects containing the sequences to align.
        match (int): score for a match.
        subst (int): score for a mismatch.
        op_gap (int): score for opening a gap.
        ex_gap (int): score for extending a gap.
        
    Returns:
        int: the maximum global alignment score.    
    """
    seq1 = records[0].seq
    seq2 = records[1].seq
    
    # Create a pairwise aligner with the specified scoring parameters
    aligner = PairwiseAligner()
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = op_gap
    aligner.extend_gap_score = ex_gap
    
    # Perform the alignment
    aligner.mode = 'global'
    alignment = aligner.align(seq1, seq2)
    
    score = int(alignment.score)
    
    return score


if __name__ == "__main__":
    genbank_ids = None
    with open("./Armory/data/rosalind_need.txt", "r") as f:
        genbank_ids = ", ".join(f.readline().strip().split())
    
    # Input validation
    if not genbank_ids:
        raise ValueError("GenBank IDs cannot be empty.")
    if not isinstance(genbank_ids, str):
        raise TypeError("GenBank IDs must be a string.")
    if len(genbank_ids.split(',')) != 2:
        raise ValueError("The number of GenBank IDs must be equal to 2.")

    records = retrieve_records(genbank_ids)
    score = str(score_alignment(records=records))
    
    with open('./Armory/output/output_need.txt', 'w') as output_file:
        print(score)
        output_file.write(score)