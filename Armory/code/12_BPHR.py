"""
Given: FASTQ file, quality threshold q

Return: Number of positions where mean base quality falls below given threshold
"""

from Bio import SeqIO
from io import StringIO
import statistics

def average_bq(q, records):
    """
    Calculate the number of positions where the average base quality is below a given threshold.
    
    Parameters:
        q (int): the quality threshold.
        records (iterable): an iterable of SeqIO records (FASTQ entries).
        
    Returns:
        result (int): the count of positions where the average base quality is below q.
    """
    qualities, avg_qualities = [], []
    # Collect quality scores for each record
    for record in records:
        if 'phred_quality' not in record.letter_annotations:
            continue
        
        qualities.append(record.letter_annotations['phred_quality'])
    # Check if all sequences have the same length
    ref_len = len(qualities[0])
    
    if not all(len(quality) == ref_len for quality in qualities):
        raise ValueError("The length of the sequences must be the same.")
    # Calculate the average quality for each position
    for l in range(ref_len):
        avg = statistics.mean([quality[l] for quality in qualities])
        avg_qualities.append(avg)
    # Count the number of positions where the average quality is below q
    result = sum(num < q for num in avg_qualities)

    return result
        

if __name__ == "__main__":
    lines = None
    with open('./Armory/data/rosalind_bphr.txt', 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    q = int(lines[0]) 
    fastq_data = "\n".join(lines[1:])
    
    handle = StringIO(fastq_data)
    records = SeqIO.parse(handle, 'fastq')
    
    num = average_bq(q, records)
    
    with open('./Armory/output/output_bphr.txt', 'w') as output_file:
        print(num)
        output_file.write(str(num))
    