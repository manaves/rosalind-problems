"""
Given: A quality threshold, along with FASTQ entries for multiple reads.

Return: The number of reads whose average quality is below the threshold.
"""

from Bio import SeqIO
from io import StringIO

def quality_control(records, threshold):
    """
    Counts the number of reads with an average quality below a given threshold.
    
    Parameters:
        records (iterable): an iterable of SeqIO records (FASTQ entries).
        threshold (int): the quality threshold to compare against.
    
    Returns:
        num (int): the number of reads with an average quality below the threshold.
    """
    num = 0
    # Iterate through each record in the FASTQ data
    for record in records:
        # Check if the record has letter annotations for quality
        if 'phred_quality' not in record.letter_annotations:
            continue
        # Calculate the average quality score
        qualities = record.letter_annotations['phred_quality']
        avg_quality = sum(qualities) / len(qualities)
        # If the average quality is below the threshold, increment the count
        if avg_quality < threshold:
            num += 1
    return num

if __name__ == "__main__":
    lines = None
    with open('./Armory/data/rosalind_phre.txt', 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    threshold = int(lines[0]) 
    fastq_data = "\n".join(lines[1:])
    
    handle = StringIO(fastq_data)
    records = SeqIO.parse(handle, 'fastq')
    
    num = quality_control(records, threshold)
    
    with open('./Armory/output/output_phre.txt', 'w') as output_file:
        print(num)
        output_file.write(str(num))