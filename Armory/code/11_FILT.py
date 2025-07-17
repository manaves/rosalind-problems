from Bio import SeqIO
from io import StringIO

def filter_sequences(p, q, records):
    """
    Filter sequences based on the percentage of bases with quality scores below a threshold.
    
    Parameters:
        p (int): the percentage threshold for quality scores.
        q (int): the quality score threshold.
        records (iterable): an iterable of SeqIO records (FASTQ entries).
        
    Returns:
        result (int): the count of sequences where the percentage of bases with quality scores below q
                      is greater than or equal to p.
    """
    result = 0
    for record in records:
        if 'phred_quality' not in record.letter_annotations:
            continue
        # Get the quality score by base
        qualities = record.letter_annotations['phred_quality']
        # Count the number of bases with quality scores under q
        num_under_q = sum(num >= q for num in qualities)
        # Calculate the percentage of bases with quality scores under q
        perc = num_under_q/len(qualities)*100
        print(f"Perc: {perc}")
        # Check if the percentage meets or exceeds the threshold p
        if perc >= p:
            result += 1
           
    return result

if __name__ == "__main__":
    lines = None
    with open('./Armory/data/rosalind_filt.txt', 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    q, p = map(int, lines[0].split(' ')) 
    fastq_data = "\n".join(lines[1:])
    
    handle = StringIO(fastq_data)
    records = SeqIO.parse(handle, 'fastq')
    
    num = filter_sequences(p, q, records)
    
    with open('./Armory/output/output_filt.txt', 'w') as output_file:
        print(num)
        output_file.write(str(num))