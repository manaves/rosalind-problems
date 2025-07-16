"""
Given: FASTQ file

Return: Corresponding FASTA records
"""

from Bio import SeqIO

def convert_fastq_into_fasta(record):
    """
    Converts a FASTQ record into a FASTA format and writes it to an output file.
    
    Parameters:
        record (SeqIO.index): a SeqIO index of the FASTQ file.
        
    Returns:
        None: The function writes the FASTA records to an output file.
    """
    id_records = list(record)
    with open('./Armory/output/output_tfsq.txt', 'w') as output_file:
        for id in id_records:
            fasta = record[id].format('fasta')
            
            print(fasta)
            output_file.write(fasta)

if __name__ == "__main__":
    
    record = SeqIO.index('./Armory/data/rosalind_tfsq.txt', 'fastq')
   
    # Input validation
    if not record:
        raise ValueError("The FASTQ file is empty or not found.")
    
    convert_fastq_into_fasta(record)