"""
Given: A collection of n (n≤10) GenBank entry IDs.

Return: The shortest of the strings associated with the IDs in FASTA format.

Example:
    Input: FJ817486 JX069768 JX469983
    Output:
        >JX469983.1 Zea mays subsp. mays clone UT3343 G2-like transcription factor mRNA, partial cds
        ATGATGTATCATGCGAAGAATTTTTCTGTGCCCTTTGCTCCGCAGAGGGCACAGGATAATGAGCATGCAA
        GTAATATTGGAGGTATTGGTGGACCCAACATAAGCAACCCTGCTAATCCTGTAGGAAGTGGGAAACAACG
        GCTACGGTGGACATCGGATCTTCATAATCGCTTTGTGGATGCCATCGCCCAGCTTGGTGGACCAGACAGA
        GCTACACCTAAAGGGGTTCTCACTGTGATGGGTGTACCAGGGATCACAATTTATCATGTGAAGAGCCATC
        TGCAGAAGTATCGCCTTGCAAAGTATATACCCGACTCTCCTGCTGAAGGTTCCAAGGACGAAAAGAAAGA
        TTCGAGTGATTCCCTCTCGAACACGGATTCGGCACCAGGATTGCAAATCAATGAGGCACTAAAGATGCAA
        ATGGAGGTTCAGAAGCGACTACATGAGCAACTCGAGGTTCAAAGACAACTGCAACTAAGAATTGAAGCAC
        AAGGAAGATACTTGCAGATGATCATTGAGGAGCAACAAAAGCTTGGTGGATCAATTAAGGCTTCTGAGGA
        TCAGAAGCTTTCTGATTCACCTCCAAGCTTAGATGACTACCCAGAGAGCATGCAACCTTCTCCCAAGAAA
        CCAAGGATAGACGCATTATCACCAGATTCAGAGCGCGATACAACACAACCTGAATTCGAATCCCATTTGA
        TCGGTCCGTGGGATCACGGCATTGCATTCCCAGTGGAGGAGTTCAAAGCAGGCCCTGCTATGAGCAAGTC
        A
"""

from Bio import SeqIO, Entrez

def search_shortest(records):
    lengths = []
    for record in records:
        sequence = record.seq
        lengths.append(len(sequence))
    idx_min = lengths.index(min(lengths))
    
    min_record = records[idx_min]
    
    return min_record

def get_sequences(ids):
    Entrez.email = "manaves@uoc.edu"
    
    # Construct the search term with the genus name and date range
    term = ""
    for id in ids:
        if id == ids[-1]:
            term += f'{id}'
        else:
            term += f'{id},'
    search_term = [term]
    handle = Entrez.efetch(db='nucleotide', id=search_term, rettype='fasta')
    records = list (SeqIO.parse(handle, "fasta"))
    handle.close()
    
    return records

if __name__ == "__main__":
    ids = []
    
    with open('./Armory/data/rosalind_frmt.txt', 'r') as file:
        ids = file.read().strip().split(' ')
    
    # Input validation
    if len(ids) > 10:
        raise ValueError("The number of IDs must be less than or equal to 10.")
    
    records = get_sequences(ids)
    result = search_shortest(records)
    
    SeqIO.write(result, handle='./Armory/output/output_frmt.txt', format='fasta')
    
    