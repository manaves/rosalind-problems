"""
The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. 
    For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA 
    string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string 
    labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', 
    followed by some labeling information. Subsequent lines contain the string itself; the first line to begin 
    with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" 
    denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. 
    Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated.
"""


def gc_content(sequence):
    """
    Calculate the GC content of a given DNA sequence.
    
    Parameters:
    sequence (str): A string representing the DNA sequence (only contains 'A', 'T', 'C', 'G').
    
    Returns:
    float: The GC content as a percentage of the total bases in the sequence.
    """
    if not sequence or len(sequence) == 0:
        return 0.0
    
    gc_count = sum(1 for base in sequence if base in 'GC')
    length = len(sequence)
    
    return (gc_count / length) * 100

if __name__ == "__main__":
    dna_sequence = """
    """
    with open('./data/rosalind_gc.txt', 'r') as file:
        dna_sequence = file.read()
        
    # Run all the FASTa sequences through the gc_content function
    sequences = {}
    current_id = None
    for line in dna_sequence.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_id is not None:
                sequences[current_id] = ''.join(sequences[current_id])
            current_id = line[1:]
            sequences[current_id] = []
        else:
            if current_id is not None:
                sequences[current_id].append(line)
        
    if current_id is not None:
        sequences[current_id] = ''.join(sequences[current_id])
        
    # Calculate GC content for each sequence and validate the input
    gc_contents = {}
    for seq_id, seq in sequences.items():
        if len(seq) > 1000:
            raise ValueError(f"Sequence {seq_id} is too long, must be 1000 characters or less.")
        if not all(nt in 'ACGT' for nt in seq.upper()):
            raise ValueError(f"Sequence {seq_id} contains invalid characters, only A, C, G, T are allowed.")
        
        # Calculate GC content for each sequence
        gc_contents[seq_id] = gc_content(seq)
        
    # Find the sequence with the highest GC content
    max_gc_id = max(gc_contents, key=gc_contents.get)
    max_gc_content = gc_contents[max_gc_id]
    
    # Print the result
    print(f"{max_gc_id}\n{max_gc_content}")