"""
To run this code it's necessary to have Trimmomatic installed, and you have to change the path.
Given: FASTQ file, quality cut-off value q, Phred33 quality score assumed.

Return: FASTQ file trimmed from the both ends (removed leading and trailing bases with quality lower than q)
"""

import subprocess

def run_trimmomatic(q):
    """
    Run Trimmomatic to trim low-quality bases from a FASTQ file.
    
    Parameters:
        q (int): the quality threshold for trimming.
        
    Returns:
        None
    """
    command = f"java -jar ../../trimmomatic/trimmomatic-0.39/dist/jar/trimmomatic-0.39.jar SE -phred33 ./Armory/data/rosalind_bfil_fastq.fastq ./Armory/output/output_bfil.fastq LEADING:{q} TRAILING:{q}"
    
    try:
        subprocess.run(command, shell=True, check=True)
        print("Trimmomatic completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running Trimmomatic: {e}")
        
if __name__ == "__main__":
    # Run txt file to get the quality threshold and fastq file
    with open('./Armory/data/rosalind_bfil.txt', 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    q = int(lines[0])  # Quality threshold
    # Save the FASTQ data to a file
    fastq_data = "\n".join(lines[1:])
    with open('./Armory/data/rosalind_bfil_fastq.fastq', 'w') as fastq_file:
        fastq_file.write(fastq_data)
    # Run Trimmomatic with the specified quality threshold
    run_trimmomatic(q)