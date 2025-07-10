"""
For this code, it's necessary to have MEME installed (https://meme-suite.org/meme/doc/install.html#installinggit).
Given: A set of protein strings in FASTA format that share some motif with minimum length 20.

Return: Regular expression for the best-scoring motif.
"""

from Bio import SeqIO
import xml.etree.ElementTree as ET
import subprocess
import os

def run_meme(file_path, output):
    """
    Run the MEME tool on the provided FASTA file and save the output to the specified directory.
    
    Parameters:
        file_path (str): path to the input FASTA file.
        output (str): path to the output directory where MEME results will be saved.
    
    Returns:
        None
    """
    try:
        # Options for MEME
        options = ["-protein", "-mod", "anr", "-nmotifs", "2", "-minw", "20"]
        command = ["meme", file_path, "-oc", output] + options
        # Run the MEME command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        # Check for errors
        if stderr:
                print("Error while running MEME:", stderr)
    except Exception as e:
        print("An error occurred:", e)

def get_best_motif_regex(file):
    """
    Find the motif with the lowest E-value in a MEME XML output file and return its consensus sequence.
    
    Parameters:
        file (str): path to the MEME XML output file.
    Returns:
        dict: a dictionary containing the motif ID, E-value, and consensus sequence of the best motif.
        If no valid motif is found, returns None.
    """
    lowest_e_value = float('inf')
    best_motif_info = None

    try:
        # Parse the XML file
        tree = ET.parse(file)
        root = tree.getroot()

        # Iterate through all motifs in the XML
        for motif in root.findall('.//motif'):
            motif_id = motif.get('id')
            
            # Searching for the E-value
            e_value_str = motif.get('e_value')
            if not e_value_str:
                # Search the e-value in the 'regular_expression' element
                # This is a common place where MEME stores the E-value for the motif
                e_value_element = motif.find('.//regular_expression')
                if e_value_element is not None:
                    e_value_str = e_value_element.get('e_value')
                
                # If not found, try to find it in the 'model' element
                if not e_value_str:
                    model_element = motif.find('.//model')
                    if model_element is not None:
                        e_value_str = model_element.get('e_value')

            if e_value_str:
                try:
                    e_value = float(e_value_str)
                except ValueError:
                    print(f"Warning: E-value for motif '{motif_id}' is not a valid float: {e_value_str}. Skipping this motif.")
                    continue
            else:
                print(f"Warning: No E-value found for motif '{motif_id}'. Skipping this motif.")
                continue

            # Searching for the consensus sequence
            # The consensus sequence is typically found in the 'regular_expression' element
            consensus_element = motif.find('regular_expression') 
            consensus_sequence = None
            if consensus_element is not None:
                # The consensus sequence is usually the second line in the text of the 'regular_expression' element
                consensus_sequence = consensus_element.text.split('\n')[1]
            
            if not consensus_sequence:
                pass

            if consensus_sequence and e_value < lowest_e_value:
                lowest_e_value = e_value
                best_motif_info = {
                    'motif_id': motif_id,
                    'e_value': e_value,
                    'consensus_sequence': consensus_sequence
                }

    except FileNotFoundError:
        print(f"Error: file '{file}' not found.")
        return None
    except ET.ParseError as e:
        print(f"Error parsing XML file: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

    return best_motif_info

if __name__ == "__main__":
    file_path = './Armory/data/rosalind_meme.txt'
    output_path = './Armory/output/meme'
    # Ensure the output directory exists, if not, create it
    os.makedirs(output_path, exist_ok=True)
    # Read the FASTA file    
    sequences = SeqIO.to_dict(SeqIO.parse(file_path, 'fasta'))
    
    # Input validation
    for id, information in sequences.items():
        if len(information.seq) < 20:
            raise ValueError("The length of the sequences must be upper than or equal to 20.")
    
    run_meme(file_path, output_path)
    output_file = os.path.join(output_path, 'meme.xml')
    motif = get_best_motif_regex(output_file)
   
    # Save the motif to a file
    with open('./Armory/output/output_meme.txt', 'w') as f:
        print(motif['consensus_sequence'])
        f.write(motif['consensus_sequence'])