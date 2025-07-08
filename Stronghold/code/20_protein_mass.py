"""
Given: A protein string P of length at most 1000 aa.

Return: The total weight of P. Consult the monoisotopic mass table.
"""

WEIGHTS = {
    'A':    71.03711,
    'C':    103.00919,
    'D':    115.02694,
    'E':    129.04259,
    'F':    147.06841,
    'G':    57.02146,
    'H':    137.05891,
    'I':    113.08406,
    'K':    128.09496,
    'L':    113.08406,
    'M':    131.04049,
    'N':    114.04293,
    'P':    97.05276,
    'Q':    128.05858,
    'R':    156.10111,
    'S':    87.03203,
    'T':    101.04768,
    'V':    99.06841,
    'W':    186.07931,
    'Y':    163.06333,
    'H2O':  18.01056
}

def calculate_prot_weight(protein):
    """
    Calculate the total weight of a protein string based on the monoisotopic mass table.
    
    Parameters:
        protein (str): a protein string consisting of amino acids represented by their one-letter codes.
        
    Returns:
        float: the total weight of the protein string.
    """
    
    result = 0
    
    for aa in protein:
        result += WEIGHTS[aa]
    
    return result

if __name__ == "__main__":
    protein = ""
    
    with open('./data/rosalind_prtm.txt', 'r') as file:
        protein = file.read().strip()
        
    # Input validation
    if len(protein) > 1000:
        raise ValueError("The protein length must be less than or equal to 1000 aa.")
    
    w = calculate_prot_weight(protein)
    
    with open('./output/output_prtm.txt', 'w') as f:
        print(w)
        f.write(f"{w}\n")