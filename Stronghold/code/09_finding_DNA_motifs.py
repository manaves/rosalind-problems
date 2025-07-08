"""
Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
"""

def finding_motifs(s, t):
    """
    Finds all occurrences of the motif t as a substring of s.
    
    Parameters:
        s (str): DNA sequence.
        t (str): DNA motif.
        
    Return:
        positions (str): s string containing the 1-based starting positions of 
            each occurrence of t in s, separated by spaces.
    """
    len_t = len(t)
    positions = []
    for m in range(len(s) - len_t + 1):
        # Extract a substring of s starting at position m with length equal to t
        motif = s[m:(m+len_t)]
        # If the substring matches t, add its 1-based position to the list
        if motif == t:
            positions.append(str(m+1))
    
    return ' '.join(positions)

if __name__ == "__main__":
    s = ""
    t = ""
    
    with open('./data/rosalind_subs.txt', 'r') as file:
        s, t = map(str, file.read().strip().split())
    
    # Input validation
    if len(s) > 1000 or len(t) > 1000:
        raise ValueError("Sequences are too long, must be 1000 characters or less.")
    if not all(nt in 'ACGT' for nt in s.upper()):
        raise ValueError("Sequence contains invalid characters, only A, C, G, T are allowed.")
    if len(t) > len(s):
        raise ValueError("The length of the motif must be shorter than that of the sequence.")
    
    positions = finding_motifs(s, t)
    print(positions)