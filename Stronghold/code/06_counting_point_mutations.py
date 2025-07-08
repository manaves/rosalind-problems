"""
Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), 
is the number of corresponding symbols that differ in s and t.

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t).
"""

def hamming_distance(s, t):
    differences = 0
    for c in range(len(s)):
        if s[c] != t[c]:
            differences += 1
            
    return differences

if __name__ == "__main__":
    s = ""
    t = ""
    with open("data/rosalind_hamm.txt", "r") as file:
        s = file.readline().strip()
        t = file.readline().strip()
    
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length")
    
    if len(s) > 1000 or len(t) > 1000:
        raise ValueError("Strings must not exceed 1 kbp")
    
    if not all(nt in 'ACGT' for nt in s.upper()) or not all(nt in 'ACGT' for nt in t.upper()):
        raise ValueError(f"Sequence contains invalid characters, only A, C, G, T are allowed.")
    
    # Calculate the Hamming distance
    distance = hamming_distance(s, t)
    
    print(distance)