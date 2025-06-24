"""
Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: 
    k individuals are homozygous dominant for a factor, m  are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a 
    dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
"""

def mendels_first_law(m, n, k):
    """
    Calculate the probability that two randomly selected mating organisms will produce an individual possessing a 
        dominant allele.
    
    Parameters:
        m (int): number of heterozygous.
        n (int): number of homozygous recessive.
        k (int): number of homozygous dominant.
        
    Return:
        p_dominant (float): probability that two randomly selected mating organisms will produce an individual possessing a 
            dominant allele.
    """
    
    # Total of organisms
    total = m+n+k
    # Calculation of the probability of being recessive
    p_recessive = ((m*(m-1)*1/4) + m*n + n*(n-1))/(total*(total-1))
    # Calculation of the probability of being dominant
    p_dominant = 1 - p_recessive
    
    return p_dominant

if __name__ == "__main__":
    m = 0
    n = 0
    k = 0
    
    with open('./data/rosalind_iprb.txt', 'r') as file:
        k, m, n = map(int, file.read().strip().split())
    
    # Input validation
    if m < 0 or n < 0 or k < 0:
        raise ValueError("Invalid input: m, n and k must be positive numbers.")
    
    result = mendels_first_law(m, n, k)
    print(result)