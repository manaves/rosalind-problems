"""
Given: Two positive integers k (k≤7) and N (N≤2k). In this problem, we begin with Tom, 
    who in the 0th generation has genotype Aa Bb. Tom has two children in the 1st generation, 
    each of whom has two children, and so on. Each organism always mates with an organism 
    having genotype Aa Bb.

Return: The probability that at least N Aa Bb organisms will belong to the k-th generation of 
    Tom's family tree (don't count the Aa Bb mates at each level). Assume that Mendel's second 
    law holds for the factors.
"""

def factorial(n):
    """
    Calculate the factorial of a number n.
    
    Parameters:
        n (int): the number to calculate the factorial for.
        
    Returns:
        int: the factorial of n.
    """
    fact = 1
    for num in range(2, n + 1):
        fact *= num
    return fact

def independent_alleles(k, N):
    """
    Calculate the probability of having at least N AaBb organisms in the k-th generation.
    
    Parameters:
        k (int): the generation number (k ≤ 7).
        N (int): the minimum number of AaBb organisms (N ≤ 2^k).
        
    Returns:
        float: the probability of having at least N AaBb organisms in the k-th generation.
    """
    # Num of organisms
    n = 2**k
    # AaBb probability
    p = 1/4
    result = 0
    
    for i in range(N, n+1):
        # Binomial distribution
        result += ((factorial(n)/(factorial(i)*factorial(n-i)))) * (p ** i) * ((1-p) ** (n - i))
    return result

if __name__ == "__main__":
    k = 0
    N = 0
    
    with open('./data/rosalind_lia.txt', 'r') as file:
        k, N = map(int, file.read().strip().split())
        
    # Input validation
    if N > 2**k:
        raise ValueError("Error: N must be equal to or less than 2^k.")
    
    if k > 7:
        raise ValueError("Error: k must be equal to or less than 7.")
    
    result = independent_alleles(k,N)
    print(result)
    
    # Save result into txt
    with open('./output/rosalind_lia.txt', 'w') as file:
        file.write(str(result))