"""
Given: A positive integer n≤7.

Return: The total number of permutations of length n, followed by 
    a list of all such permutations (in any order).
"""

from itertools import permutations as perm

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

def permutations(n):
    """
    Calculate the total number of permutations of length n and return a list of all such permutations.
    
    Parameters:
        n (int): the length of the permutations (n ≤ 7).
        
    Returns:
        tuple: a tuple containing the total number of permutations and a list of all permutations.
    """
    
    # Generate all permutations of length n
    perm_list = list(perm(range(1, n + 1)))
    total_permutations = len(perm_list)
    
    return total_permutations, perm_list

if __name__ == "__main__":
    n = 0
    
    with open('./data/rosalind_perm.txt', 'r') as file:
        n = int(file.read().strip())
        
    # Input validation
    if  n > 7:
        raise ValueError("Error: n must be less than or equal to 7.")
    
    total_permutations, perm_list = permutations(n)
    
    # Output the result
    with open('./output/rosalind_perm.txt', 'w') as file:
        print(total_permutations)
        file.write(f"{total_permutations}\n")
        for p in perm_list:
            perm = ' '.join(map(str, p))
            print(perm)
            file.write(f"{perm}\n")