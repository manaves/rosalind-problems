"""
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, 
    which followed the recurrence relation Fn=Fn−1+Fn−2 and assumed that each pair of 
    rabbits reaches maturity in one month and produces a single pair of offspring 
    (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming 
    solution in the case that all rabbits die out after a fixed number of months. 

Given: Positive integers n≤100 and m≤20.

Return: The total number of pairs of rabbits that will remain after the n-th month if 
    all rabbits live for m months.
"""
MEM = {}
def fibonacci_mortality(n, m):
    """
    Calculate the number of rabbit pairs after n months with mortality of m months.
    
    Parameters:
        n (int): the number of months.
        m (int): the lifespan of each pair of rabbits in months.
        
    Returns:
        int: the total number of rabbit pairs after n months.
    """
    if n in MEM:
        return MEM[n]
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        # New rabbit pairs: previous month + births
        result = fibonacci_mortality(n - 1, m) + fibonacci_mortality(n - 2, m)
        
        # Subtract rabbits that die this month (born m months ago)
        if n - m >= 1:
            result -= fibonacci_mortality(n - m - 1, m)
    
    MEM[n] = result
    return result
    
if __name__ == "__main__":
    n = 0
    m = 0
    with open('./data/rosalind_fibd.txt', 'r') as file:
        n, m = map(int, file.read().strip().split())
    
    # Validate inputs
    if n > 100 or m > 20:
        raise ValueError("Invalid input: n must be less or equal to 100 and m must be less of equal to 20.")
    
    result = fibonacci_mortality(n, m)
    print(result)   