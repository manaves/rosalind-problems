
def fibonacci(n):
    """
    Calculate the nth Fibonacci number using recursion.
    
    Parameters:
        n (int): the position in the Fibonacci sequence (0-indexed).
    
    Returns:
        int: the nth Fibonacci number.
    """
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        return fibonacci(n-1)+fibonacci(n-2)

if __name__ == "__main__":
    n = None
    
    with open('./Algorithms/data/rosalind_fibo.txt', 'r') as file:
        n = int(file.read().strip())
        
    # Input validation
    if n > 25:
        raise ValueError("The value of n must be less than or equal to 25.")
    
    r = fibonacci(n)
    
    with open('./Algorithms/output/output_fibo.txt', 'w') as output_file:
        print(r)
        output_file.write(str(r))