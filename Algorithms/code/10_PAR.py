"""
Given: A positive integer n≤10^5 and an array A[1..n] of integers from −10^5 to 10^5.

Return: A permuted array B[1..n] such that it is a permutation of A and there is an 
    index 1≤q≤n such that B[i]≤A[1] for all 1≤i≤q−1, B[q]=A[1], and B[i]>A[1] for all q+1≤i≤n.

Example:
    Input:
    9
    7 2 5 6 1 3 9 4 8
    
    Output: 
    5 6 3 4 1 2 7 9 8
"""

def permute_array(A, n):
    """
    Permute an array such that the first element is at a specific index
    and all elements before it are less than or equal to it, and all elements after it
    are greater than it.
    
    Parameters:
        A (list): a list of integers to be permuted.
        n (int): the number of elements in the array.
        
    Returns:
        B (list): a permuted list of integers.
    """
    
    # Find the first element to be the pivot
    pivot = A[0]
    less_than_or_equal, greater_than = [], []
    
    # Partition the array into two parts
    for i in range(1, n):
        if A[i] <= pivot:
            less_than_or_equal.append(A[i])
        else:
            greater_than.append(A[i])

    # Combine the parts with the pivot in the middle    
    B = less_than_or_equal + [pivot] + greater_than
            
    return B

if __name__ == "__main__":
    n, A = None, []
    
    with open('./Algorithms/data/rosalind_par.txt', 'r') as file:
        lines = file.read().splitlines()
        
    n = int(lines[0])
    A = list(map(int, lines[1].split()))
    
    # Input validation
    if n < 1 or n > 100000:
        raise ValueError("n must be a positive integer between 1 and 100,000.")
    if any(num < -100000 or num > 100000 for num in A):
        raise ValueError("Each number in the array must be between -100,000 and 100,000.")
    if len(A) != n:
        raise ValueError("The length of the array must be equal to n.")
    
    # Get the permuted array
    permuted_array = permute_array(A, n)
    str_permuted_array = ' '.join(map(str, permuted_array))
    
    with open('./Algorithms/output/output_par.txt', 'w') as output_file:
        output_file.write(str_permuted_array)
        print(str_permuted_array)