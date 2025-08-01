"""
The merging procedure is an essential part of “Merge Sort” (which is considered in one of the next problems).

Given: A positive integer n≤10^5 and a sorted array A[1..n] of integers from −10^5 to 10^5, a positive integer m≤10^5
    and a sorted array B[1..m] of integers from −10^5 to 10^5.

Return: A sorted array C[1..n+m] containing all the elements of A and B.

Example:
    Input:
    4
    2 4 10 18
    3
    -5 11 12
    
    Output:
    -5 2 4 10 11 12 18
"""

def merge_sorted_arrays(A, n, B, m):
    """
    Merge two sorted arrays into one sorted array.
    
    Parameters:
        A (list): a sorted list of integers.
        n (int): the number of elements in A.
        B (list): another sorted list of integers.
        m (int): the number of elements in B.
        
    Returns:
        C (list): a merged and sorted list containing all elements from A and B.
    """
    # Initialize an empty list to store the merged result
    C = []
    
    i, j = 0, 0
    
    # Merge the two arrays
    while i < n and j < m:
        if A[i] <= B[j]:
            C.append(A[i])
            i += 1
        else:
            C.append(B[j])
            j += 1
            
    # If there are remaining elements in A
    while i < n:
        C.append(A[i])
        i += 1
        
    # If there are remaining elements in B
    while j < m:
        C.append(B[j])
        j += 1
        
    return C


if __name__ == "__main__":
    n, A, m, B = None, [], None, []
    
    with open('./Algorithms/data/rosalind_mer.txt', 'r') as file:
        lines = file.read().splitlines()
    
    n = int(lines[0])
    A = list(map(int, lines[1].split()))
    m = int(lines[2])
    B = list(map(int, lines[3].split()))
    
    # Input validation
    if n < 1 or n > 100000:
        raise ValueError("n must be a positive integer between 1 and 100,000.")
    if m < 1 or m > 100000:
        raise ValueError("m must be a positive integer between 1 and 100,000.")
    if len(A) != n:
        raise ValueError("The length of array A must be equal to n.")
    if len(B) != m:
        raise ValueError("The length of array B must be equal to m.")
    if not all(-100000 <= x <= 100000 for x in A + B):
        raise ValueError("All elements in arrays A and B must be between -100,000 and 100,000.")
    
    result = merge_sorted_arrays(A, n, B, m)
    
    with open('./Algorithms/output/output_mer.txt', 'w') as output_file:
        str_results = ' '.join(map(str, result))
        print(str_results)
        output_file.write(str_results)