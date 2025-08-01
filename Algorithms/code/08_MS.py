"""
The problem of sorting a list of numbers lends itself immediately to a divide-and-conquer strategy: 
    split the list into two halves, recursively sort each half, and then merge the two sorted sublists 
    (recall the problem “Merge Two Sorted Arrays”).

Given: A positive integer n≤10^5 and an array A[1..n] of integers from −10^5 to 10^5.

Return: A sorted array A[1..n].

Example:
    Input:
    10
    20 19 35 -18 17 -20 20 1 4 4
    
    Output: 
    -20 -18 1 4 4 17 19 20 20 35
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

def merge_sort(arr, n):
    """
    Sort an array using the merge sort algorithm.
    
    Parameters:
        arr (list): a list of integers to be sorted.
        n (int): the number of elements in the array.
        
    Returns:
        list: a sorted list containing all elements from arr.
    """
    # Base case: if the array has one or no elements, it is already sorted
    if n <= 1:
        return arr
    # Split the array into two halves and recursively sort each half
    mid = n // 2
    left_half = merge_sort(arr[:mid], len(arr[:mid]))
    right_half = merge_sort(arr[mid:], len(arr[mid:]))
    
    return merge_sorted_arrays(left_half, len(left_half), right_half, len(right_half))

if __name__ == "__main__":
    arr, n = None, None
    
    with open('./Algorithms/data/rosalind_ms.txt', 'r') as file:
        lines = file.read().splitlines()
    
    n = int(lines[0])
    arr = list(map(int, lines[1].split()))
    
    # Input validation
    if n < 1 and n > 100000:
        raise ValueError("n must be a positive integer between 1 and 100,000.")
    if any(num < -100000 or num > 100000 for num in arr):
        raise ValueError("Each number in the array must be between -100,000 and 100,000.")
    if len(arr) != n:
        raise ValueError("The length of the array must be equal to n.")
    
    sorted_arr = merge_sort(arr, n)
    str_sorted_arr = ' '.join(map(str, sorted_arr))
    
    with open('./Algorithms/output/output_ms.txt', 'w') as output_file:
        output_file.write(str_sorted_arr)
        print(str_sorted_arr)