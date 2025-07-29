"""
BINARY SEARCH

Given: Two positive integers n≤10^5 and m≤10^5 , a sorted array A[1..n] of integers from −10^5 
    to 10^5 and a list of m  integers −10^5≤k1,k2,…,km≤10^5 
 
Return: For each ki, output an index 1≤j≤n  s.t. A[j]=ki  or "-1" if there is no such index.
"""

def binary_search(arr, x):
    """
    Perform binary search on a sorted array to find the index of a given element.
    
    Parameters:
        arr (list): A sorted list of integers.
        x (int): The integer to search for in the array.
    
    Returns:
        int: The index of x in arr if found, otherwise -1.
    """
    left, right = 0, len(arr) - 1
    # While loop to perform binary search
    while left <= right:
        # Calculate the middle index
        mid = (left + right) // 2
        
        if arr[mid] == x: # Check if the middle element is the target
            return mid + 1
        elif arr[mid] < x: # If the middle element is less than x, search in the right half
            left = mid + 1
        else: # If the middle element is greater than x, search in the left half
            right = mid - 1
    # If the element is not found, return -1
    return -1
        

if __name__ == "__main__":
    n, m = None, None
    arr = []
    ks = []
    
    with open('./Algorithms/data/rosalind_bins.txt', 'r') as file:
        lines = file.read().splitlines()

    n = int(lines[0])
    m = int(lines[1])
    arr = list(map(int, lines[2].split()))
    ks = list(map(int, lines[3].split()))
            
    # Input validation
    if n > 100000 or m > 100000:
        raise ValueError("n and m must be less than or equal to 100000.")
    if len(arr) != n or len(ks) != m:
        raise ValueError("The length of arr must be equal to n and the length of ks must be equal to m.")
    if any(x < -100000 or x > 100000 for x in arr + ks):
        raise ValueError("All elements in arr and ks must be between -100000 and 100000.")
    
    results = [binary_search(arr, k) for k in ks]
    
    with open('./Algorithms/output/output_bins.txt', 'w') as output_file:
        output_file.write(' '.join(map(str, results)))
        print(' '.join(map(str, results)))