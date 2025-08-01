"""
Given: A positive integer n≤10^3 and an array A[1..n] of integers.

Return: The number of swaps performed by insertion sort algorithm on A[1..n].
"""

def insertion_sort(arr, n):
    """
    Perform insertion sort on the array and count the number of swaps.
    
    Parameters:
        arr (list): a list of integers to be sorted.
        n (int): the number of elements in the array.
        
    Returns:
        counts (int): the number of swaps performed during the sorting.
    """
    counts = 0
    for i in range(1, n):
        k = i
        # Move the current element to its correct position in the sorted part of the array
        while k > 0 and arr[k] < arr[k-1]:
            # Swap elements
            arr[k], arr[k-1] = arr[k-1], arr[k]
            k -= 1
            counts += 1
    
    return counts

if __name__ == "__main__":
    arr, n = None, None
    
    with open('./Algorithms/data/rosalind_ins.txt', 'r') as file:
        lines = file.read().splitlines()
    
    n = int(lines[0])
    arr = list(map(int, lines[1].split()))
    
    # Input validation
    if n > 1000:
        raise ValueError("n must be less than or equal to 1000.")
    if n < 1:
        raise ValueError("n must be a positive integer.")
    if len(arr) != n:
        raise ValueError("The length of the array must be equal to n.")
    
    swaps = insertion_sort(arr, n)
    
    with open('./Algorithms/output/output_ins.txt', 'w') as output_file:
        print(swaps)
        output_file.write(str(swaps))