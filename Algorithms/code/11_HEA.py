"""
A binary heap is a binary tree based data structure that is often used to implement priority queues. 
    Binary heaps, in turn, can be easily implemented using an array if the underlying tree is a 
    complete binary tree. The tree nodes have a natural ordering: row by row, starting at the root 
    and moving left to right within each row. If there are n nodes, this ordering specifies their 
    positions 1,2,…,n within the array. Moving up and down the tree is easily simulated on the array, 
    using the fact that node number j has parent ⌈j/2⌉ and children 2j and 2j+1.

The goal of this problem is to build a heap from the given array. For this, go from the end to the 
    beginning of a given array and let each element "bubble up".
    
Given: A positive integer n≤10^5 and array A[1..n] of integers from −10^5 to 10^5.

Return: A permuted array A satisfying the binary max heap property: for any 2≤i≤n, A[⌊i/2⌋]≥A[i].

Example:
    Input:
    5
    1 3 5 7 2
    
    Output:
    7 5 1 3 2
"""

def max_heap(array, n):
    """
    Build a max heap from the given array.
    
    Parameters:
        array (list): a list of integers to be transformed into a max heap.
        n (int): the number of elements in the array.
    
    Returns:
        A (list): a permuted list of integers that satisfies the binary max heap property.
    """
    
    # Heapify the array by bubbling up each element
    def heapify(i):
        """
        Heapify the subtree rooted at index i.
        
        Parameters:
            i (int): the index of the root of the subtree to be heapified.
        
        Returns:
            None: modifies the array in place.
        """
        
        while 2 * i <= n:
            # Find the largest among the root, left child, and right child
            left = 2 * i
            right = 2 * i + 1
            largest = i

            # Check if the left exists and is greater than the largest found so far
            if left <= n and A[left] > A[largest]:
                largest = left
            # Check if the right exists and is greater than the largest found so far
            if right <= n and A[right] > A[largest]:
                largest = right
            # If the largest is not the root, swap and continue heapifying
            if largest != i:
                A[i], A[largest] = A[largest], A[i]
                i = largest  # Continue sifting down
            else:
                break
    
    # Initialize the array with a dummy value at index 0 for easier calculations
    A = [None] + array[:]
    
    # Heapify the array starting from the last non-leaf node
    for i in range(n//2, 0, -1):
        heapify(i)
    
    return A[1:] # Return the array without the dummy value

if __name__ == "__main__":
    n, A = None, []
    
    with open('./Algorithms/data/rosalind_hea.txt', 'r') as file:
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
    
    result = max_heap(A, n)
    
    with open('./Algorithms/output/output_hea.txt', 'w') as output_file:
        print(' '.join(map(str, result)))
        output_file.write(' '.join(map(str, result)))