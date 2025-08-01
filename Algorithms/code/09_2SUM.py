"""
Given: A positive integer k≤20, a positive integer n≤10^4, and k arrays of size n containing integers from −10^5 to 10^5.

Return: For each array A[1..n], output two different indices 1≤p<q≤n such that A[p]=−A[q] if exist, and "-1" otherwise.

Example:
    Input:
    4 5
    2 -3 4 10 5
    8 2 4 -2 -8
    -5 2 3 2 -4
    5 4 -5 6 8
    
    Output:
    -1
    2 4
    -1
    1 3
"""

def find_two_sum_pairs(n, arrays):
    """
    Find pairs of indices in each array such that the elements at those indices sum to zero.
    
    Parameters:
        n (int): the size of each array.
        arrays (list of lists): a list containing k arrays, each of size n.
        
    Returns:
        results (list): a list where each element is a string representing the indices of the pair
            that sums to zero, or "-1" if no such pair exists.
    """
    results = []
    
    for arr in arrays:
        # Initialize a dictionary to keep track of seen numbers and their indices
        num_seen = {}
        found = False
        # Iterate through the array to find pairs
        for i in range(n):
            # Check if the negative of the current number exists in the seen numbers
            num_check = -arr[i]
            # If it exists, we have found a pair
            if num_check in num_seen:
                results.append(f"{num_seen[num_check] + 1} {i + 1}")
                found = True
                break
            # Store the current number and its index
            num_seen[arr[i]] = i
            
        if not found:
            results.append("-1")
            
    return results

if __name__ == "__main__":
    k, n, arrays = None, None, []
    
    with open('./Algorithms/data/rosalind_2sum.txt', 'r') as file:
        lines = file.read().splitlines()
        
    k, n = map(int, lines[0].split())
    arrays = [list(map(int, line.split())) for line in lines[1:k+1]]
    
    # Input validation
    if n < 1 or n > 10000:
        raise ValueError("n must be a positive integer between 1 and 10,000.")
    if k < 1 or k > 20:
        raise ValueError("k must be a positive integer between 1 and 20.")
    
    if any(num < -100000 or num > 100000 for arr in arrays for num in arr):
        raise ValueError("Each number in the arrays must be between -100,000 and 100,000.")
    
    results = find_two_sum_pairs(n, arrays)
    
    with open('./Algorithms/output/output_2sum.txt', 'w') as output_file:
        output_file.write('\n'.join(results))
        print('\n'.join(results))