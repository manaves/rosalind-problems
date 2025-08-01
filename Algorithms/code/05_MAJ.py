"""
An array A[1..n] is said to have a majority element if more than half of its entries are the same.

Given: A positive integer k≤20, a positive integer n≤10^4, and k arrays of size n 
    containing positive integers not exceeding 10^5.

Return: For each array, output an element of this array occurring strictly more than 
    n/2 times if such element exists, and "-1" otherwise.
"""

def majority_element(arr_list, n):
    """
    Calculate the majority element for each array in the list.
    
    Parameters:
        arr_list (list of lists): a list containing k arrays, each of size n.
        n (int): the size of each array.
        
    Returns:
        results (list): a list where each element is the majority element of the corresponding array,
            or -1 if no majority element exists.
    """
    
    results = []
    
    for arr in arr_list:
        count = {}
        
        for num in arr:
            if num not in count.keys():
                count[num] = 1
            else:
                count[num] += 1
        
        max_num = max(count.values())
        
        if max_num <= (n//2):
            results.append(-1)
            continue
        else:
            key_max = [key for key in count if count[key] == max_num]
            results.append(key_max[0])
        
    return results
    
if __name__ == "__main__":
    k, n, arr_list = None, None, []
    
    with open('./Algorithms/data/rosalind_maj.txt', 'r') as file:
        lines = file.read().splitlines()
    
    k, n = map(int, lines[0].split())
    
    arr_list = [list(map(int, line.split())) for line in lines[1:k+1]]
    
    # Input validation
    if n < 1 or n > 10000:
        raise ValueError("n must be a positive integer between 1 and 10,000.")
    if k < 1 or k > 20:
        raise ValueError("k must be a positive integer between 1 and 20.")
    if len(arr_list) != k:
        raise ValueError("The number of arrays must be equal to k.")
    if any(len(arr) != n for arr in arr_list):
        raise ValueError("Each array must have a length equal to n.")
    if any(num < 1 or num > 100000 for arr in arr_list for num in arr):
        raise ValueError("Each number in the arrays must be between 1 and 100,000.")
    
    results = majority_element(arr_list, n)
    str_results = ' '.join(map(str, results))
    
    with open('./Algorithms/output/output_maj.txt', 'w') as output_file:
        print(str_results)
        output_file.write(str_results)