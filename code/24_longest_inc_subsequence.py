"""
A subsequence of a permutation is a collection of elements of the permutation in the order that they appear. 
    For example, (5, 3, 4) is a subsequence of (5, 1, 3, 4, 2).

A subsequence is increasing if the elements of the subsequence increase, and decreasing if the elements decrease. 
    For example, given the permutation (8, 2, 1, 6, 5, 7, 4, 3, 9), an increasing subsequence is (2, 6, 7, 9), 
    and a decreasing subsequence is (8, 6, 5, 4, 3). You may verify that these two subsequences are as long as possible.

Given: A positive integer n≤10000 followed by a permutation π of length n.

Return: A longest increasing subsequence of π, followed by a longest decreasing subsequence of π.
"""

def longest_subsequence(n, permutation, inc):
    # dp stores the length of the longest subsequence
    dp = [1]*n
    # prev stores the previous index in the subsequence
    prev = [-1]*n
    
    for i in range(n):
        for j in range(i):
            # Check if seq[j] can come before seq[i] in the subsequence
            if (inc == True and permutation[j] < permutation[i]) or (inc == False and permutation[j] > permutation[i]):
                # If connecting j to i forms a longer subsequence, update dp[i] and prev[i]
                if dp[j]+1 > dp[i]:
                    dp[i] = dp[j]+1
                    prev[i] = j
                    
    # Find the index where the longest subsequence ends
    max_index = dp.index(max(dp))
    
    # Reconstruct the subsequence by following prev backwards
    result = []
    while max_index != -1:
        result.append(permutation[max_index])
        max_index = prev[max_index]

    # Reverse the result to get the correct order
    return result[::-1]

if __name__ == "__main__":
    data = []
    
    with open('./data/rosalind_lgis.txt', 'r') as file:
        data = file.read().strip().split('\n')
        
    n = int(data[0])
    permutation = tuple(map(int, data[1].split(' ')))
    
    # Input validation
    if n > 10000:
        raise ValueError("n must be less than or equal to 10000.")
    
    if n != len(permutation):
        raise ValueError("n must be the length of permutation.")
    
    subseq_inc = ' '.join(map(str, longest_subsequence(n, permutation, inc=True)))
    subseq_dec = ' '.join(map(str, longest_subsequence(n, permutation, inc=False)))
    
    with open('./output/output_lgis.txt', 'w') as output_file:
        
        print(subseq_inc)
        output_file.write(subseq_inc+'\n')
        print(subseq_dec)
        output_file.write(subseq_dec)
    