"""
Given: A simple graph with n≤10^3 vertices in the edge list format.

Return: An array D[1..n] where D[i] is the sum of the degrees of i's neighbors.

Example:
    Input:
    5 4 # n and m, number of vertices and edges
    1 2
    2 3
    4 3
    2 4
    
    Output:
    3 5 5 5 0
"""

def degree_of_vertices(edges, n):
    """
    Calculate the degree of each vertex in an undirected graph.
    
    Parameters:
        edges (list of tuples): a list of edges where each edge is represented as a tuple
            (u, v) indicating an edge between vertices u and v.
        n (int): the number of vertices in the graph.
        
    Returns:
        list: a list where the i-th element is the degree of vertex i+1.
    """
    # Initialize a list to store the degree of each vertex
    degrees = [0] * n
    
    # Iterate through each edge and increment the degree for both vertices
    for edge in edges:
        u, v = edge
        # Increment the degree for both vertices u and v
        degrees[u - 1] += 1
        degrees[v - 1] += 1
        
    return degrees

def neighbor_degrees_sum(edges, n):
    """
    Calculate the sum of degrees of neighbors for each vertex in an undirected graph.
    
    Parameters:
        edges (list of tuples): a list of edges where each edge is represented as a tuple
            (u, v) indicating an edge between vertices u and v.
        n (int): the number of vertices in the graph.
        
    Returns:
        neighbor_sums (list): a list where the i-th element is the sum of degrees of neighbors of vertex i+1.
    """
    
    # Calculate the degree of each vertex
    degrees = degree_of_vertices(edges, n)
    
    # Initialize a list to store the sum of degrees of neighbors for each vertex
    neighbor_sums = [0] * n
    
    # Iterate through each edge and add the degree of the neighbor to the sum
    for edge in edges:
        u, v = edge
        # Add the degree of v to the sum for u and vice versa
        neighbor_sums[u - 1] += degrees[v - 1]
        neighbor_sums[v - 1] += degrees[u - 1]
        
    return neighbor_sums

if __name__ == "__main__":
    
    n, m, edges = None, None, []
    
    with open('./Algorithms/data/rosalind_ddeg.txt', 'r') as file:
        lines = file.read().splitlines()
        
    n, m = map(int, lines[0].split())
    edges = [tuple(map(int, line.split())) for line in lines[1:m+1]]
    
    # Input validation
    if n < 1 or n > 1000:
        raise ValueError("n must be a positive integer between 1 and 1000.")
    if len(edges) != m:
        raise ValueError("The number of edges must be equal to m.")
    
    neighbor_sums = neighbor_degrees_sum(edges, n)
    str_results = ' '.join(map(str, neighbor_sums))
    
    with open('./Algorithms/output/output_ddeg.txt', 'w') as output_file:
        output_file.write(str_results)
        print(str_results)
