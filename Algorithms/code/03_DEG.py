"""
In an undirected graph, the degree d(u) of a vertex u is the number of neighbors u has,
    or equivalently, the number of edges incident upon it.

Given: A simple graph with n≤10^3 vertices in the edge list format.

Return: An array D[1..n] where D[i] is the degree of vertex i.

Example:

Input:
    6 7 # n and m, number of vertices and edges
    1 2
    2 3
    6 3
    5 6
    2 5
    2 4
    4 1

Return:
    2 4 2 2 2 2
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

if __name__ == "__main__":
    n, m = None, None
    edges = []
    
    with open('./Algorithms/data/rosalind_deg.txt', 'r') as file:
        lines = file.read().splitlines()
    
    n, m = map(int, lines[0].split())
    edges = [tuple(map(int, line.split())) for line in lines[1:m+1]]
    
    # Input validation
    if n > 1000:
        raise ValueError("n must be less than or equal to 1000.")
    if len(edges) != m:
        raise ValueError("The number of edges must be equal to m.")
    
    degrees = degree_of_vertices(edges, n)
    
    with open('./Algorithms/output/output_deg.txt', 'w') as output_file:
        output_file.write(' '.join(map(str, degrees)))
        print(' '.join(map(str, degrees)))