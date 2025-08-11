"""
The task is to use breadth-first search to compute single-source shortest distances in an unweighted directed graph.

Given: A simple directed graph with n≤10^3 vertices in the edge list format.

Return: An array D[1..n] where D[i] is the length of a shortest path from the vertex 1 to the vertex i (D[1]=0). 
    If i is not reachable from 1 set D[i] to −1.
    
Example:
    Input:
    6 6 Number of vertices and edges
    4 6
    6 5
    4 3
    3 5
    2 1
    1 4
    
    Output:
    0 -1 2 1 3 2
"""

def bfs_shortest_path(edges, n):
    """
    Compute single-source shortest distances in an unweighted directed graph using BFS.
    
    Parameters:
        edges (list of tuples): a list of edges where each edge is represented as a tuple
            (u, v) indicating a directed edge from vertex u to vertex v.
        n (int): the number of vertices in the graph.
        
    Returns:
        D (list): a list where the i-th element is the length of the shortest path from vertex 1 to vertex i.
    """
    # Initialize a directed graph as an adjacency list
    graph = {i: [] for i in range(1, n + 1)}
    
    # Build the directed graph from the edge list
    for u, v in edges:
        graph[u].append(v)
        
    # Initialize distances with -1, except for the starting vertex
    D = [-1] * n
    D[0] = 0
    # Initialize a queue for BFS starting from vertex 1 (index 0)
    queue = [0]
    current = 0
    
    # Perform BFS to compute shortest paths
    # Iterate through the queue until all reachable vertices are processed
    while current < len(queue):
        # Get the current vertex from the queue
        u = queue[current]+1
        # Iterate through all neighbors of the current vertex
        for v in graph[u]:
            # If the neighbor has not been visited (distance is -1)
            if D[v-1] == -1:
                # Set the distance to the neighbor and add it to the queue
                D[v-1] = D[u-1] + 1
                queue.append(v-1)
                
        current += 1
    return D

if __name__ == "__main__":
    n, m, edges = None, None, []
    
    with open('./Algorithms/data/rosalind_bfs.txt', 'r') as file:
        lines = file.readlines()
        
    n, m = map(int, lines[0].strip().split())
    edges = [tuple(map(int, line.strip().split())) for line in lines[1:m + 1]]
    
    # Input validation
    if n < 1 or n > 1000:
        raise ValueError("n must be a positive integer between 1 and 1000.")
    
    result = bfs_shortest_path(edges, n)
    result_str = ' '.join(map(str, result))
    
    with open('./Algorithms/output/output_bfs.txt', 'w') as output_file:
        output_file.write(result_str)
        print(result_str)
    