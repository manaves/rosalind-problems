"""
The task is to use depth-first search to compute the number of connected components in a given undirected graph.

Given: A simple graph with n≤10^3 vertices in the edge list format.

Return: The number of connected components in the graph.

Example:
    Input:
    12 13
    1 2
    1 5
    5 9
    5 10
    9 10
    3 4
    3 7
    3 8
    4 8
    7 11
    8 11
    11 12
    8 12
    
    Output:
    3
"""

def count_cc(n, edges):
    """
    Count the number of connected components in an undirected graph.
    
    Parameters:
        n (int): the number of vertices in the graph.
        edges (list): a list of tuples representing the edges in the graph.
        
    Returns:
        count (int): the number of connected components in the graph.
    """

    def dfs(vertex):
        """
        Perform depth-first search to visit all vertices in the connected component.
        
        Parameters:
            vertex (int): the current vertex to start DFS from.
            
        Returns:
            None: modifies the visited list in place.
        """
        # Initialize a stack for DFS
        stack = [vertex]
        # While there are vertices to visit
        while stack:
            # Pop a vertex from the stack
            node = stack.pop()
            # If the node not exists in the graph, continue to the next iteration
            if node not in graph.keys():
                continue
            # Visit all unvisited neighbors of the current node
            for neighbor in graph[node]:
                # If the neighbor has not been visited, mark it as visited and add it to the stack
                if not visited[neighbor]:
                    visited[neighbor] = True
                    stack.append(neighbor)
    
    # Initialize the graph as an adjacency list
    # Each vertex will map to a list of its neighbors
    graph = {}
    # Initialize visited array to keep track of visited vertices
    # 1-based indexing, so we use n + 1 to accommodate vertex indices from 1 to n
    visited = [False] * (n + 1)
        
    # Build the graph from edges
    for u, v in edges:
        if u not in graph:
            graph[u] = []
        if v not in graph:
            graph[v] = []
        # Append edges in both directions since the graph is undirected
        graph[u].append(v)
        graph[v].append(u)

    count = 0
    for vertex in range(1, n + 1):
        # If the vertex has not been visited, it starts a new connected component
        if not visited[vertex]:
            visited[vertex] = True
            dfs(vertex)
            count += 1

    return count

if __name__ == "__main__":
    n, m = None, None
    edges = []
    
    with open('./Algorithms/data/rosalind_cc.txt', 'r') as file:
        lines = file.read().splitlines()
        
    n, m = map(int, lines[0].split())
    for line in lines[1:m + 1]:
        u, v = map(int, line.split())
        edges.append((u, v))
    
    # Input validation
    if n < 1 or n > 1000:
        raise ValueError("n must be a positive integer between 1 and 1000.")
    if m < 0 or m > (n * (n - 1)) // 2:
        raise ValueError("m must be a non-negative integer not exceeding the maximum number of edges in a simple graph.")

    result = count_cc(n, edges)
    
    with open('./Algorithms/output/output_cc.txt', 'w') as output_file:
        output_file.write(str(result))
        print(result)