"""
Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population 
    possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the 
    following genotypes:

    AA-AA (100% dominant phenotype)
    AA-Aa (100% dominant phenotype)
    AA-aa (100% dominant phenotype)
    Aa-Aa (75% dominant phenotype)
    Aa-aa (50% dominant phenotype)
    aa-aa (0% dominant phenotype)
Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple 
    has exactly two offspring.
"""
# Coefficients for the expected number of offspring for each genotype pairing
COEFFICIENTS = [1, 1, 1, 0.75, 0.5, 0]

def expected_offspring(integers, of = 2):
    """
    Calculate the expected number of offspring displaying the dominant phenotype.
    
    Parameters:
        integers (list of int): list of six integers representing the number of couples for each genotype
        of (int): number of offspring per couple (default is 2)
        
    Returns:
        float: expected number of offspring displaying the dominant phenotype
    """

    result = sum([of*integers[i]*COEFFICIENTS[i] for i in range(len(integers))])
    return result

if __name__ == "__main__":
    integers = []
    
    with open('./data/rosalind_iev.txt', 'r') as file:
        integers = list(map(int, file.read().split(' ')))
    
    # Input validation
    if len(integers) != 6:
        raise ValueError("The number of integers must be 6.")
    for i in integers:
        if i > 20000:
            raise ValueError("Integers must be less or equal to 20000.")
        
    offspring = expected_offspring(integers)
    print(offspring)
    # Save offspring to file
    with open('./output/output_iev.txt', 'w') as file:
        file.write(str(offspring))