"""
Given: A genus name, followed by two dates in YYYY/M/D format.

Return: The number of Nucleotide GenBank entries for the given genus that were published between the dates specified.
"""

from Bio import Entrez 

def find_entries(name, start, end):
    """
    Finds the number of Nucleotide GenBank entries for a given genus name that were published between two dates.
    
    Parameters:
        name (str): the genus name to search for.
        start (str): the start date in YYYY/M/D format.
        end (str): the end date in YYYY/M/D format.
    
    Returns:
        int: the number of entries found.
    """
    # Construct the search term with the genus name and date range
    search_term = f'{name}[Organism] AND ("{start}"[Publication Date] : "{end}"[Publication Date])'
    handle = Entrez.esearch(db='nucleotide', term=search_term)
    record = Entrez.read(handle)
    handle.close()
    # Return the count of entries found
    return record['Count']

if __name__ == "__main__":
    name, start, end = None, None, None
    
    with open('./Armory/data/rosalind_gbk.txt', 'r') as file:
        name, start, end = map(str, file.read().strip().split('\n'))
    
    # Input validation
    # Check start and end dates format based in the pattern YYYY/M/D
    if not (start and end):
        raise ValueError("Start and end dates must be provided.")
    if not (start.count('/') == 2 and end.count('/') == 2):
        raise ValueError("Start and end dates must be in the format YYYY/M/D.")
    try:
        start_year, start_month, start_day = map(int, start.split('/'))
        end_year, end_month, end_day = map(int, end.split('/'))
    except ValueError:
        raise ValueError("Start and end dates must be valid integers in the format YYYY/M/D.")
    
    if not (1 <= start_month <= 12 and 1 <= end_month <= 12):
        raise ValueError("Months must be between 1 and 12.")
    if not (1 <= start_day <= 31 and 1 <= end_day <= 31):
        raise ValueError("Days must be between 1 and 31.")
    if start_year > end_year or (start_year == end_year and (start_month > end_month or (start_month == end_month and start_day > end_day))):
        raise ValueError("Start date must be earlier than or equal to end date.")
    
    # Check name format
    if not name.isalpha():
        raise ValueError("Genus name must contain only alphabetic characters.")
    
    # Call the function to find entries
    Entrez.email = "manaves@uoc.edu"
    entries_count = find_entries(name, start, end)
    
    with open('./Armory/output/output_gbk.txt', 'w') as output_file:
        print(entries_count)
        output_file.write(str(entries_count))