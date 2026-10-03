# TODO:
# Implement BLAST integration.
#
# Questions to reason through:
# - What sequence are we sending to BLAST?
# - Which BLAST program is appropriate?
# - Which database are we searching?
# - What does percent identity mean?
# - What does query coverage mean?
# - What does E-value mean?
# - How should BLAST results be parsed?
def BLAST_search(sequence, program='blastn', database='nt'):
    """
    Perform a BLAST search for the given sequence.

    Parameters:
    - sequence: The nucleotide or protein sequence to search.
    - program: The BLAST program to use (e.g., 'blastn', 'blastp').
    - database: The database to search against (e.g., 'nt', 'nr').

    Returns:
    - A list of BLAST results with relevant information.
    """
    