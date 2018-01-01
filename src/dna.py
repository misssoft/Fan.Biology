# Input:  Two strings, Pattern and Genome
# Output: A list containing all starting positions where Pattern appears as a substring of Genome
def PatternMatching(pattern, genome):
    positions = [] # output variable
    # your code here
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions


# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(genomo):
    rev_genemo = genomo[::-1]
    rev_comp = '' # output variable
    for i in rev_genemo:
        comp = Complement(i)
        rev_comp = rev_comp + comp
    return rev_comp


# Copy your reverse function from the previous step here.


# HINT:   Filling in the following function is optional, but it may come in handy when solving ReverseComplement
# Input:  A character Nucleotide
# Output: The complement of Nucleotide
def Complement(nucleotide):
    comp = '' # output variable
    # your code here
    if nucleotide == 'A':
        comp= 'T'
    if nucleotide == 'T':
        comp=  'A'
    if nucleotide == 'C':
        comp=  'G'
    if nucleotide == 'G':
        comp=  'C'
    return comp

if __name__ == '__main__':
    print(Complement('A'))
    print(Complement('G'))
    print(Complement('C'))
    print(Complement('T'))
    print(ReverseComplement('AAAACCCGGT'))
    print(PatternMatching('ATAT','GATATATGCATATACTT'))

    print(ReverseComplement('TTGTGTC'))