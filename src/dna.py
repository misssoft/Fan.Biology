# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(genomo):
    rev_genemo = genomo[::-1]
    rev_comp = '' # output variable
    for i in rev_genemo:
        comp = Complement(i)
        rev_comp = rev_comp + comp
    return rev_comp

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
    print(ReverseComplement('TTGTGTC'))