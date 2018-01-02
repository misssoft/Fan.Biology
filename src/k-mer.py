# Input:  Two strings, Pattern and Genome
# Output: A list containing all starting positions where Pattern appears as a substring of Genome
def PatternMatching(pattern, genome):
    positions = [] # output variable
    # your code here
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
    return positions

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(pattern, genome, d):
    positions = [] # output variable
    for i in range(len(genome)-len(pattern)+1):
        variants = HammingDistance(pattern,genome[i:i+len(pattern)])        
        if variants <= d:
            positions.append(i)
    return positions

'''
Input:  Strings Pattern and Text
Output: The number of times Pattern appears in Text
PatternCount("ACTAT", "ACAACTATGCATACTATCGGGAACTATCCT") = 3.
'''
def PatternCount(Pattern, Text):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        variants = HammingDistance(Pattern,Text[i:i+len(Pattern)])        
        if variants <= d:
            count = count + 1
    return count

'''
Input:  Strings Text and integer k
Output: The dictionary of the count of every k-mer string  
CountDict("CGATATATCCATAG",3) = {0: 1, 1: 1, 2: 3, 3: 2, 4: 3, 5: 2, 6: 1, 7: 1, 8: 1, 9: 1, 10: 3, 11: 1}.
'''
def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

'''
Input:  a list of items
Output: the non-duplicate items 
Remove_duplicates([1,2,3,3]) = [1,2,3]
'''
def Remove_duplicates(Items):
    ItemsNoDuplicates = [] # output variable
    ItemsNoDuplicates = list(set(Items))
    return ItemsNoDuplicates

'''
Input:  Strings Text and integer k
Output: The freqent k-mer 
FrequentWords("GATCCAGATCCCCATAC",2) = "CC".
'''
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m= max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = Remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates 

# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

# Input:  A String Genome
# Output: Skew(Genome)
def Skew(Genome):
    n = len(Genome)
    skew = {} #initializing the dictionary
    skew[0] = 0
    for i in range(1,n+1):
        skew[i] = skew[i-1]
        if Genome[i-1] == 'G':
            skew[i] = skew[i] + 1
        if Genome[i-1] == 'C':
            skew[i] = skew[i] - 1         
    return skew

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = [] # output variable
    # your code here
    a_skew = Skew(Genome)
    values = list(a_skew.values())
    min_value = min(values)
    for i in range(len(a_skew)):
        if (a_skew[i] == min_value):
            positions.append(i)
    return positions

# Input:  Two strings p and q (same length)
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    m = len(p)
    n = len(q)
    result = 0
    if (m != n):
        return result
    else:
        for i in range(n):
            if (p[i] != q[i]):
                result = result + 1
    return result



### DO NOT MODIFY THE CODE BELOW THIS LINE ###
if __name__ == '__main__':
    target ='TGATCA'
    gene =str('ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC')
    #print(PatternCount(target, gene))
    #print (CountDict(gene,10))
    #print (FrequentWords(gene,10))
    #print(PatternCount("TGT", "ACTGTACGATGATGTGTGTCAAAG"))
    #print(PatternMatching('ATAT','GATATATGCATATACTT'))
    print(SymbolArray('AAAAGGGG','A'))
    print(FasterSymbolArray('AAAAGGGG','A'))
    print(Skew('CATGGGCATCGGCCATACGCC'))
    print(MinimumSkew('CATTCCAGTACTTCGATGATGGCGTGAAGA'))
    print(MinimumSkew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'))
    print(HammingDistance('GGGCCGTTGGT','GGACCGTTGAC'))
    print(ApproximatePatternMatching('ATTCTGGA','CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',3))
    print(ApproximatePatternCount('GAGG','TTTAGAGCCTTCAGAGG',2))
    print(HammingDistance('CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG','ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'))