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


### DO NOT MODIFY THE CODE BELOW THIS LINE ###
if __name__ == '__main__':
    target ='TGATCA'
    gene =str('ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC')
    print(PatternCount(target, gene))
    print (CountDict(gene,10))
    print (FrequentWords(gene,10))