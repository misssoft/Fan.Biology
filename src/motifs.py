# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts_Worked(Motifs):
    count = Count(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    for symbol in "ACGT":
        counts = count[symbol]
        newcounts = []
        for j in counts:
            jplus = j + 1
            newcounts.append(jplus)
        count[symbol] = newcounts
    return count
            
    # insert your code here

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += (1/t)
    return count

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    consensus = Consensus(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    score = []
    for i in range(k):
        count = 0
        for j in range(t):
            if Motifs[j][i] != consensus[i]:
                count += 1
        score.append(count)
    return sum(score)

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    length = len(Text)
    result = 1
    for i in range(length):
        symbol = Text[i]
        probability = Profile[symbol][i]
        result = result * probability
    return result 

# Input:  String Text, an integer k, and profile matrix Profile
# Output: ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    length = len(Text)
    count = -1
    result = Text[0:k]
    for i in range(length - k + 1):
        motif = Text[i:i+k]
        probability = Pr(motif, Profile)
        #print("motif:" + motif + "prob: " + str(probability) )
        if probability > count:
            count = probability
            result = motif
        #print("result:" + result )
    return result

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


if __name__ == '__main__':
    DnaTest = ['GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT',
               'CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT',
               'GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT',
               'AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC',
               'ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG']
    print(GreedyMotifSearch(DnaTest,6,5))



