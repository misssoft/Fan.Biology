# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {}
    t = len(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1) #pseudo=1
    # insert your code here
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
            
#Input:  A set of kmers Motifs
#Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    pseudos = CountWithPseudocounts (Motifs)
    sumofsymbol = 0
    for symbol in "ACGT":
        sumofsymbol += pseudos[symbol][0]
    for symbol in "ACGT":
        counts = pseudos[symbol]
        newcounts = []
        for j in counts:
            jprofile = j / sumofsymbol
            newcounts.append(jprofile)
        pseudos[symbol] = newcounts
    profile = pseudos
    return profile

Count = CountWithPseudocounts
Profile = ProfileWithPseudocounts

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
    for i in range(length - k):
        motif = Text[i:i+k]
        probability = Pr(motif, Profile)
        #print("motif:" + motif + "prob: " + str(probability) )
        if probability > count:
            count = probability
            result = motif
        #print("result:" + result )
    return result

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
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

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs



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

def GreedyMotifSearchWithPseudocounts_Workaround(Dna, k, t):
    return GreedyMotifSearch(Dna,k,t)


if __name__ == '__main__':
    #motifs = ['AACGTA','CCCGTT','CACCTT','GGATTA','TTCCGG']
    # print(Count(motifs))
    #print(CountWithPseudocounts(motifs))
    # print(Profile(motifs))
    #testmotifs = ['GTACAACTGT','CAACTATGAA','TCCTACAGGA','AAGCAAGGGT','GCGTACGACC','TCGTCAGCGT','AACAAGGTCA','CTCAGGCGTC','GGATCCAGGT','GGCAAGTACC']
    #print(ProfileWithPseudocounts(testmotifs))
    #Dna = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
    #print(GreedyMotifSearchWithPseudocounts(Dna,3,5))

    DnaTest = [
    'GCACATCATTAAACGATTCGCCGCATTGCCTCGATAGGCG',
    'TCATAACTGACACCTGCTCTGGCACCGCTCATCCGTCGAA',
    'AAGCGGGTATAGCCAGATAGTGCCAATAATTTCCTTCGGC',
    'AGTCGGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
    'AACCGGACGGCAACTACGGTTACAACGCAGCAAGAATATT',
    'AGGCGTCTGTTGTTGCTAACACCGTTAAGCGACGGCAACT',
    'AAGCTTCCAACATCGTCTTGGCATCTCGGTGTGTGAGGCG',
    'AATTGAACATCTTACTCTTTTCGCTTTCAAAAAAAAGGCG']
    print(GreedyMotifSearchWithPseudocounts(DnaTest,5,8))
