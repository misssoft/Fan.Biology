import random

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {}
    t = len(Motifs)
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)  # pseudo=1
    # insert your code here
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    pseudos = CountWithPseudocounts(Motifs)
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


# Insert your Pr(Text, Profile) function here from Motifs.py.
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
        if probability > count:
            count = probability
            result = motif
    return result

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    k = len(Profile['A'])
    t = len(Dna)
    result = []
    for i in range(t):
        probablekmer = ProfileMostProbablePattern(Dna[i],k, Profile)
        result.append(probablekmer)
    return result

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    t = len(Dna)
    length = len(Dna[0])
    beginpos = 0
    result = []
    for i in range(0,t):
        if k == length:
            beginpos = 0
        else:
            beginpos = random.randint(1, length - k + 1)
        result.append(Dna[i][beginpos-1:beginpos+k-1])
    return result

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs
    return  BestMotifs

### DO NOT MODIFY THE CODE BELOW THIS LINE ###
def RepeatedRandomizedMotifSearch(Dna, k, t):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(1000):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

