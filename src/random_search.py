import random

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    result = {}
    counts = list(Probabilities.values())
    sumvalues = sum(counts)
    for k,v in Probabilities.items():
        key = k
        value = v / sumvalues
        result[key] = value
    return  result

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = '' # output variable
    RandomScale = {}
    start = 0
    for k,v in Probabilities.items():
        end = start + v;
        RandomScale[k] = [start,end]
        start = end
    randomnum = random.uniform(0, 1)
    for k,v in RandomScale.items():
        if randomnum >= v[0] and randomnum < v[1]:
            kmer = k
            break
    return kmer

# Insert your Pr(Text, Profile) function here from Motifs.py.
def Pr(Text, Profile):
    length = len(Text)
    result = 1
    for i in range(length):
        symbol = Text[i]
        probability = Profile[symbol][i]
        result = result * probability
    return result

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[Text[i:i + k]] = Pr(Text[i:i + k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

if __name__ == '__main__':
    input = {'A':0.45, 'C':0.63, 'G':0.09, 'T': 0.27,'N': 0.36}
    output = Normalize(input)
    print(output)

