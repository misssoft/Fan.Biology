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


