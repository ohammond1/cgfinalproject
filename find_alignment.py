#!/usr/bin/env python3
import sys
import numpy as np
import pickle

def query(p, kmer_num, index, k):
    ''' Return index hits for first k-mer of P '''
    kmer = p[kmer_num:(kmer_num+k)]  # query with k-mer number kmer_num
    return index.get(kmer, [])

def edDistDp(x, y, worstDist):
    """ Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return distance. """
    D = np.zeros((len(x)+1, len(y)+1), dtype=int)
    D[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)

    minIndex = 0
    for col in range(len(y)+1):
        if D[len(x), col] < D[len(x), minIndex]:
            minIndex = col
    if D[len(x), minIndex] > worstDist:
        return (-1, -1)
    row = len(x)
    col = minIndex
    while row > 0:
        minNext = min(D[row-1, col-1], D[row-1, col], D[row, col-1])
        if minNext == D[row-1, col-1]:
            row -= 1
            col -= 1
        elif minNext == D[row-1, col]:
            row -= 1
        else:
            col -= 1
    startIndex = col
    return (startIndex, D[len(x), minIndex])


def queryIndex(p, t, index, k, worstDist):
    # Keep track of index of best matches
    bestIndex = -1
    # Keep best edit distance
    bestDistance = -1
    for i in range(0, len(p), k):
        matches = query(p, i, index, k)
        print (matches)
        for match in matches:
            startIndex = match - i - worstDist
            endIndex = match - i + len(p) + worstDist
            if startIndex < 0:
                startIndex = 0
            if endIndex > len(t):
                endIndex = len(t)
            ed = edDistDp(p, t[startIndex:endIndex], worstDist)
            if ed[0] != -1:
                if ed[1] < bestDistance or bestDistance == -1:
                    bestDistance = ed[1]
                    bestIndex = startIndex + ed[0]
    return bestIndex

k = 20
editDist = 3
index = pickle.load(open("genomes/woollyMonkeyHepB20.p", "rb"))
t = ""
with open("genomes/woollyMonkeyHepB.fa", "r") as infile:
    infile.readline()
    t = infile.read()
    t = t.strip().replace('\n', '')
p = "GACCTTTTTGAAAGATCAATACATGCACCTTTACCCCGTT"
print(queryIndex(p, t, index, k, editDist))
