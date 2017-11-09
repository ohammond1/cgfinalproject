#!/usr/bin/env python3
import sys
import numpy as np
import pickle
import math
from pre_process import create_dict 


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

    split = partition(p, worstDist + 2)
    kmer_len = len(split[0])
    i = 0
    for kmer in split:
        #matches = query(kmer_len, kmer, index)
        matches = index.get(kmer, [])
        print(matches)
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
        i += k
    return bestIndex

def partition(p, pieces=2):
    assert len(p) >= pieces
    base, mod = int(math.ceil(len(p) / pieces)), len(p) % pieces
    idx = 0
    ps = []
    modAdjust = 1
    for i in range(0, pieces):
        if i == pieces-1:
            break
        else:
            newIdx = idx + base
            ps.append(p[idx:newIdx])
            idx = newIdx
    print(ps)
    return ps


def main():
    editDist = int(sys.argv[1])
    filename = sys.argv[2]
    p = sys.argv[3].strip()
    t = ""
    k = int(math.ceil(len(p) / (editDist + 2)))
    print(k)
    print(len(p))
    with open(filename, "r") as infile:
        infile.readline()
        t = infile.read()
        t = t.strip().replace('\n', '')
    index = create_dict(t,k)
    print(queryIndex(p, t, index, k, editDist))



if __name__ == "__main__":
    main()
