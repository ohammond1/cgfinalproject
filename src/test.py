#!/usr/bin/env python3

import sgRNA_finder

seq_file = 'woollyMonkeyHepB.fa'

sgRNA = sgRNA_finder.sgRNAFinder(seq_file)
front_can = sgRNA.get_sgRNA_front(2340) # finding front sgRNA possibilities
back_can = sgRNA.get_sgRNA_back(2341) # finding back sgRNA possibilities

print('BWT TEST:')
for i in front_can:
    out = ''
    out += i + ' '
    out += str(sgRNA.query_index_bwt(i))
    print(out)

print('\nSELF COMPLEMENT SCORE TEST:')
for i in front_can:
    out = ''
    out += i + ' '
    out += str(sgRNA.self_complement_score(i))
    print(out)
print("***")
print(sgRNA.self_complement_score('AAAAGGGGGGGGTTTT'))
print("***")
# For testing
'''
print(front_can)
print(len(front_can))
print()
print(back_can)
print(len(back_can))
print()'''

# Dealing with first sgRNA
final_front_candidates = []
for can in front_can: # for each candidate
    hits = set() # kmer hits
    off_target_hits = [] # approximate matches
    # for each kmer, find matches and add appropriate offset to index
    for hit in sgRNA.queryIndex(can[:5])[0]:
        hits.add(hit)
    for hit in sgRNA.queryIndex(can[5:10])[0]:
        hits.add(hit-5)
    for hit in sgRNA.queryIndex(can[10:15])[0]:
        hits.add(hit-10)
    for hit in sgRNA.queryIndex(can[15:])[0]:
        hits.add(hit-15)

    # checking that the rest of the pattern matches
    for i in hits:
        if i >= 0 and i+20 < len(sgRNA.ref_genome): # checking that the alignment is in bounds
            ham = sgRNA.naive_approx_hamming(can, sgRNA.ref_genome[i:i+20], 3)[0]
            if ham != -1: # checking that hamming distance is <= 3
                off_target_hits.append(i)

    # append candidate and number of off target hits to final list
    # -1 because it will always find itself
    final_front_candidates.append((can, len(off_target_hits)-1))

print() # formatting

# printing out final list of candidates
print("Possible sgRNAs for 5' end: ")
for can in final_front_candidates:
    print(str(can[0]) + " " + str(can[1]))

print() #formatting

# Dealing with second sgRNA
final_back_candidates = []
for can in back_can: # for each candidate
    hits = set() # kmer hits
    off_target_hits = [] # approximate matches
    # for each kmer, find matches and add appropriate offset to index
    for hit in sgRNA.queryIndex(can[:5])[0]:
        hits.add(hit)
    for hit in sgRNA.queryIndex(can[5:10])[0]:
        hits.add(hit-5)
    for hit in sgRNA.queryIndex(can[10:15])[0]:
        hits.add(hit-10)
    for hit in sgRNA.queryIndex(can[15:])[0]:
        hits.add(hit-15)

    # checking that the rest of the pattern matches
    for i in hits:
        if i >= 0 and i+20 < len(sgRNA.ref_genome): # checking that the alignment is in bounds
            ham = sgRNA.naive_approx_hamming(can, sgRNA.ref_genome[i:i+20], 3)[0]
            if ham != -1: # checking that hamming distance is <= 3
                off_target_hits.append(i)

    # append candidate and number of off target hits to final list
    # -1 because it will always find itself
    final_back_candidates.append((can, len(off_target_hits)-1))

# printing out final list of back candidiates
print("Possible sgRNAs for 3' end: ")
for can in final_back_candidates:
    print(str(can[0]) + " " + str(can[1]))

print()

backs = sgRNA.find_CG_composition(final_back_candidates)
for can in backs:
    print(str(can[0]) + " " + str(can[1]) + " " + str(can[2]))

b = sgRNA.find_G_pos20(backs)
for can in b:
    print(str(can[0]) + " " + str(can[1]) + " " + str(can[2]) + " " + can[3])
