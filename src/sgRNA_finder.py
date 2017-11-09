import pre_process
from find_alignment import find_alignment

class sgRNAFinder:
    def __init__(self, ref_genome_file):
        self.ref_genome = self.parse_fasta(ref_genome_file)
        self.kmer_index = pre_process.create_dict(self.ref_genome, 5)

    def parse_fasta(self, filename):
        data = open('../genomes/' + filename, 'r').read()
        data = data.strip()
        data = data.split("\n")
        data = data[1:]
        data = ''.join(data)
        return data

    def get_sgRNA_front(self, front_index, front_offset=-20, back_offset=0):
        front_candidates = []
        for i in range(front_index+front_offset+1, front_index+back_offset+1):
            if i < 0:
                continue
            if i+19 >= len(self.ref_genome):
                break
            front_candidates.append(self.ref_genome[i:i+20])
        return front_candidates

    def get_sgRNA_back(self, back_index, front_offset=-20, back_offset=0):
        back_candidates = []
        for i in range(back_index+front_offset+1, back_index+back_offset+1):
            if i+19 >= len(self.ref_genome):
                break
            if i < 0:
                continue
            back_candidates.append(self.ref_genome[i:i+20])
        return back_candidates

    def query(self, candidate):
        return self.kmer_index.get(candidate, [])

    def queryIndex(self, candidate):
        k = 5
        offsets = []
        index_hits = 0
        for i in self.query(candidate):
            index_hits += 1
            # verify that the rest of P matches
            if candidate[k:] == self.ref_genome[i+k:i+len(candidate)]: 
                offsets.append(i)
        return(offsets, index_hits)

    def naive_approx_hamming(self, candidate, t, maxDistance):
        # loop over alignments
        for i in range(len(t) - len(candidate) + 1):
            nmm = 0
            for j in range(len(candidate)): # loop over characters
                if t[i+j] != candidate[j]: # compare characters
                    nmm += 1 # mismatch
                    if nmm > maxDistance:
                        break # exceeded max hamming dist
            if nmm <= maxDistance:
                return (i, nmm) # approximate match with number mismatches
        return (-1, -1)

    def find_CG_composition(self, candidates):
        possibilities = []
        for can in candidates:
            num_CG = 0
            for char in can[0]:
                if char == "C" or char == "G":
                    num_CG += 1
            comp_CG = float(num_CG) / 20
            possibilities.append((can[0], can[1], comp_CG))
        return possibilities

    def find_G_pos20(self, candidates):
        possibilities = []
        is_G = ''
        for can in candidates:
            if can[0][19] == "G":
                is_G = 'Y'
            else:
                is_G = 'N'
            possibilities.append((can[0], can[1], can[2], is_G))
        return possibilities



def find_sgRNA(seq_file, start, end):

    sgRNA = sgRNAFinder(seq_file)
    front_can = sgRNA.get_sgRNA_front(2340) # finding front sgRNA possibilities
    back_can = sgRNA.get_sgRNA_back(2341) # finding back sgRNA possibilities

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
            if i >= 0 and i+20 < len(sgRNA.ref_genome): # checking that the alignmen t is in bounds
                ham = sgRNA.naive_approx_hamming(can, sgRNA.ref_genome[i:i+20], 3)[0 ]
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
            if i >= 0 and i+20 < len(sgRNA.ref_genome): # checking that the alignmen t is in bounds
                ham = sgRNA.naive_approx_hamming(can, sgRNA.ref_genome[i:i+20], 3)[0 ]
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
