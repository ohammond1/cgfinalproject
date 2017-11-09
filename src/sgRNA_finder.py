import pre_process

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
