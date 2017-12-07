from find_alignment import create_dict
import bwt
import operator

class sgRNAFinder:
    # Initialize object to have reference genome and 5-mer index
    def __init__(self, ref_genome_file):
        self.ref_genome = self.parse_fasta(ref_genome_file)
        self.kmer_index = create_dict(self.ref_genome, 5)

    # Parse fasta file to get reference genome as string
    def parse_fasta(self, filename):
        data = open('../genomes/' + filename, 'r').read()
        data = data.strip()
        data = data.split("\n")
        data = data[1:]
        data = ''.join(data)
        return data

    # Finding potential sgRNAs on both strands
    def get_sgRNA_cand(self, front_index, back_index, search_type):
        '''
        Search Types:
        1. Knockout
        2. Edit
        3. Activation
        4. Interference
        '''
        f_candidates = [] # Candidates on forward strand
        r_candidates = [] # Candidates on reverse strand
        # Setting offsets based on search type
        if search_type == 1 or search_type == 2:
            front_offset = front_index
            back_offset = back_index - 19
        elif search_type == 3:
            front_offset = front_index - 500
            back_offset = front_index - 50
        elif search_type == 4:
            front_offset = front_index - 300
            back_offset = front_index + 50
        else:
            raise Exception('Invalid Search Type Error')
        if front_offset < 0:
            front_offset = 0

        # Finding candidates in correct region with PAM sequence
        for i in range(front_offset, back_offset):
            if self.ref_genome[i+21:i+23] == 'GG':
                f_candidates.append(self.ref_genome[i:i+23])
            if i < 3:
                continue
            if self.ref_genome[i-3:i-1] == 'CC':
                r_candidates.append(self.ref_genome[i-3:i+20])

        return (f_candidates, r_candidates)

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
            for char in range(0,20):
                if can[0][char] == "C" or can[0][char] == "G":
                    num_CG += 1
            comp_CG = float(num_CG) / 20
            possibilities.append((can[0], can[1], comp_CG))
        return possibilities

    # Find candidates with G in 20th position for sgRNAs on forward strand
    # Determining whether the sgRNA violates parameters
    def find_G_pos20(self, candidates):
        possibilities = []
        is_G = ''
        for can in candidates:
            if can[0][19] == "G":
                is_G = 'Y'
            else:
                is_G = 'N'
            warn = "OK"
            # Checking to see if sgRNA violates parameters
            if can[1] > 0 or can[2] > .8 or can[2] < .4 or is_G == 'N' or self.self_complement_score(can[0]) > 0:
                warn = "WARNING"
            possibilities.append((can[0], can[1], can[2], is_G, warn))
        return possibilities

    # Find candidates with C in 3rd position for sgRNAs on reverse strand
    # Determining whether the sgRNA violates parameters
    def find_C_pos1(self, candidates):
        possibilities = []
        is_C = ''
        for can in candidates:
            if can[0][3] == 'C':
                is_C = 'Y'
            else:
                is_C = 'N'
            warn = "OK"
            # Checking to see if sgRNA violates parameters
            if can[1] > 0 or can[2] > .8 or can[2] < .4 or is_C == 'N' or self.self_complement_score(can[0]) > 0:
                warn = "WARNING"
            possibilities.append((can[0], can[1], can[2], is_C, warn))
        return possibilities

    # Finds off target hits using BWT algorithm
    def query_index_bwt(self, candidate):
        bwt_data = bwt.make_all(self.ref_genome)
        indices = bwt.find(candidate, self.ref_genome, mismatches=3, bwt_data=bwt_data)
        return (indices, -1)

    # Determining if the sgRNA is capable on folding in on itself
    # by seeing if it is self-complementary
    def self_complement_score(self, candidate):
        # n_map is for finding RNA complements
        n_map = {'A': 'U', 'U':'A', 'G':'C', 'C':'G'}
        first_half = self.get_RNA_complement(candidate[:int(len(candidate)/2)])
        second_half = self.get_RNA_complement(candidate[int(len(candidate)/2):])
        second_half_complement = ''
        for i in range(len(second_half)):
            second_half_complement += n_map[second_half[i]] #find complements
        second_half_complement = second_half_complement[::-1] # reverse string
        count = 0
        for i in range(0, len(second_half_complement)-3):
            if ''.join(first_half[i:i+4]) in ''.join(second_half_complement):
                count += 1
        return count

    def get_RNA_complement(self, dna):
        complement = {'A':'U', 'T':'A', 'G':'C', 'C':'G'}
        rna = [complement[i] for i in dna]
        return ''.join(rna)

    def get_reverse_complement(self, dna_original):
        dna = dna_original[:]
        dna = dna[::-1]
        rev_comp = ""

        for letter in dna:
            if letter == 'A':
                rev_comp += 'T'
            elif letter == 'T':
                rev_comp += 'A'
            elif letter == 'C':
                rev_comp += 'G'
            elif letter == 'G':
                rev_comp += 'C'
        return rev_comp


def find_sgRNA(seq_file, start, end, search_type):
    '''
    search_type:
    1 = Knockout
    2 = Edit
    3 = Activation
    4 = Interference
    '''

    '''
    Rows:
    1. sgRNA strand
    2. offset hits
    3. CG composition
    4. G 20th position (for rna)
    5. self complementarity (for rna)
    '''
    sgRNA = sgRNAFinder(seq_file)
    (front_can, back_can) = sgRNA.get_sgRNA_cand(start,end,search_type) # finding sgRNAs

    # Dealing with forward strand sgRNAs
    final_front_candidates = []
    for can in front_can: # for each candidate
        pam = can[20:23]
        can = can[0:20]
        f_hits = set() # kmer hits
        r_hits = set() # kmer hits
        off_target_hits = [] # approximate matches
        # for each kmer, find matches and add appropriate offset to index
        for hit in sgRNA.queryIndex(can[:5])[0]:
            f_hits.add(hit)
        for hit in sgRNA.queryIndex(can[5:10])[0]:
            f_hits.add(hit-5)
        for hit in sgRNA.queryIndex(can[10:15])[0]:
            f_hits.add(hit-10)
        for hit in sgRNA.queryIndex(can[15:])[0]:
            f_hits.add(hit-15)

        # checking that the rest of the pattern matches
        for i in f_hits:
            if i >= 0 and i+20 < len(sgRNA.ref_genome): # checking that the alignmen t is in bounds
                ham = sgRNA.naive_approx_hamming(can, sgRNA.ref_genome[i:i+20], 3)[0 ]
                if ham != -1: # checking that hamming distance is <= 3
                    off_target_hits.append(i)

        # Search for each kmer hit of reverse complement
        rev_can = sgRNA.get_reverse_complement(can)
        for hit in sgRNA.queryIndex(rev_can[:5])[0]:
            r_hits.add(hit)
        for hit in sgRNA.queryIndex(rev_can[5:10])[0]:
            r_hits.add(hit-5)
        for hit in sgRNA.queryIndex(rev_can[10:15])[0]:
            r_hits.add(hit-10)
        for hit in sgRNA.queryIndex(rev_can[15:20])[0]:
            r_hits.add(hit-15)

        # checking that the rest of the pattern matches
        for i in r_hits:
            # checking that the alignment is in bounds
            if i >= 0 and i+20 < len(sgRNA.ref_genome):
                ham = sgRNA.naive_approx_hamming(rev_can, sgRNA.ref_genome[i:i+20], 3)[0 ]
                if ham != -1: # checking that hamming distance is <= 3
                    off_target_hits.append(i)

        # append candidate and number of off target hits to final list
        # -1 because it will always find itself
        final_front_candidates.append((can + pam, len(off_target_hits)-1))

    # Dealing with second sgRNA
    final_back_candidates = []
    for can in back_can: # for each candidate
        pam = can[0:3]
        can = can[3:]
        f_hits = set() # kmer hits
        r_hits = set() # kmer hits
        off_target_hits = [] # approximate matches
        # for each kmer, find matches and add appropriate offset to index
        for hit in sgRNA.queryIndex(can[:5])[0]:
            f_hits.add(hit)
        for hit in sgRNA.queryIndex(can[5:10])[0]:
            f_hits.add(hit-5)
        for hit in sgRNA.queryIndex(can[10:15])[0]:
            f_hits.add(hit-10)
        for hit in sgRNA.queryIndex(can[15:20])[0]:
            f_hits.add(hit-15)

        # checking that the rest of the pattern matches
        for i in f_hits:
            if i >= 0 and i+20 < len(sgRNA.ref_genome): # checking that the alignmen t is in bounds
                ham = sgRNA.naive_approx_hamming(can, sgRNA.ref_genome[i:i+20], 3)[0 ]
                if ham != -1: # checking that hamming distance is <= 3
                    off_target_hits.append(i)


        # Search for each kmer hit of reverse complement
        rev_can = sgRNA.get_reverse_complement(can)
        for hit in sgRNA.queryIndex(rev_can[:5])[0]:
            r_hits.add(hit)
        for hit in sgRNA.queryIndex(rev_can[5:10])[0]:
            r_hits.add(hit-5)
        for hit in sgRNA.queryIndex(rev_can[10:15])[0]:
            r_hits.add(hit-10)
        for hit in sgRNA.queryIndex(rev_can[15:])[0]:
            r_hits.add(hit-15)

        # checking that the rest of the pattern matches
        for i in r_hits:
            # checking that the alignment is in bounds
            if i >= 0 and i+20 < len(sgRNA.ref_genome):
                ham = sgRNA.naive_approx_hamming(rev_can, sgRNA.ref_genome[i:i+20], 3)[0 ]
                if ham != -1: # checking that hamming distance is <= 3
                    off_target_hits.append(i)

        # append candidate and number of off target hits to final list
        # -1 because it will always find itself
        final_back_candidates.append((pam + can, len(off_target_hits)-1))

    # Adding information about parameters to sgRNAs
    backs = sgRNA.find_CG_composition(final_back_candidates)
    fronts = sgRNA.find_CG_composition(final_front_candidates)
    f = sgRNA.find_G_pos20(fronts)
    b = sgRNA.find_C_pos1(backs)

    f.sort(key=operator.itemgetter(1, 4)) #sorts list by off-target hits, warning
    b.sort(key=operator.itemgetter(1, 4)) #sorts list by off-target hits, warning
    print("FRONT sgRNA")
    # Format output
    out = 'Sequence\t\tOff-site hits\tCG%    G-20th\tSelf–Complementarity\tWarning\n'
    f_set = set()
    b_set = set()
    for count, seq in enumerate(f):
        if f[count] in f_set:
            continue
        f_set.add(f[count])
        out += f[count][0] + '\t' + str(f[count][1]) + ' \t\t' + str(f[count][2]) + '   \t' + f[count][3]
        out += '\t' + str(sgRNA.self_complement_score(f[count][0]))
        out += '\t\t\t' + f[count][4]
        out += '\n'
    print(out)
    print("-------------------------------------------------------------------------")
    print("BACK sgRNA")
    # Format output
    out = 'Sequence\t\tOff-site hits\tCG%    G-20th\tSelf–Complementarity\tWarning\n'
    for count, seq in enumerate(b):
        if b[count] in b_set:
            continue
        b_set.add(b[count])
        out += b[count][0] + '\t' + str(b[count][1]) + ' \t\t' + str(b[count][2]) + '   \t' + b[count][3]
        out += '\t' + str(sgRNA.self_complement_score(b[count][0]))
        out += '\t\t\t' + b[count][4]
        out += '\n'
    print(out)
