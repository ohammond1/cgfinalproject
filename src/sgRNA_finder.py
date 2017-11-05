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
