# Computation Genomics Final project
# Oliver, Allie, Sanat, Jack
#
# pre_process.py 
#
# Program to process a reference genome and stores
# the resulting dictionary into a pickle dump


import pickle
import sys

def main():
    kmer_length = int(sys.argv[1])
    reference_file = sys.argv[2]
    with open(reference_file, "r") as infile:
        infile.readline() # skip the ID line
        sequence = infile.read()
        sequence = sequence.strip().replace('\n','')
        index = create_dict(sequence, kmer_length)

    pickle_file = reference_file[:len(reference_file)-3] + sys.argv[1] + ".p"
    pickle.dump(index, open(pickle_file, "wb"))
    


def create_dict(t, k):
    ''' Create index from all substrings of size 'length' '''
    index = {}
    for i in range(len(t) - k + 1):  # for each k-mer
        kmer = t[i:i+k]
        if kmer not in index:
            index[kmer] = [i]
        else:
            index[kmer].append(i)
        # could also have used collections.defaultdict
    return index


if __name__ == "__main__":
    main()
