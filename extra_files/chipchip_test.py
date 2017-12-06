import unittest
from sgRNA_finder import sgRNAFinder


class Chipchip_Tester(unittest.TestCase):
    def test_front_candidates(self):
        # seq_file = 'herpesvirusHG52.fna'
        # sgRNA = sgRNAFinder(seq_file)
        # start, end = 0, len(sgRNA.ref_genome)
        # candidates = sgRNA.get_sgRNA_front(start, end, 1)
        # #print(candidates[-3])
        # test = "CGGGGGAAAAGAGGCGGGGCGGG"
        assert True

    def test_candidates(self):
        seq_file = 'herpesvirusHG52.fna'
        sgRNA = sgRNAFinder(seq_file)
        #sgRNA.ref_genome = sgRNA.ref_genome[:80]
        start, end = 0, len(sgRNA.ref_genome)
        candidates = sgRNA.get_sgRNA_front(start, end, 1)
        for i in candidates:
            print(i)
        assert True






if __name__ == '__main__':
    unittest.main()