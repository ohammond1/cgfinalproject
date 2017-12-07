# cgfinalproject

Tool designed to find CRIPSR/Cas9 sgRNA for specified sequences.

# How to run
To run our code, execute:
```
./chipchip.py
```
and input parameters through the prompts.

Search Types:
1. Knockout
2. Edit
3. Activation
4. Interference

Genome filename is the path from the current working directory to the genome's
.fna file.

The sequence to search for should be all uppercase and consist of A, C, G, T.

# Our files

chipchip.py
Implements a menu to make user input easy for ChipChip.

find_alignment.py
Finds the start and end indices of the best approximate alignment of a sequence
within a genome.

sgRNA_finder.py
Finds the possible target sequences for sgRNA in a genome based on a sequence's
alignment in the genome, as well as the type of search.

timer_chipchip.py
Implements ChipChip with command line arguments as opposed to a menu. Is used
for timing the difference between BWT and kmer index as methods for finding off-
target hits.

timer.sh
Timer file.


# Running Tests

To run tests referenced in our paper, use the following commands

./chipchip.py
1
../genomes/herpesvirusHG52.fna
CCTCGGGTTCCCAAGACCTATCACGTGTGCGCAGGGGAGGGGAGGACGCGGGGGAGGGGAGGACGC
GGGGGAGGGGAGGACGCGGGGGATATATAAAGCGGTAGAAAGCGCGGG 


In order to get the results for CHOPCHOP as shown in the paper, go to the website 
http://chopchop.cbu.uib.no/
Once at this website change the search type from gene target to fasta target
by clicking on the Fasta target button above Find Target. Next enter the target
sequence above, and click Find Targets.


Timing tests
To run the timing tests described in the paper enter the following commands.
On the ugrad server enter:
time ./timer.sh timer_chipchip.py ../genomes/herpesvirusHG52.fna <target sequence>

where in the case of our testing target sequence was 100, 250, and 500 base
pair prefixes of the following sequence.
GCCGGGCGGGGGCGCGCGGCGGCCGGGCGGGGGCGCGCGGCGGCCGGGCGGGGGCGCGCTTTCCCCGCGTCGCCCCT
CGGGTTCCCAAGACCTATCACGTGTGCGCAGGGGAGGGGAGGACGCGGGGGAGGGGAGGACGCGGGGGAGGGGAGGA
CGCGGGGGATATATAAAGCGGTAGAAAGCGCGGGAATGGGCATATTGGACCCGCGTGATTCGGTTGCTCGCGGTTGT
CTTGTTTGGACGTTTTTTATGCGGGAACAAGGGGGCTTACCGGTTACACTGTCCGCTCGCTATGGGGTTCGTCTGTC
TGTTTGGGCTTGTCGTTATGGGAGCCTGGGGGGCGTGGGGTGGGTCACAGGCAACCGAATATGTTCTTCGTAGTGTT
ATTGCCAAAGAGGTGGGGGACATACTAAGAGTGCCTTGCATGCGGACCCCCGCGGACGATGTTTCTTGGCGCTACGA
GGCCCCGTCCGTTATTGACTATGCCCGCATAGACGGAA
