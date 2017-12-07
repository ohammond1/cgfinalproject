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
