# cgfinalproject

Tool designed to find CRIPSR/Cas9 sgRNA for specified sequences.

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
