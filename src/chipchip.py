#!/usr/bin/env python3

import sys
from sgRNA_finder import find_sgRNA
from find_alignment import find_alignment


def main():
    filename = sys.argv[1]
    p = sys.argv[2].strip()
    t = ""
    with open(filename, "r") as infile:
        infile.readline()
        t = infile.read()
        t = t.strip().replace('\n', '')
    indices = find_alignment(p, t)
    if indices[0] == -1:
        print("No valid alignments found")
        return
    find_sgRNA(filename, indices[0], indices[1])


if __name__ == "__main__":
    main()
