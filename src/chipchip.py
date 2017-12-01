#!/usr/bin/env python3

import sys, os
from sgRNA_finder import find_sgRNA
from find_alignment import find_alignment


def main():
    print("Search type:")
    print("1. Knockout")
    print("2. Edit")
    print("3. Activation")
    print("4. Interference")
    search = input("Enter search type [1-4]: ")
    while search != "1" and search != "2" and search != "3" and search != "4":
        search = input("Please enter a valid input [1-4]: ")
    filename = input("Enter genome filename: ")
    while not os.path.isfile(os.getcwd() + "/" + filename):
        filename = input("Please enter a valid genome filename: ")
    t = ""
    with open(filename, "r") as infile:
        infile.readline()
        t = infile.read()
        t = t.strip().replace('\n', '')
    p = input("Enter sequence to search for: ")
    p = p.strip()
    indices = find_alignment(p, t)
    if indices[0] == -1:
        print("No valid alignments found")
        return

    find_sgRNA(filename, indices[0], indices[1], int(search))


if __name__ == "__main__":
    main()
