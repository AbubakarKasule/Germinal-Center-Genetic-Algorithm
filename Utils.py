"""
Author: Abubakar Kasule
Description: Utility functions for GCGA 
Note: 
"""

# Imports
import numpy
import math
import random

# Constants
base_to_bitpair = dict()
bitpair_to_base = dict()
bases = ["A", "T", "C", "G"]

base_to_bitpair["A"] = "00"
base_to_bitpair["T"] = "10"
base_to_bitpair["C"] = "01"
base_to_bitpair["G"] = "11"

bitpair_to_base["00"] = "A"
bitpair_to_base["10"] = "T"
bitpair_to_base["01"] = "C"
bitpair_to_base["11"] = "G"

# Functions

# Convert DNA sequence into Bitstring
def dna_to_bitstring(dna):
        bits = ""

        for basepair in dna:
            bits += base_to_bitpair[basepair]

        return bits

# Convert Bitstring sequence into DNA
def bitstring_to_dna(bitstring):
        dna = ""

        for i in range(0, len(bitstring) - 1, 2):
            dna += bitpair_to_base[bitstring[i] + bitstring[i + 1]]

        return dna

# Calculate Hamming Distance. 
def get_hamming_distance(bitstr1, bitstr2):
    count = 0

    n = len(bitstr1)
    m = len(bitstr2)

    if m != n:
        return math.inf # No hamming distance. String do not match


    for i in range(n):
        if bitstr1[i] != bitstr2[i]:
            count += 1

    return count

# Return base that is not the same as the given base
def return_different_base(excluded_base):
    l = list(bases)
    l.remove(excluded_base)

    return random.choice(l)

# Generate random dna sequence
def generate_random_dna(length):
    res = ""

    for i in range(length):
        res += random.choice(bases)

    return res




