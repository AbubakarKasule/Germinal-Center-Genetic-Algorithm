"""
Author: Abubakar Kasule
Description: Utility functions for GCGA 
Note: 
"""

# Imports
import numpy
import math
import random
import networkx as nx

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


# Copied from https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3
def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5):

    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.  
    Licensed under Creative Commons Attribution-Share Alike 
    
    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.
    
    G: the graph (must be a tree)
    
    root: the root node of current branch 
    - if the tree is directed and this is not given, 
      the root will be found and used
    - if the tree is directed and this is given, then 
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given, 
      then a random choice will be used.
    
    width: horizontal space allocated for this branch - avoids overlap with other branches
    
    vert_gap: gap between levels of hierarchy
    
    vert_loc: vertical location of root
    
    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, pos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''
    
        if pos is None:
            pos = {root:(xcenter,vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            dx = width/len(children) 
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap, 
                                    vert_loc = vert_loc-vert_gap, xcenter=nextx,
                                    pos=pos, parent = root)
        return pos

            
    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)


