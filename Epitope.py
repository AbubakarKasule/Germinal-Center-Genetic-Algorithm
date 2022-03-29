"""
Author: Abubakar Kasule
Description: Class used to represent epitopes 
Note: 
"""

# Imports
import Utils
import random


# Epitope class
class Epitope:
    # Static constant
    MUTATION_RATE = 0.1
    DNA_SEQUENCE_LENGTH = 24

    def __init__(self, dna_sequence, parent, generation):
        self.generation = generation
        self.parent = parent
        self.dna_sequence = dna_sequence
        self.bit_string = Utils.dna_to_bitstring(dna_sequence)

    # Run on children after creation
    def mutate(self):
        new_dna = ""

        for base in self.dna_sequence:
            if random.random() < Epitope.MUTATION_RATE:
                new_dna += Utils.return_different_base(base)
            else:
                new_dna += base

        self.dna_sequence = new_dna
        self.bit_string = Utils.dna_to_bitstring(self.dna_sequence)

    # Function to handle reproduction and crossover
    @staticmethod
    def create_mutated_clone(parent, generation):

        child = Epitope(parent.dna_sequence, parent, generation)

        child.mutate()
        return child

    # Function to create random epitope
    @staticmethod
    def create_random_epitope():
        return Epitope(Utils.generate_random_dna(Epitope.DNA_SEQUENCE_LENGTH), None, 0)

      
