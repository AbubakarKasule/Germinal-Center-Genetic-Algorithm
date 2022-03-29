"""
Author: Abubakar Kasule
Description: Class used to represent antibodies 
Note: Realistically this should be an entire B-cell but we are abstracting to just the relevant parts AKA the antibody
"""

# Imports
import Utils
import random


# Antibody class
class Antibody:
    # Static constant
    MUTATION_RATE = 0.1    # One in Five will mutate by default
    DNA_SEQUENCE_LENGTH = 24
    
    # Static variable
    ID = 1

    def __init__(self, dna_sequence, parent1, parent2, generation):
        self.generation = generation
        self.parents = (parent1, parent2)
        self.dna_sequence = dna_sequence
        self.bit_string = Utils.dna_to_bitstring(dna_sequence)
        self.id = Antibody.ID

        Antibody.ID += 1

    def __str__(self) -> str:
        return str("N" + str(self.id) + "G" + str(self.generation))
        # return str("Node #" + str(self.id) + " Generation #" + str(self.generation))

    # Run on children after creation
    def mutate(self):
        new_dna = ""

        for base in self.dna_sequence:
            if random.random() < Antibody.MUTATION_RATE:
                new_dna += Utils.return_different_base(base)
            else:
                new_dna += base

        self.dna_sequence = new_dna
        self.bit_string = Utils.dna_to_bitstring(new_dna)

    # Function to handle reproduction and crossover
    @staticmethod
    def create_children(parent1, parent2, generation):
        
        crossover_point = random.randint(0, len(parent1.dna_sequence) - 1)

        # Two possible children
        child_dna_1 = parent1.dna_sequence[:crossover_point] + parent2.dna_sequence[crossover_point:]
        child_dna_2 = parent2.dna_sequence[:crossover_point] + parent1.dna_sequence[crossover_point:]

        child_1 = Antibody(child_dna_1, parent1, parent2, generation)
        child_2 = Antibody(child_dna_2, parent2, parent1, generation)

        child_1.mutate()
        child_2.mutate()

        if (len(child_1.dna_sequence) != Antibody.DNA_SEQUENCE_LENGTH or len(child_2.dna_sequence) != Antibody.DNA_SEQUENCE_LENGTH):
            print("SOMETHING IS VERY WRONG: AB57")

        return (child_1, child_2)

    # Function to create random antibody
    @staticmethod
    def create_random_antibody():
        return Antibody(Utils.generate_random_dna(Antibody.DNA_SEQUENCE_LENGTH), None, None, 0)

    # Function to partially binding antibody
    @staticmethod
    def create_partially_binding_antibody(epitope_dna, percentage_of_conserved_region):
        random_dna = Utils.generate_random_dna(Antibody.DNA_SEQUENCE_LENGTH)

        if random.randint(0, 1) == 0:
            # Conserved region is in the front
            cronserved_region_threshold = round(len(random_dna) * percentage_of_conserved_region)

            result = str(epitope_dna[:cronserved_region_threshold])

            curr_len = len(result)

            while(curr_len < Antibody.DNA_SEQUENCE_LENGTH):
                result += random_dna[curr_len]
                curr_len += 1

                if result == epitope_dna:
                    print("SOMETHING IS VERY WRONG: AB84")
        
        else:
            # Conserved region is in the back
            cronserved_region_threshold = len(random_dna) - round(len(random_dna) * percentage_of_conserved_region)

            result = "" 

            curr_len = len(result)

            while(curr_len < cronserved_region_threshold):
                result += random_dna[curr_len]
                curr_len += 1

            result += str(epitope_dna[cronserved_region_threshold:])

            if result == epitope_dna:
                print("SOMETHING IS VERY WRONG: AB101")

        if (len(result) != Antibody.DNA_SEQUENCE_LENGTH):
            print("SOMETHING IS VERY WRONG: AB104")

        

        return Antibody(result, None, None, 0)

    


