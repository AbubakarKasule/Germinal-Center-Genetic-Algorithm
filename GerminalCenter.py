"""
Author: Abubakar Kasule
Description: Class used to represent Germinal Center
Note: Each GC starts with one epitope from the antigen. 
"""

# Imports
import functools
import numpy
import Utils
import random
from Antibody import Antibody
from Epitope import Epitope

class GerminalCenter:
    def __init__(self, antibody_population_size, initial_epitope, antibody_muatation_rate, epitope_mutation_rate):
        self.antibody_population_size = antibody_population_size
        self.antibody_population = []   # Empty until population is initialized
        self.epitope_population = [initial_epitope]
        self.generation = 1
        self.performance_history = []
 
        # Set static mutation rates
        Antibody.MUTATION_RATE = antibody_muatation_rate
        Epitope.MUTATION_RATE = epitope_mutation_rate
    
    def initialize_antibody_population(self):
        self.add_naive_pop()
        self.add_partially_binding_pop()

    def add_naive_pop(self):
        for i in range(self.antibody_population_size//2):
            self.antibody_population.append(Antibody.create_random_antibody())

    def add_partially_binding_pop(self):
        for i in range(self.antibody_population_size//2):
            self.antibody_population.append(Antibody.create_partially_binding_antibody(self.epitope_population[0].dna_sequence, 0.2))  # Partially binds to 20% of epitope

    def get_fitness(self, antibody):
        fitness = Antibody.DNA_SEQUENCE_LENGTH * 2 * len(self.epitope_population)

        for epitope in self.epitope_population:
            fitness -= Utils.get_hamming_distance(antibody.bit_string, epitope.bit_string)

        return fitness / len(self.epitope_population)

    def fitness_comparator(self, item1, item2):
        if self.get_fitness(item1) < self.get_fitness(item2):
            return -1
        elif self.get_fitness(item1) > self.get_fitness(item2):
            return 1
        else:
            return 0

    def selection_and_reproduction(self):
        print(0, end='')
        # Sort by fitness
        sorted(self.antibody_population, key=functools.cmp_to_key(self.fitness_comparator))

        # get rid of low performing Antibodies
        self.antibody_population = list(self.antibody_population[len(self.antibody_population)//2:])

        # create new final population by breeding remaining antibodies
        n = len(self.antibody_population) - 1
        print(1, end='') 
        for i in range(0, n, 2):
            children = Antibody.create_children(self.antibody_population[i], self.antibody_population[i + 1], self.generation)
            self.antibody_population.append(children[0])
            self.antibody_population.append(children[1])

        if (len(self.antibody_population) != self.antibody_population_size):
            print("SOMETHING IS VERY WRONG: GC70")
            
        print(2, end='')
        # Randomly add new mutated epitope
        if random.randint(0, 10) == 3:         # 3 is the magic number
            for ep in self.epitope_population:
                self.epitope_population.append(Epitope.create_mutated_clone(ep, self.generation))

        self.generation += 1

        print(3, end='')
        # Record fitness of best antibody
        self.performance_history.append(self.get_fitness(self.return_best_antibody()))

        print(4)

    def return_best_antibody(self):
        # Sort by fitness
        sorted(self.antibody_population, key=functools.cmp_to_key(self.fitness_comparator))

        return self.antibody_population[-1]
