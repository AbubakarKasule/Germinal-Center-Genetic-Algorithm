"""
Author: Abubakar Kasule
Description: Class used to represent Germinal Center
Note: Each GC starts with one epitope from the antigen. 
"""

# Imports
import functools
import Utils
import random
from Antibody import Antibody
from Epitope import Epitope
import matplotlib.pyplot as plt
import networkx as nx
#from networkx.drawing.nx_agraph import graphviz_layout

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
        Antibody.ID = 1
        self.add_naive_pop()
        #self.add_naive_pop()
        self.add_partially_binding_pop()

        

        

    def add_naive_pop(self):
        for i in range(self.antibody_population_size//2):
            self.antibody_population.append(Antibody.create_random_antibody())

    def add_partially_binding_pop(self):
        for i in range(self.antibody_population_size//2):
            self.antibody_population.append(Antibody.create_partially_binding_antibody(self.epitope_population[0].dna_sequence, 0.05))  # Partially binds to 5% of epitope guaranteed

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
        # print(0, end='')
        # Sort by fitness
        self.antibody_population = sorted(self.antibody_population, key=functools.cmp_to_key(self.fitness_comparator))

        # get rid of low performing Antibodies
        self.antibody_population = list(self.antibody_population[len(self.antibody_population)//2:])

        # create new final population by breeding remaining antibodies
        n = len(self.antibody_population) - 1
        
        for i in range(0, n, 2):
            children = Antibody.create_children(self.antibody_population[i], self.antibody_population[i + 1], self.generation)
            self.antibody_population.append(children[0])
            self.antibody_population.append(children[1])

        if (len(self.antibody_population) != self.antibody_population_size):
            print("SOMETHING IS VERY WRONG: GC70")
            
        
        # Randomly add new mutated epitope
        if random.randint(0, 1000) == 3:         # 3 is the magic number
            self.epitope_population.append(Epitope.create_mutated_clone(self.epitope_population[0], self.generation))
            # for ep in list(self.epitope_population):
            #     self.epitope_population.append(Epitope.create_mutated_clone(ep, self.generation))

        self.generation += 1

        
        # Record fitness of best antibody
        self.performance_history.append(self.get_fitness(self.return_best_antibody()))



    def return_best_antibody(self):
        # Sort by fitness
        self.antibody_population = sorted(self.antibody_population, key=functools.cmp_to_key(self.fitness_comparator))
        with open('pop.txt', 'w') as f:
            print([self.get_fitness(x) for x in self.antibody_population], file=f)

        return self.antibody_population[-1]

    # Graphics functions
    def generate_fitness_to_generation_graph(self, filename):
        gens = [int(x + 1) for x in range(len(self.performance_history))]

        plt.xlabel("Generation")
        plt.ylabel("Fitness")
        plt.title("Fitness Comparison between Generations")
        #plt.xticks(gens, gens[::2])

        plt.gca().xaxis.set_major_locator(plt.MultipleLocator(100))

        plt.plot(gens, self.performance_history)
        plt.savefig(filename)
        plt.clf()

    def generate_lineage_tree_for_best_antibody(self, filename):
        g = nx.DiGraph()

        nodes = [self.return_best_antibody()]

        curr = 0

        # Add nodes until out of range error
        try:
            while True:
                parents = nodes[curr].parents

                if parents[0] is not None and parents[1] is not None:
                    if parents[0] not in nodes:
                        nodes.append(parents[0])

                    if parents[1] not in nodes:
                        nodes.append(parents[1])

                curr += 1
        except:
            print("No more nodes to add")

        # Add nodes to graph/tree
        for node in nodes:
            parents = node.parents
            
            for prnt in list(parents):
                if prnt is not None:
                    if str(str(prnt) + "\n") in g:
                        i = 1

                        while str(str(prnt) + "-" + str(i) + "\n") in g:
                            i += 1

                        g.add_edge(str(node) + "\n", str(str(prnt) + "-" + str(i) + "\n"))
                    
                    else:
                        g.add_edge(str(node) + "\n", str(prnt) + "\n")

        if (len(nodes) == 0):
            print("yyyyy")

        plt.figure(3,figsize=(12, 5))
        plt.title("Lineage tree for best performing Antibody")
        

        pos=Utils.hierarchy_pos(g, str(nodes[0]) + "\n", 100)
        nx.draw(g, pos=pos, node_size=50, with_labels=True, arrows=False)

        
        plt.savefig(filename)
        plt.clf()

        


