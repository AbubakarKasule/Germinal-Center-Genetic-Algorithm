"""
Author: Abubakar Kasule
Description: Driver code for my genetic algorithm
Note: Each GC starts with one epitope from the antigen. 
"""
# Imports
from numpy import record
from GerminalCenter import GerminalCenter
from Epitope import Epitope

# Constants
NUMBER_OF_GERMINAL_CENTERS = 15
ANTIBODY_POPULATION_SIZE = 20
ANTIBODY_MUTATION_RATE = 0.1
EPITOPE_MUTATION_RATE = 0.1
NUMBER_OF_GENERATIONS = 10

# variables
germinal_centers = []
epitopes = []

# Initialize epitopes and germinal centers
for i in range(NUMBER_OF_GERMINAL_CENTERS):
    ep = Epitope.create_random_epitope()
    epitopes.append(ep)
    germinal_centers.append(GerminalCenter(ANTIBODY_POPULATION_SIZE, ep, ANTIBODY_MUTATION_RATE, EPITOPE_MUTATION_RATE))

# Initialize antibody populations
for gc in germinal_centers:
    gc.initialize_antibody_population()

# Simulate evolution
curr_generation = 1
print(0)
while (curr_generation <= NUMBER_OF_GENERATIONS):
    print("-------------")
    germinal_centers[0].selection_and_reproduction()
    # for gc in germinal_centers:
    #     gc.selection_and_reproduction()

    curr_generation += 1
print(1)

with open('results.txt', 'w') as f:
    for i in range(len(germinal_centers)):
        gc = germinal_centers[i]
        print("Germinal Center #" + str(i + 1) + ": ", end="", file=f)
        print(gc.performance_history, file=f)


print("End of Code")