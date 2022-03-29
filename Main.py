"""
Author: Abubakar Kasule
Description: Driver code for my genetic algorithm
Note: Each GC starts with one epitope from the antigen. 
"""
# Imports
from GerminalCenter import GerminalCenter
from Epitope import Epitope
from Antibody import Antibody
import matplotlib.pyplot as plt

# Constants
NUMBER_OF_GERMINAL_CENTERS = 15
ANTIBODY_POPULATION_SIZE = 20
ANTIBODY_MUTATION_RATE = 0.3
EPITOPE_MUTATION_RATE = 0.1
NUMBER_OF_GENERATIONS = 5

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

print(0)

for gc in germinal_centers:
    curr_generation = 1
    while (curr_generation <= NUMBER_OF_GENERATIONS):
        Antibody.ID = 21  # reset id count for each gc
        gc.selection_and_reproduction()

        curr_generation += 1
print(1)

with open('results.txt', 'w') as f:
    for i in range(len(germinal_centers)):
        gc = germinal_centers[i]
        print("Germinal Center #" + str(i + 1) + ": ", end="", file=f)
        print(gc.performance_history[0], gc.performance_history[-1], file=f)


# Graphics stuff

# line graph fitness x generation for each gc
for i in range(len(germinal_centers)):
    gc = germinal_centers[i]
    gc.generate_fitness_to_generation_graph("generation_fitness_graphs/GC" + str(i + 1))

# lineage of best antibody demarcated by generation
for i in range(len(germinal_centers)):
    gc = germinal_centers[i]
    gc.generate_lineage_tree_for_best_antibody("lineage_trees/GC" + str(i + 1))


# bar graph comparing gc success

gc_ids = ["GC#" + str(x + 1) for x in range(NUMBER_OF_GERMINAL_CENTERS)]
gc_best_fitness = [x.get_fitness(x.return_best_antibody()) for x in germinal_centers]
plt.figure(3,figsize=(15, 10))
plt.bar(gc_ids, gc_best_fitness)

plt.xlabel("Germinal Center")
plt.ylabel("Fitness")
plt.title("Fitness Comparison between Germinal Centers")
plt.savefig("GC results")
plt.clf()


print("End of Code")