import dimod

import savanna
import networkx as nx
import random

G = nx.grid_graph(dim=[3,3])
print(G.nodes())
G.remove_node((2,2))
G.remove_edge((0,2),(1,2))
G.remove_edge((2,0),(2,1))
print(G.edges())
#G = nx.grid_graph(dim=[2,2])
#G.remove_node((1,1))
#G.add_edge((0,1),(1,0))

h = {}
#J = {e: -1 for e in G.edges()} #random.randint(0,1) for e in G.edges()}
J = {e: 2*random.randint(0,1)-1 for e in G.edges()}
bqm = dimod.BinaryQuadraticModel.from_ising(h, J)

pos = {v: v for v in G.nodes()}
rotation_system = savanna.rotation_system_from_coordinates(G, pos)


#bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): -1})
print(bqm)
#rotation_system = {'a': ['b'], 'b': ['a']}

sample = savanna.ground_state_sample(bqm, rotation_system)
print(sample)
print(bqm.energy(sample))

response = dimod.ExactSolver().sample(bqm)
for ground_sample in response.samples(1):
    print(ground_sample)
    print(bqm.energy(ground_sample))
