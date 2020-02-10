#In [1]: 

import numpy as np
import networkx as nx
import math
import csv

#In [2]: 

def compute_trophic_levels(DG):
adj_mtrx = nx.to_numpy_matrix(DG)
Kv = np.array([max(DG.in_degree(node), 1) for node in DG.nodes()])
A = np.diag(Kv) - np.transpose(adj_mtrx)
S = None
if np.linalg.det(A) == 0:
print('Error: Singular matrix')
return None, None
else:
S = np.linalg.solve(A, Kv)
trophic_levels = []
idx = 0
for node in DG.nodes():
trophic_levels.append([node, S[idx]])
DG.nodes()[node]['tr_lvl'] = float(S[idx])
idx += 1
return S, DG
def trophic_diff(DG):
nodes = DG.nodes(data=True)
trophic_diff = []
for edge in DG.edges():
t_diff = abs(nodes[edge[0]]['tr_lvl'] - nodes[edge[1]]['tr_lvl'])
trophic_diff.append([edge, t_diff])
return trophic_diff
def coherence_parameter(DG, trophic_diff):
q = 0
for x in trophic_diff:
q += x[1] * x[1]
q /= float(DG.number_of_edges())
q = math.sqrt(q - 1)
return q

#In [3]: 

def basal_ensemble_expectation(DG):
# number of edges connected to basal nodes
L_b = 0
basal_nodes, _ = get_basal_and_top_nodes(DG)
for node in basal_nodes:
L_b += len(DG.out_edges(node))
q_b = math.sqrt(DG.number_of_edges() / float(L_b) - 1)
return q_b
def get_basal_and_top_nodes(G):
basal_nodes = []
for node in G.nodes():
if G.in_degree(node) == 0:
basal_nodes.append(node)
top_nodes = []
for node in G.nodes():
if G.out_degree(node) == 0:
top_nodes.append(node)
return basal_nodes, top_nodes

#In [4]: 

def load_file(filename):
G = nx.DiGraph()
with open(filename, 'r') as csvfile:
file = csv.reader(csvfile, delimiter=',', quotechar='|')
for row in file:
if len(row) > 1 and row[1] != "" and row[1] != "Source":
G.add_edge(row[1], row[2])#, weight = row[3])
nx.draw_shell(G, with_labels=True, node_color="r", alpha = 0.7, arrowsize=20,
node_size=500, font_weight="bold")
return G
def add_trophic_levels(G):
_, G = compute_trophic_levels(G)
labels = {}
for node in G.nodes(data=True):
labels[node[0]] = round(node[1]['tr_lvl'], 1)
nx.draw_shell(G, labels = labels, with_labels=True, node_color="r", alpha = 0.7, arrowsize=20,
node_size=500, font_weight="bold")
return G
def print_coherence(G):
coherence = coherence_parameter(G, trophic_diff(G))
print("Coherence of the network (q): ", coherence)
basal_ensemble = basal_ensemble_expectation(G)
print("Basal ensemble (q'): ", basal_ensemble)
print("Normalized coherence (q/q'):", coherence / basal_ensemble)
