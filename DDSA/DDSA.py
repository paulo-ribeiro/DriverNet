# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import random
import math
QUANT_DRIVERS = 50





# FUNÇÕES UTIILZADAS PELOS ALGORITMOS DO DDSA
def patient_mutation_exist(gene, patient):
  return (mutation_matrix.at[patient,gene] == 1)

def verify_interaction(gi, gj):
  return influence_graph.has_edge(gi, gj)

def calc_coverage(neighbors_lists):
  all_neighbors_set = set()
  for neighbors_list in neighbors_lists:
    all_neighbors_set.update(neighbors_list)
  return len(all_neighbors_set)

def get_mutation_matrix_from_maf(maf_file_name):
  maf = pd.read_csv(maf_file_name, sep="\t", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
  mut = pd.crosstab(maf.Tumor_Sample_Barcode, maf.Hugo_Symbol).clip(upper=1)
  return mut





# LEITURA DOS ARQUIVOS E CONSTRUÇÃO DOS CONJUNTOS DE DADOS
mutation_matrix = get_mutation_matrix_from_maf("tcga_data_mutations_gbm.txt") # matriz de mutação do TCGA
influence_graph = nx.read_edgelist("Reactome_FIsInGene_2021.txt", delimiter='\t') # rede de influência do Reactome





# CONSTRUÇÃO DO GRAFO BIPARTIDO
bipartite_graph=nx.Graph()
patients = mutation_matrix.index.to_list() # Todos os dados da primeira coluna da matriz de mutação
genes = mutation_matrix.columns.to_list(); # Nomes das colunas da matriz de mutação 
# Adicionando os nós da partição da esquerda
left_partition = genes
bipartite_graph.add_nodes_from(left_partition, bipartite=0)
# Adicionando os nós da partição da direita
right_partition = []
for p in patients:
  for gene in genes:
    right_partition.append(gene + "_" + p)
bipartite_graph.add_nodes_from(right_partition, bipartite=1)
# Criação das arestas do grafo bipartido
green_nodes = set()
for gi in left_partition:
  for pk in patients:
    if patient_mutation_exist(gi, pk):
      green_nodes.add(gi)
      for gj in genes:
        if verify_interaction(gi, gj):
          label_node_right_partition = gj + "_" + pk
          bipartite_graph.add_edge(gi, label_node_right_partition)





# ALGORITMO SIMULATED ANNEALING
# Criando o dicionário com a lista de vizinhos de cada nó da partição esquerda que possuem vizinhos (green_nodes)
neighbors_list = []
for g in green_nodes: 
  neighbors_list.append(list(bipartite_graph.neighbors(g)))
neighbors_dictionary = dict(zip(green_nodes, neighbors_list))
# Definindo os parâmetros do algoritmo simulated annealing
initial_temperature = 10**3 # temperatura inicial 
cooling_factor = 1 - (10**-2)  # coeficiente de resfriamento
quant_iterations = 10**5
# Criação da lista dos genes da partição esquerda com grau > 0
left_partition_degree = list(green_nodes)
# Seleção aleatória dos possíveis genes drivers
best_set = []
random_number_list = random.sample(range(0, len(left_partition_degree)-1), QUANT_DRIVERS)
for i in range(0, QUANT_DRIVERS):
  g = left_partition_degree[random_number_list[i]]
  best_set.append(g)
# Criação da lista dos vizinhos dos genes existentes em "best_set"
neighbors_lists_best_set = []
for g in best_set:
  neighbors_lists_best_set.append(list(bipartite_graph.neighbors(g)))
best_coverage = calc_coverage(neighbors_lists_best_set)
current_temperature = initial_temperature
for i in range(0, quant_iterations):
  current_set = best_set.copy()
  neighbors_lists_current_set = neighbors_lists_best_set.copy()
  # Trocando os genes em "current_set"
  g_j = left_partition_degree[random.randint(0, len(left_partition_degree)-1)]
  if (not (g_j in best_set)):
    index = random.randint(0,len(current_set)-1)
    del current_set[index]
    del neighbors_lists_current_set[index]
    current_set.append(g_j)
    neighbors_lists_current_set.append(neighbors_dictionary.get(g_j))
  # Comparando "current_set" e "best_set" e verificando possível substituição de "best_set" por "current_set"  
  if (current_set != best_set):
    current_coverage = calc_coverage(neighbors_lists_current_set)
    if (current_coverage > best_coverage):
      best_set = current_set
      best_coverage = current_coverage
      neighbors_lists_best_set = neighbors_lists_current_set
    else:
      diff = abs(current_coverage - best_coverage)
      if (random.random() < math.exp(-diff/current_temperature)):
        best_set = current_set
        best_coverage = current_coverage
        neighbors_lists_best_set = neighbors_lists_current_set
  current_temperature = current_temperature * cooling_factor
# Ordenação pelo grau do nó
drivers_list1 = []
for g in best_set:
  drivers_list1.append((g, bipartite_graph.degree[g]))
drivers_order_list1 = sorted(drivers_list1, key=lambda x: x[1], reverse=True)
drivers_order_list1 = list(map(lambda x: x[0], drivers_order_list1))
# Ordenação pelo grau do nó (removendo os nós da partição direita)
drivers_order_list2 = []
for g in drivers_order_list1:
  drivers_order_list2.append((g, bipartite_graph.degree[g]))
  # Removendo os nós da partição direita do gene de maior grau
  g_neighbors = list(bipartite_graph.neighbors(g))
  for neighbor in g_neighbors:
    bipartite_graph.remove_node(neighbor)
drivers_order_list2 = sorted(drivers_order_list2, key=lambda x: x[1], reverse=True)
ddsa_drivers_list = list(map(lambda x: x[0], drivers_order_list2))
print("DDSA - Lista dos " + str(QUANT_DRIVERS) + " possíveis genes drivers melhores classificados:")
print('   '.join(ddsa_drivers_list[0:QUANT_DRIVERS]))
