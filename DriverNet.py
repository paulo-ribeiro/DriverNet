# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
QUANT_DRIVERS = 50





# FUNÇÕES UTIILZADAS PELOS ALGORITMOS DO DRIVERNET
def patient_mutation_exist(gene, patient):
  return (mutation_matrix.at[patient,gene] == 1)

def verify_gene_expression(gene, patient):
  return gene_expression_matrix.at[patient,gene]

def verify_interaction(gi, gj):
  return (influence_graph.at[gi,gj] == 1)





# LEITURA DOS ARQUIVOS CSV E CONSTRUÇÃO DOS CONJUNTOS DE DADOS
mutation_matrix = pd.read_csv('samplePatientMutationMatrix.csv', index_col=0)
influence_graph = pd.read_csv('sampleInfluenceGraph.csv', index_col=0)
gene_expression_matrix = pd.read_csv('samplePatientOutlierMatrix.csv', index_col=0)





# CONSTRUÇÃO DO GRAFO BIPARTIDO
bipartite_graph=nx.Graph()
patients = mutation_matrix.index.to_list() # Todos os dados da primeira coluna da matriz de mutação
genes = mutation_matrix.columns.to_list(); # Nomes das colunas da matriz de mutação (excluindo a 1ª coluna)
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
          if (verify_gene_expression(gj, pk)):
            label_node_right_partition = gj + "_" + pk
            bipartite_graph.add_edge(gi, label_node_right_partition)





# DRIVERNET - ALGORITMO GULOSO
bipartite_graph_alg2 = bipartite_graph.copy()
genes_drivers = []
while (bipartite_graph_alg2.size() > 0):
  # Buscando gene de maior grau
  highest_degree_gene = ""
  highest_degree = 0
  for g in left_partition:
    degree_g = bipartite_graph_alg2.degree[g]
    if (degree_g > highest_degree):
      highest_degree_gene = g
      highest_degree = degree_g
  # Inserindo gene de maior grau na lista de genes drivers
  genes_drivers.append(highest_degree_gene)
  # Removendo os nós da partição direita do gene de maior grau
  g_neighbors = list(bipartite_graph_alg2.neighbors(highest_degree_gene))
  for neighbor in g_neighbors:
    bipartite_graph_alg2.remove_node(neighbor)
# Imprimindo lista de possíveis genes drivers
print("DriverNet - Lista dos " + str(QUANT_DRIVERS) + " possíveis genes drivers melhores classificados:")
print('   '.join(genes_drivers[0:QUANT_DRIVERS]))