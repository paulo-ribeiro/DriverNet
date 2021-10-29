import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
QUANT_DRIVERS = 50




#FUNÇÕES UTIILZADAS PELOS ALGORITMOS DO DRIVERNET

def existe_mutacao_paciente(gene, paciente):
  indexPaciente = pacientesMatrizExpressaoGenica.index(paciente)
  if (matrizMutacao.at[indexPaciente,gene] == 1):
    return True
  else:
    return False

def verifica_expressao_genica(gene, paciente):
  indexPaciente = pacientesMatrizExpressaoGenica.index(paciente)
  return matrizExpressaoGenica.at[indexPaciente,gene]

def verifica_interacao(gi, gj):
  indexGi = genesRedeInteracao.index(gi)
  return redeInteracaoGenes.at[indexGi,gj]

def calc_quant_mutacoes(gene):
  index_gene = matrizMutacao.columns.to_list().index(gene)
  mutacoes = matrizMutacao.iloc[:,index_gene].to_list() # Todos os dados da primeira coluna da matriz de mutação*
  quant_mutacoes = 0
  for m in mutacoes:
    if (m == 1):
      quant_mutacoes += 1
  return quant_mutacoes





#LEITURA DOS ARQUIVOS CSV E CONSTRUÇÃO DOS CONJUNTOS DE DADOS

matrizMutacao = pd.read_csv('samplePatientMutationMatrix.csv')
matrizExpressaoGenica = pd.read_csv('samplePatientOutlierMatrix.csv')
redeInteracaoGenes = pd.read_csv('sampleInfluenceGraph.csv')
canonicalCancerDrivers = pd.read_csv('canonical_cancer_drivers.csv').iloc[:,0].to_list()

genes = matrizMutacao.columns.delete(0).to_list(); # Nomes das colunas da matriz de mutação (excluindo a 1ª coluna)
pacientes = matrizMutacao.iloc[:,0].to_list() # Todos os dados da primeira coluna da matriz de mutação
pacientesMatrizExpressaoGenica = matrizExpressaoGenica.iloc[:,0].to_list() # Todos os dados da primeira coluna da matriz
genesRedeInteracao = redeInteracaoGenes.iloc[:,0].to_list() # Todos os dados da primeira coluna da matriz






# DRIVERNET - ALGORITMO 1 

# CRIAÇÃO DOS NÓS DE CADA PARTIÇÃO DO GRAFO BIPARTIDO
BP=nx.Graph()
# Adicionando os nós da partição da esquerda
particaoEsquerda = genes
BP.add_nodes_from(particaoEsquerda, bipartite=0)
# Adicionando os nós da partição da direita
particaoDireita = []
for paciente in pacientes:
  for gene in genes:
    particaoDireita.append(gene + "_" + paciente)
BP.add_nodes_from(particaoDireita, bipartite=1)

#CRIAÇÃO DAS ARESTAS DO GRAFO BIPARTIDO
nosVerdes = []
nosVermelhos = []
for gi in particaoEsquerda:
  for pk in pacientes:
    if existe_mutacao_paciente(gi, pk):
      nosVerdes.append(gi)
      for gj in genes:
        if (verifica_expressao_genica(gj, pk)):
          if (verifica_interacao(gi, gj) == 1):
            labelNoParticaoDireita = gj + "_" + pk
            BP.add_edge(gi, labelNoParticaoDireita)
            nosVermelhos.append(labelNoParticaoDireita)






# DRIVERNET - ALGORITMO 2

BP_alg2 = BP.copy()
genes_drivers = []
while (BP_alg2.size() > 0):
  # Buscando gene de maior grau
  gene_maior_grau = ""
  maior_grau = 0
  for g in particaoEsquerda:
    grau_g = len(list(BP_alg2.neighbors(g)))
    if (grau_g > maior_grau):
      gene_maior_grau = g
      maior_grau = grau_g

  # Inserindo gene de maior grau na lista de genes drivers
  genes_drivers.append(gene_maior_grau)

  # Removendo os nós da partição direita do gene de maior grau
  vizinhos_g = list(BP_alg2.neighbors(gene_maior_grau))
  for vizinho in vizinhos_g:
    BP_alg2.remove_node(vizinho)

  #print("Gene: " + gene_maior_grau + " - Grau atual:" + str(maior_grau))


# Imprimindo lista de possíveis genes drivers
print("DriverNet - Lista dos " + str(QUANT_DRIVERS) + " possíveis genes drivers melhores classificados:")
print('   '.join(genes_drivers[0:QUANT_DRIVERS]))




