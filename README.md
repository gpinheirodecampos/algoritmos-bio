# Algoritmos em Bioinformática

Repositório criado para armazenar as atividades desenvolvidas na Unidade Curricular **Algoritmos em Bioinformática** do curso de graduação da Universidade Federal de São Paulo (UNIFESP).

## Estrutura do Repositório

O repositório está organizado por atividades práticas:

### Atividade 1: Contador de Nucleotídeos Simples

Implementação de um algoritmo básico para contar nucleotídeos (A, C, G, T) em sequências de DNA em formato FASTA.

**Arquivos:**
- `contador_nucleotideos_simples.py`: Script Python para contagem de nucleotídeos
- `teste1.fasta`: Arquivo de exemplo com sequência de DNA para teste
- `algoritmo_contagem_nucleotideos.txt`: Documentação do algoritmo
- `resultado_contagem.txt`: Saída da contagem de nucleotídeos

### Atividade 2: Análise do Genoma de E. coli Sakai

Análise mais avançada do genoma da bactéria *Escherichia coli* O157:H7 Sakai, incluindo seu cromossomo principal e dois plasmídeos.

**Arquivos:**
- `nucleotide_counter.py`: Script Python melhorado utilizando Biopython para análise das sequências
- `Ecoli Sakai sequence.fasta`: Sequência cromossômica da *E. coli* Sakai
- `Ecoli Sakai Plasmid 1 sequence (1).fasta`: Sequência do plasmídeo 1 da *E. coli* Sakai
- `Ecoli Sakai plasmid 2 sequence (2).fasta`: Sequência do plasmídeo 2 da *E. coli* Sakai
- `Ecoli_Sakai_genome_counts.csv`: Resultados da contagem para o genoma principal
- `Ecoli_Sakai_plasmid1_counts.csv`: Resultados da contagem para o plasmídeo 1
- `Ecoli_Sakai_plasmid2_counts.csv`: Resultados da contagem para o plasmídeo 2

### Atividade 3: Temperatura de Melting e Conteúdo GC

Análise de sequências de DNA da *Escherichia coli* O157:H7 Sakai para calcular as temperaturas de melting (anelamento) e sua relação com o conteúdo GC.

**Arquivos:**
- `calculo_temperatura_melting.py`: Script Python para análise das sequências
- `calculo_temperatura_melting.ipynb`: Notebook Jupyter com análises interativas e visualizações
- `Algoritmos.txt`: Documentação do algoritmo (narrativo e pseudocódigo)
- `Ecoli_Sakai_cds_from_genomic.fna`: Arquivo FASTA com 5155 potenciais genes da *E. coli* Sakai
- `Dados das sequencias.csv`: Resultados da contagem de nucleotídeos
- `Conteudo_GC.csv`: Resultados da análise do conteúdo GC
- `Temperatura_x_GC.csv`: Resultados da temperatura de melting e conteúdo GC
- `GC_x_Temperatura.png`: Gráfico de dispersão mostrando a relação entre conteúdo GC e temperatura de melting
- `BioMath_Linear Functions Applications.pdf`: Documento com informações sobre o cálculo da temperatura de melting
- `Figura 1.png`: Imagem de referência para curvas de melting de DNA

## Sobre a Disciplina

A Unidade Curricular de Algoritmos em Bioinformática aborda o desenvolvimento e aplicação de algoritmos computacionais para resolver problemas biológicos, com foco na análise de sequências de DNA, RNA e proteínas. As atividades práticas visam desenvolver habilidades na manipulação e análise de dados biológicos usando linguagens de programação como Python.

## Tecnologias Utilizadas

- Python 3
- BioPython
- Matplotlib e pandas para visualização de dados
- Manipulação de arquivos FASTA
- Análise de sequências genômicas
- Cálculo de propriedades biofísicas (temperatura de melting, conteúdo GC)

---

**Autor:** Gabriel Pinheiro de Campos  
**Instituição:** Universidade Federal de São Paulo (UNIFESP)
