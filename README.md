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

### Atividade 4: Transcrição e Tradução de Sequências Genômicas

Implementação de algoritmos para a transcrição de DNA para RNA e a tradução de RNA para proteínas, aplicados ao genoma da *Escherichia coli* O157:H7 Sakai.

**Arquivos:**
- `transcricao_traducao.py`: Script Python para transcrição e tradução de sequências
- `transcricao_traducao.ipynb`: Notebook Jupyter com processos de transcrição e tradução detalhados
- `Ecoli_Sakai_cds_from_genomic.fna`: Arquivo FASTA com sequências codificantes da *E. coli* Sakai
- `Sakai_RNA.fasta`: Sequências de RNA após processo de transcrição
- `Frame1.fasta` a `Frame6.fasta`: Sequências traduzidas para proteínas em cada um dos seis frames de leitura
- `Tabela do código genético.png`: Imagem da tabela de código genético utilizada na tradução
- `Figura 1.png`, `Figura 2.png`, `Figura 3a.png`, `Figura 3b.png`: Imagens referentes ao processo de transcrição e tradução

### Atividade 5: Anotação Genômica - Identificação de ORFs

Desenvolvimento de algoritmo para identificação de regiões potencialmente codificadoras de proteínas (ORFs - Open Reading Frames) no genoma da *Escherichia coli* O157:H7 Sakai, incluindo o genoma principal e os dois plasmídeos.

**Arquivos:**
- `protein_finder.py`: Script Python para identificação de ORFs nos seis frames de leitura
- `protein_finder.ipynb`: Notebook Jupyter com implementação detalhada e análise do processo de identificação de ORFs
- `Ecoli Sakai sequence.fasta`: Sequência do genoma principal da *E. coli* Sakai
- `Ecoli Sakai Plasmid 1 sequence.fasta`: Sequência do plasmídeo pO157 da *E. coli* Sakai
- `Ecoli Sakai plasmid 2 sequence.fasta`: Sequência do plasmídeo pOSAK1 da *E. coli* Sakai
- `gc_frame1_ativ5_156315.fasta` a `gc_frame6_ativ5_156315.fasta`: Potenciais proteínas do genoma principal nos 6 frames
- `p1_frame1_ativ5_156315.fasta` a `p1_frame6_ativ5_156315.fasta`: Potenciais proteínas do plasmídeo pO157 nos 6 frames
- `p2_frame1_ativ5_156315.fasta` a `p2_frame6_ativ5_156315.fasta`: Potenciais proteínas do plasmídeo pOSAK1 nos 6 frames
- `README.md`: Documento detalhando a relação da atividade com o processo de anotação genômica

### Atividade 7: Montagem de Genoma com Algoritmo Guloso

Implementação de um algoritmo guloso para montagem de contigs a partir de reads curtas de DNA, simulando o processo de montagem de genomas.

**Arquivos:**
- `montagem_genoma.py`: Script Python implementando o algoritmo guloso para montagem de genoma
- `explicacao_algoritmo_guloso.md`: Documentação detalhada sobre o funcionamento do algoritmo guloso
- `reads4.fasta`: Arquivo FASTA com 4 reads curtas de DNA para montagem
- `contig.fasta`: Resultado da montagem - contig final obtido pelo algoritmo guloso
- `montagem_genoma_notebook.ipynb`: Notebook Jupyter com análises interativas


## Sobre a Disciplina

A Unidade Curricular de Algoritmos em Bioinformática aborda o desenvolvimento e aplicação de algoritmos computacionais para resolver problemas biológicos, com foco na análise de sequências de DNA, RNA e proteínas. As atividades práticas visam desenvolver habilidades na manipulação e análise de dados biológicos usando linguagens de programação como Python.

## Tecnologias Utilizadas

- Python 3
- BioPython
- Matplotlib e pandas para visualização de dados
- Manipulação de arquivos FASTA
- Análise de sequências genômicas
- Cálculo de propriedades biofísicas (temperatura de melting, conteúdo GC)
- Algoritmos gulosos para montagem de genoma

---

**Autor:** Gabriel Pinheiro de Campos  
**Instituição:** Universidade Federal de São Paulo (UNIFESP)
