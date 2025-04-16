# Nome: Gabriel Pinheiro de Campos
# RA: 156315
# Turma: NB
# Programa para contar nucleotídeos em um arquivo FASTA
# Lê o nome do gene na primeira linha e a sequência de DNA na segunda linha

# Passo 1: Abrir e ler o arquivo FASTA
with open("teste1.fasta", "r") as arquivo_fasta:
    linhas = arquivo_fasta.readlines()

# Passo 2: Extrair o nome do gene e a sequência de DNA
nome_gene = linhas[0].strip()
if nome_gene.startswith('>'):
    nome_gene = nome_gene[1:].strip()  # Remove o '>' do início e espaços em branco

sequencia_dna = linhas[1].strip().upper()  # Normaliza para maiúsculas

# Passo 3: Contar os nucleotídeos
contador_A = 0
contador_C = 0
contador_G = 0
contador_T = 0

for nucleotideo in sequencia_dna:
    if nucleotideo == 'A':
        contador_A += 1
    elif nucleotideo == 'C':
        contador_C += 1
    elif nucleotideo == 'G':
        contador_G += 1
    elif nucleotideo == 'T':
        contador_T += 1

# Passo 4: Calcular o total
total = contador_A + contador_C + contador_G + contador_T

# Passo 5: Exibir os resultados
print(f"Nome do gene: {nome_gene}")
print(f"A: {contador_A}")
print(f"C: {contador_C}")
print(f"G: {contador_G}")
print(f"T: {contador_T}")
print(f"Total: {total} nucleotideos")

# Salvar os resultados em um arquivo de saída
with open("resultado_contagem.txt", "w") as arquivo_saida:
    arquivo_saida.write(f"Nome do gene: {nome_gene}\n")
    arquivo_saida.write(f"A: {contador_A}\n")
    arquivo_saida.write(f"C: {contador_C}\n")
    arquivo_saida.write(f"G: {contador_G}\n")
    arquivo_saida.write(f"T: {contador_T}\n")
    arquivo_saida.write(f"Total: {total} nucleotideos\n")
