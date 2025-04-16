import os
from Bio import SeqIO

import functools
imprimir = functools.partial(print, flush=True)

arquivo_genoma = "Ecoli Sakai sequence.fasta"
arquivo_plasmideo1 = "Ecoli Sakai Plasmid 1 sequence (1).fasta"
arquivo_plasmideo2 = "Ecoli Sakai plasmid 2 sequence (2).fasta"

imprimir(f"Executando script a partir de: {os.getcwd()}")
imprimir(f"Arquivos no diretório atual: {os.listdir('.')}")

def contar_nucleotideos(seq):
    contagem_a = seq.count('A')
    contagem_c = seq.count('C')
    contagem_g = seq.count('G')
    contagem_t = seq.count('T')
    total = len(seq)
    return {'A': contagem_a, 'C': contagem_c, 'G': contagem_g, 'T': contagem_t, 'Total': total}

def processar_fasta(caminho_arquivo):
    with open(caminho_arquivo, "r") as manipulador:
        for registro in SeqIO.parse(manipulador, "fasta"):
            sequencia = registro.seq
            descricao = registro.description
            contagens = contar_nucleotideos(sequencia)
            return registro, contagens
    return None, None

def escrever_csv(nome_arquivo, registro, contagens):
    arquivo_saida = f"{nome_arquivo.split('.')[0]}_counts.csv"
    with open(arquivo_saida, "w") as f:
        f.write(f"{registro.description}\n")
        f.write("A,C,G,T,Total\n")
        f.write(f"{contagens['A']},{contagens['C']},{contagens['G']},{contagens['T']},{contagens['Total']}\n")
    return arquivo_saida

def principal():
    imprimir("Processando arquivos FASTA...")
    
    try:
        imprimir(f"Processando {arquivo_genoma}...")
        if os.path.exists(arquivo_genoma):
            imprimir(f"Arquivo encontrado: {arquivo_genoma}")
            registro_genoma, contagens_genoma = processar_fasta(arquivo_genoma)
            imprimir(f"Arquivo de genoma processado: {arquivo_genoma}, registro obtido: {registro_genoma is not None}")
            if registro_genoma:
                imprimir(f"Contagens do genoma: {contagens_genoma}")
                arquivo_saida = escrever_csv("Ecoli_Sakai_genome", registro_genoma, contagens_genoma)
                imprimir(f"Contagens escritas em {arquivo_saida}")
            else:
                imprimir(f"Falha ao processar arquivo do genoma: {arquivo_genoma}")
        else:
            imprimir(f"Arquivo do genoma não encontrado: {arquivo_genoma}")

        imprimir(f"Processando {arquivo_plasmideo1}...")
        if os.path.exists(arquivo_plasmideo1):
            registro_plasmideo1, contagens_plasmideo1 = processar_fasta(arquivo_plasmideo1)
            if registro_plasmideo1:
                arquivo_saida = escrever_csv("Ecoli_Sakai_plasmid1", registro_plasmideo1, contagens_plasmideo1)
                imprimir(f"Contagens escritas em {arquivo_saida}")
            else:
                imprimir(f"Falha ao processar arquivo do plasmídeo 1: {arquivo_plasmideo1}")
        else:
            imprimir(f"Arquivo do plasmídeo 1 não encontrado: {arquivo_plasmideo1}")
        
        imprimir(f"Processando {arquivo_plasmideo2}...")
        if os.path.exists(arquivo_plasmideo2):
            registro_plasmideo2, contagens_plasmideo2 = processar_fasta(arquivo_plasmideo2)
            if registro_plasmideo2:
                arquivo_saida = escrever_csv("Ecoli_Sakai_plasmid2", registro_plasmideo2, contagens_plasmideo2)
                imprimir(f"Contagens escritas em {arquivo_saida}")
            else:
                imprimir(f"Falha ao processar arquivo do plasmídeo 2: {arquivo_plasmideo2}")
        else:
            imprimir(f"Arquivo do plasmídeo 2 não encontrado: {arquivo_plasmideo2}")
            
    except Exception as e:
        imprimir(f"Ocorreu um erro: {str(e)}")
        with open("log_de_erro.txt", "w") as log:
            log.write(f"Erro: {str(e)}\n")

if __name__ == "__main__":
    principal()
